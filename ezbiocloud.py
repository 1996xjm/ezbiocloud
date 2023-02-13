import http.cookiejar, urllib.request, urllib.parse
import json
import math
from urllib.error import HTTPError
import time
from Bio import SeqIO
import os
from os import path
import re
import optparse
import getpass
import datetime
class EZ:
    cj = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
    login_url = "https://www.ezbiocloud.net/user/login"
    submit_data_url = "https://www.ezbiocloud.net/cl16s/submit_identify_data"
    get_jobs_url = "https://www.ezbiocloud.net/cl16s/get_user_jobs"
    delete_jobs_url = "https://www.ezbiocloud.net/cl16s/delete_jobs"
    def __init__(self, baseDir):
        self.baseDir = baseDir
        outputTemp = path.join(self.baseDir,"ez_temp")
        if not path.exists(outputTemp):
            os.makedirs(outputTemp)
        self.baseDir = outputTemp

    def login(self,username,password):
        formData = {
            "txtID": username,
            "txtPWD": password
        }
        formDataBytes = urllib.parse.urlencode(formData).encode("utf-8")
        req = urllib.request.Request(self.login_url, formDataBytes)
        req.add_header('User-Agent',
                       'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.159 Safari/537.36')
        req.add_header('Content-Type', 'application/x-www-form-urlencoded')

        res = self.opener.open(req)
        json_data = json.loads(res.read().decode())
        print(json_data)

    def submitData(self, fastaFile):

        records = SeqIO.parse(fastaFile, "fasta")
        records = list(records)

        count = 0
        flag = 0

        if path.exists(path.join(self.baseDir,"ez_status.json")):
            with open(path.join(self.baseDir,"ez_status.json"), "r") as f:
                statusDict = json.load(f)
                count = statusDict["count"]
                flag = statusDict["flag"]
        while count < len(records):
            formData = {
                "jsonStr": []
            }
            for i in range(flag*10, (flag+1)*10):


                rec = records[i]
                seq_obj = {
                    "strain_name":rec.id,
                    "ssurrn_seq":str(rec.seq)
                }
                formData["jsonStr"].append(seq_obj)
                count += 1
                print(count)
                if count == len(records):
                    break
            formDataBytes = urllib.parse.urlencode(formData).encode("utf-8")
            req = urllib.request.Request(self.submit_data_url, formDataBytes)
            req.add_header('User-Agent',
                           'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.159 Safari/537.36')
            req.add_header('Content-Type', 'application/x-www-form-urlencoded')

            res = self.opener.open(req)
            json_data = json.loads(res.read().decode())
            print(datetime.datetime.now(), json_data)
            flag += 1
            if count%150 == 0 or count == len(records):
                idList = self.getJobs(count)
                with open(path.join(self.baseDir,"ez_status.json"),"w") as f:
                    json.dump({"count":count,"flag":flag},f)
                if count != len(records):
                    self.deleteJobs(idList)

    def getJobs(self, count):
        req = urllib.request.Request(self.get_jobs_url)
        req.add_header('User-Agent',
                       'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.159 Safari/537.36')
        req.add_header('Content-Type', 'application/x-www-form-urlencoded')

        res = self.opener.open(req)
        json_data = json.loads(res.read().decode())["data"]
        with open(path.join(self.baseDir, f"{(math.ceil(count/150)-1)*150 + 1}-{count}.ez.json"), "w") as f:
            json.dump(json_data, f, indent=4)
        return [i["sge_job_id"] for i in json_data ]


    def deleteJobs(self, idList):
        formData = {
            "jobs": idList

        }
        formDataBytes = urllib.parse.urlencode(formData).encode("utf-8")
        req = urllib.request.Request(self.delete_jobs_url, formDataBytes)
        req.add_header('User-Agent',
                       'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.159 Safari/537.36')
        req.add_header('Content-Type', 'application/x-www-form-urlencoded')
        while True:
            try:
                res = self.opener.open(req)
                json_data = json.loads(res.read().decode())
                print(json_data)
                break
            except HTTPError as e:
                print(e)
                print(e.fp.read().decode())
                time.sleep(1)


    def deleteAllJobs(self):
        req = urllib.request.Request(self.get_jobs_url)
        req.add_header('User-Agent',
                       'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.159 Safari/537.36')
        req.add_header('Content-Type', 'application/x-www-form-urlencoded')

        try:
            res = self.opener.open(req)
        except HTTPError as e:
            print(e.fp.read().decode())
        json_data = json.loads(res.read().decode())["data"]

        idList =  [i["sge_job_id"] for i in json_data]
        self.deleteJobs(idList)


    def mergeData(self):
        """合并json文件，生成qiime2需要的格式"""
        jsonList = os.listdir(self.baseDir)
        jsonList = [i for i in jsonList if "ez.json" in i]
        final_data_f = open(path.join(self.baseDir, "ez_taxonomy.tsv"), "w")
        for j in jsonList:
            with open(path.join(self.baseDir, j), "r", encoding='utf-8') as f:
                data = json.load(f)

                for i in data:
                    index = i["result_taxonomy"].rfind(";;")
                    tax = i["result_taxonomy"][:index]
                    tax = re.sub(r"(;;;)|(;;)|(;)", "; ", tax)
                    print(i["strain_name"] + "\t" + tax + "\t" + str(i["result_similarity"]), file=final_data_f)
        final_data_f.close()

    def mergeData2(self):
        """合并json文件，生成包含全部信息的tsv"""
        jsonList = os.listdir(self.baseDir)
        jsonList = [i for i in jsonList if "ez.json" in i]
        final_data_f = open(path.join(self.baseDir, "id_taxonomy.tsv"), "w")
        print("\t".join(["strain_name", "tax", "result_taxon", "result_strain", "result_similarity", "completeness"]), file=final_data_f)
        for j in jsonList:
            with open(path.join(self.baseDir,j), "r", encoding='utf-8') as f:
                data = json.load(f)

                for i in data:
                    index = i["result_taxonomy"].rfind(";;")
                    tax = i["result_taxonomy"][:index]
                    tax = re.sub(r"(;;;)|(;;)", ";", tax)
                    strain_name = i["strain_name"]
                    completeness = i["completeness"]
                    result_taxon = i["result_taxon"]
                    result_strain = i["result_strain"]
                    result_similarity = i["result_similarity"]
                    col_list =  [strain_name, tax, result_taxon, result_strain, str(result_similarity), str(completeness)]
                    print("\t".join(col_list), file=final_data_f)
        final_data_f.close()



if __name__ == '__main__':
    usage = "python %prog -f/--fastaFile <fasta file> -o/--outputPath <output annotation file path>"
    parser = optparse.OptionParser(usage,
                                   description="Example: python %prog -f example.fasta -o ./")  ## 写入上面定义的帮助信息
    parser.add_option('-f', '--fastaFile', dest='fastaFile', type='string', default="", help='fastaFile')
    parser.add_option('-o', '--outputPath', dest='outputPath', type='string', default="", help='output annotation file path')
    options, args = parser.parse_args()
    ez = EZ(options.outputPath)
    print("\033[0;34;1m{}\033[0m".format(">>>>Please input your account of ezbiocloud:"))
    username = input("username:")
    password = getpass.getpass("password:")
    ez.login(username,password)
    ez.deleteAllJobs()
    ez.submitData(options.fastaFile)
    ez.mergeData()


