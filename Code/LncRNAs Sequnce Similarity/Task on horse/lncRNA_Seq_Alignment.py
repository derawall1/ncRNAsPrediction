from __future__ import print_function

import itertools
import numpy as np
from Bio import Align
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import csv
import pandas as pd
from Bio.pairwise2 import format_alignment
import math
#from tqdm import tqdm


# In[2]:



import os
import sys
import time
import requests
import platform
from xmltramp2 import xmltramp
from optparse import OptionParser
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
    from urllib.request import __version__ as urllib_version
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError
    from urllib2 import __version__ as urllib_version

# allow unicode(str) to be used in python 3
try:
    unicode('')
except NameError:
    unicode = str


# In[29]:


# Base URL for service
baseUrl = u'https://www.ebi.ac.uk/Tools/services/rest/emboss_water'
version = u'2019-07-03 12:51'
# scratch1 path
scratch1_path = '/scratch1/Mushtaq/LncRNA_Similarity/'
# Set interval for checking status
pollFreq = 3
# Output level
outputLevel = 1
# Debug level
debugLevel = 0
# print degug message
def printDebugMessage(functionName, message, level):
    if (level <= debugLevel):
        print(u'[' + functionName + u'] ' + message, file=sys.stderr)

# User-agent for request (see RFC2616).
def getUserAgent():
    printDebugMessage(u'getUserAgent', u'Begin', 11)
    # Agent string for urllib2 library.
    urllib_agent = u'Python-urllib/%s' % urllib_version
    clientRevision = version
    # Prepend client specific agent string.
    try:
        pythonversion = platform.python_version()
        pythonsys = platform.system()
    except ValueError:
        pythonversion, pythonsys = "Unknown", "Unknown"
    user_agent = u'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
        clientRevision, os.path.abspath(''),
        pythonversion, pythonsys, urllib_agent)
    printDebugMessage(u'getUserAgent', u'user_agent: ' + user_agent, 12)
    printDebugMessage(u'getUserAgent', u'End', 11)
    return user_agent

# Get available result types for job
def serviceGetResultTypes(jobId):
    printDebugMessage(u'serviceGetResultTypes', u'Begin', 1)
    printDebugMessage(u'serviceGetResultTypes', u'jobId: ' + jobId, 2)
    requestUrl = baseUrl + u'/resulttypes/' + jobId
    printDebugMessage(u'serviceGetResultTypes', u'requestUrl: ' + requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    printDebugMessage(u'serviceGetResultTypes', u'End', 1)
    return doc[u'type':]

# Wrapper for a REST (HTTP GET) request
def restRequest(url):
    printDebugMessage(u'restRequest', u'Begin', 11)
    printDebugMessage(u'restRequest', u'url: ' + url, 11)
    try:
        # Set the User-agent.
        user_agent = getUserAgent()
        http_headers = {u'User-Agent': user_agent}
        req = Request(url, None, http_headers)
        # Make the request (HTTP GET).
        reqH = urlopen(req)
        resp = reqH.read()
        contenttype = reqH.info()

        if (len(resp) > 0 and contenttype != u"image/png;charset=UTF-8"
                and contenttype != u"image/jpeg;charset=UTF-8"
                and contenttype != u"application/gzip;charset=UTF-8"):
            try:
                result = unicode(resp, u'utf-8')
            except UnicodeDecodeError:
                result = resp
        else:
            result = resp
        reqH.close()
    # Errors are indicated by HTTP status codes.
    except HTTPError as ex:
        result = requests.get(url).content
    printDebugMessage(u'restRequest', u'End', 11)
    return result


# Get job status
def serviceGetStatus(jobId):
    printDebugMessage(u'serviceGetStatus', u'Begin', 1)
    printDebugMessage(u'serviceGetStatus', u'jobId: ' + jobId, 2)
    requestUrl = baseUrl + u'/status/' + jobId
    printDebugMessage(u'serviceGetStatus', u'requestUrl: ' + requestUrl, 2)
    status = restRequest(requestUrl)
    printDebugMessage(u'serviceGetStatus', u'status: ' + status, 2)
    printDebugMessage(u'serviceGetStatus', u'End', 1)
    return status

# Client-side poll
def clientPoll(jobId):
    printDebugMessage(u'clientPoll', u'Begin', 1)
    result = u'PENDING'
    while result == u'RUNNING' or result == u'PENDING':
        result = serviceGetStatus(jobId)
        if outputLevel > 0:
            print(result)
        if result == u'RUNNING' or result == u'PENDING':
            time.sleep(pollFreq)
    printDebugMessage(u'clientPoll', u'End', 1)

# Get result
def serviceGetResult(jobId, type_):
    printDebugMessage(u'serviceGetResult', u'Begin', 1)
    printDebugMessage(u'serviceGetResult', u'jobId: ' + jobId, 2)
    printDebugMessage(u'serviceGetResult', u'type_: ' + type_, 2)
    requestUrl = baseUrl + u'/result/' + jobId + u'/' + type_
    result = restRequest(requestUrl)
    printDebugMessage(u'serviceGetResult', u'End', 1)
    return result


def getResult(jobId,outfile=None):
    printDebugMessage(u'getResult', u'Begin', 1)
    printDebugMessage(u'getResult', u'jobId: ' + jobId, 1)
    if outputLevel > 1:
        print("Getting results for job %s" % jobId)
    # Check status and wait if necessary
    clientPoll(jobId)
    # Get available result types
    resultTypes = serviceGetResultTypes(jobId)

    for resultType in resultTypes:
        # Derive the filename for the result
        if outfile:
            filename = scratch1_path +(outfile + u'.' + unicode(resultType[u'identifier']) +
                        u'.' + unicode(resultType[u'fileSuffix']))
        else:
            filename = scratch1_path +(jobId + u'.' + unicode(resultType[u'identifier']) +
                        u'.' + unicode(resultType[u'fileSuffix']))
        # Write a result file

        outformat_parm = str('aln').split(',')
        for outformat_type in outformat_parm:
            outformat_type = outformat_type.replace(' ', '')
            print(unicode(resultType[u'identifier']))

            if outformat_type == 'None':
                outformat_type = None

            if not outformat_type or outformat_type == unicode(resultType[u'identifier']):
                if outputLevel > 1:
                    print("Getting %s" % unicode(resultType[u'identifier']))
                # Get the result
                result = serviceGetResult(jobId, unicode(resultType[u'identifier']))
                if (unicode(resultType[u'mediaType']) == u"image/png"
                        or unicode(resultType[u'mediaType']) == u"image/jpeg"
                        or unicode(resultType[u'mediaType']) == u"application/gzip"):
                    fmode = 'wb'
                else:
                    fmode = 'w'

                try:
                    fh = open(filename, fmode)
                    fh.write(result)
                    fh.close()
                except TypeError:
                    fh.close()
                    fh = open(filename, "wb")
                    fh.write(result)
                    fh.close()
                if outputLevel > 0:
                    print("Creating result file: " + filename)
    printDebugMessage(u'getResult', u'End', 1)
    
    
    # Submit job
def serviceRun(email, title, params):
    printDebugMessage(u'serviceRun', u'Begin', 1)
    # Insert e-mail and title into params
    params[u'email'] = email
    if title:
        params[u'title'] = title
    requestUrl = baseUrl + u'/run/'
    printDebugMessage(u'serviceRun', u'requestUrl: ' + requestUrl, 2)

    # Get the data for the other options
    requestData = urlencode(params)

    printDebugMessage(u'serviceRun', u'requestData: ' + requestData, 2)
    # Errors are indicated by HTTP status codes.
    try:
        # Set the HTTP User-agent.
        user_agent = getUserAgent()
        http_headers = {u'User-Agent': user_agent}
        req = Request(requestUrl, None, http_headers)
        # Make the submission (HTTP POST).
        reqH = urlopen(req, requestData.encode(encoding=u'utf_8', errors=u'strict'))
        jobId = unicode(reqH.read(), u'utf-8')
        reqH.close()
    except HTTPError as ex:
        print(xmltramp.parse(unicode(ex.read(), u'utf-8'))[0][0])
        quit()
    printDebugMessage(u'serviceRun', u'jobId: ' + jobId, 2)
    printDebugMessage(u'serviceRun', u'End', 1)
    getResult(jobId,title)
    return jobId    




def get_alignment(filename):
    inputFile = open(scratch1_path+filename,"r")
    text = inputFile.readlines()
    inputFile.close()
    align_dict = {}
    for line in text:
        #print(line)
        if ':' in line:
            line = line.replace('#','')
            key, value = line.split(":",1)
            align_dict[key.strip()] = value.strip()

    return align_dict






query_seq = list(SeqIO.parse("NP_V2_lncRNA_seq_clean_616.fasta", "fasta"))
target_seq = list(SeqIO.parse("NP_V2_lncRNA_seq_clean2.fasta", "fasta"))



import decimal

def round_down(value, decimals):
    with decimal.localcontext() as ctx:
        d = decimal.Decimal(value)
        ctx.rounding = decimal.ROUND_DOWN
        return round(d, decimals)


# In[39]:


def calculate_similarity(query_seq,target_seq):
    print("started task...")
    lnc_set = set()
    with open(scratch1_path + "LncRNA-LncRNA_Similarity_616_1923.csv", "w") as csvFile:
        # csv file header
        fieldnames = ['Query_Seq_ID','Target_Seq_ID', 'Align_Score','Normalize_Score', 'Identity', 'Similarity']
        writer = csv.DictWriter(csvFile, fieldnames=fieldnames)
        writer.writeheader()
        # loop on query seq
        count=0
        for q in range(0,len(query_seq),1):
                        
            try:
                
                seq_q = query_seq[q]
                # self alignment with query seq
                params = {}
                params[u'asequence'] = seq_q.seq
                params[u'bsequence'] = seq_q.seq
                params['stype'] = 'dna'
                #params['matrix'] = 'dnafull'
                params['gapopen'] = '10'
                params['gapext'] = '0.5'
                params['format'] = 'pair'
                #params['format'] = options.format
                job_id_q = serviceRun('derawall@gmail.com',seq_q.id + '-' + seq_q.id,params)
                seq_q_align = get_alignment(seq_q.id + '-' + seq_q.id+".aln.txt")
                #seq_q_align['Similarity']
                
               
                for t in range(0,len(target_seq),1):
                    seq_t = target_seq[t]
                    lncRna_set = (seq_q.id,seq_t.id)
                    lncRna_set2 = (seq_t.id,seq_q.id)
                              
                              
                              
                    if (seq_q.seq != seq_t.seq) and ((lncRna_set not in lnc_set) or (lncRna_set2 not in lnc_set)):
                        
                        try:
                            lnc_set.add(lncRna_set)
                            lnc_set.add(lncRna_set2)

                            # self alignment with traget seq
                            params = {}
                            params[u'asequence'] = seq_t.seq
                            params[u'bsequence'] = seq_t.seq
                            params['stype'] = 'dna'
                            #params['matrix'] = 'dnafull'
                            params['gapopen'] = '10'
                            params['gapext'] = '0.5'
                            params['format'] = 'pair'
                            #params['format'] = options.format
                            job_id_t = serviceRun('derawall@gmail.com',seq_t.id + '-' + seq_t.id,params)
                            seq_t_align = get_alignment(seq_t.id + '-' + seq_t.id+".aln.txt")
                              # getting the score of query seq self align
                            seq_q_align_score=float(seq_q_align['Score'])
                            print(seq_q_align_score)
                            # getting the score of target seq self align
                            seq_t_align_score=float(seq_q_align['Score'])
                            print(seq_t_align_score)


                            # the alignment b/w qry and trg seq
                            params = {}
                            params[u'asequence'] = seq_q.seq
                            params[u'bsequence'] = seq_t.seq
                            params['stype'] = 'dna'
                            #params['matrix'] = 'dnafull'
                            params['gapopen'] = '10'
                            params['gapext'] = '0.5'
                            params['format'] = 'pair'
                            #params['format'] = options.format
                            job_id_align = serviceRun('derawall@gmail.com',seq_q.id + '-' + seq_t.id,params)
                            seq_q_t_align = get_alignment(seq_q.id + '-' + seq_t.id+".aln.txt")
                            # getting the score of query and target seq
                            score =float(seq_q_t_align['Score'])
                            print(score)
                            # normalizing the score
                            Normalize_Score = round_down(float(score)/(math.sqrt(seq_q_align_score)*math.sqrt(seq_t_align_score)),2)
                            print(Normalize_Score)

                            Identity = seq_q_t_align['Identity']
                            #print(Identity)

                            #Identity =str(round_down(matches/(matches + mismatches) * 100,2)) +'%' 
                            # calculating the similarity in %
                            Similarity=seq_q_t_align['Similarity']
                            #Similarity=(matches/seq_len * 100)
                            #print(Similarity)
                            writer.writerow({'Query_Seq_ID': seq_q.id,                                      
                                             'Target_Seq_ID': seq_t.id,
                                             'Align_Score': score, 'Normalize_Score':Normalize_Score,
                                             'Identity':Identity, 'Similarity':Similarity })
                            writer.writerow({'Query_Seq_ID': seq_t.id,                                      
                                             'Target_Seq_ID': seq_q.id,
                                             'Align_Score': score, 'Normalize_Score':Normalize_Score,
                                             'Identity':Identity, 'Similarity':Similarity })
                            count +=1
                            print(count)

                            print(seq_q.id)
                            print('.............................')

                            print(seq_t.id)
                            #print('length: ' + str(seq_len))
                            #print('length2: ' + str(fn_format_alignment(*align[0])[4]))
                            print('Score: ' + str(score))
                            

                            print('Identity: ' + str(Identity))

                            print('Normalize_Score: ' + str(Normalize_Score))
                            print('Similarity: ' + str(Similarity))



                            print('.............................')

                            #print(align[0])
                            #print('---------------------------------------')
                            #print(align_score)
                            #print(align_score[1])
                            #print(align_score[2])
                            #print(align_score[3])
                            #print(type(align_score))  

                        except Exception as e:
                            print(e)
                            #break
                            continue # doing nothing on exception
            #count+=1
            #print(count)
            #print(seq_q.id.split('|')[1] + '  Completed')
            except Exception as e:
                    print(e)
                    continue # doing nothing on exception

    print("Task Completed.....")




print("started....")
calculate_similarity(query_seq,target_seq)
print("completed.....")





