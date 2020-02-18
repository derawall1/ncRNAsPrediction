#!/usr/bin/env python
# coding: utf-8

# ### Importing Packages

# In[15]:


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


# ### Reading fast files of protein sequence

# In[35]:



query_seq = list(SeqIO.parse("../Data/NpInter2 Interacting LncRNAs/NP_V2_lncRNA_seq_clean.fasta", "fasta"))
target_seq = list(SeqIO.parse("../Data/NpInter2 Interacting LncRNAs/NP_V2_lncRNA_seq_clean.fasta", "fasta"))
#query_seq = list(SeqIO.parse("../Data/NpInter2 Interacting LncRNAs/target_sample.fa", "fasta"))
#target_seq = list(SeqIO.parse("../Data/NpInter2 Interacting LncRNAs/target_sample.fa", "fasta"))


# In[36]:


print(len(target_seq))
print(len(target_seq))


# ### custom format_alignment to get the match, mismatch, gaps and score

# In[18]:


def fn_format_alignment(align1, align2, score, begin, end): 
    #s = [] 
    #s.append("%s\n" % align1) 
    #s.append(" " * begin) 
    #matrix=blosum62
    gap=0
    mismatch=0
    match=0
    seq_len=0
    mismatch_positive=0
    #print(align1[begin:end])
    #print(align2[begin:end])
    for a, b in zip(align1[begin:end], align2[begin:end]): 
        if a == b: 
            #s.append("|")  # match 
            match+=1
        elif a == "-" or b == "-":
            
            #s.append(" ")  # gap
            gap+=1
        elif a != b:
#             if (a,b) in blosum62:
#                 if blosum62[(a,b)] >0:
#                     mismatch_positive+=1
#             else:
#                 if blosum62[(b,a)] >0:
#                     mismatch_positive+=1
            
                                                     
            #s.append(".")  # mismatch 
            mismatch+=1
        seq_len = min(len(align1[begin:end]), len(align2[begin:end]))
    #s.append("\n") 
    #s.append("%s\n" % align2) 
    #s.append("  Score=%g\n" % score)
    #s.append("  Gap=%g\n" % gap)
    #s.append("  Mismatch=%g\n" % mismatch)
    return score,match,mismatch,gap,seq_len,mismatch_positive


# In[19]:


# round decimal number

import decimal

def round_down(value, decimals):
    with decimal.localcontext() as ctx:
        d = decimal.Decimal(value)
        ctx.rounding = decimal.ROUND_DOWN
        return round(d, decimals)


# ### Creating Empty CSV file to store the similarity results

# #### Calculating similarity score and saving in csv

# In[33]:



def calculate_similarity(query_seq,target_seq):
    print("started task...")
    with open("../Data/NpInter2 Interacting LncRNAs/LncRNA-LncRNA Sequence Similarity/LncRNA-LncRNA_Similarity.csv", "a", newline='') as csvFile:
        # csv file header
        fieldnames = ['Query_Seq_ID','Target_Seq_ID', 'Align_Score','Normalize_Score', 'Identity', 'Similarity']
        writer = csv.DictWriter(csvFile, fieldnames=fieldnames)
        writer.writeheader()
        # loop on query seq
        count=0
        for seq_q in query_seq:
            try:
                # self alignment with query seq
                seq_q_align = pairwise2.align.localms(seq_q.seq, seq_q.seq, 5, -4, -10, -0.5, one_alignment_only=True)
                #print(seq_q_align)
                # loop on target seq
                for seq_t in target_seq:
                    if seq_q.seq != seq_t.seq:
                        
                        try:
                            # self alignment with traget seq
                            seq_t_align = pairwise2.align.localms(seq_t.seq, seq_t.seq, 5, -4, -10, -0.5, one_alignment_only=True)
                            # getting the score of query seq self align
                            #print(seq_t_align)
                            #print(fn_format_alignment(*seq_q_align[0]))
                            seq_q_align_score=fn_format_alignment(*seq_q_align[0])[0]
                            #print(seq_q_align_score)
                            # getting the score of target seq self align
                            seq_t_align_score=fn_format_alignment(*seq_t_align[0])[0]
                            #print(seq_t_align_score)
                            

                            # the alignment b/w qry and trg seq
                            align = pairwise2.align.localms(seq_q.seq, seq_t.seq, 5, -4, -10, -0.5, one_alignment_only=True)
                            # getting the score of query and target seq
                            score =fn_format_alignment(*align[0])[0]
                            #print(score)
                            # normalizing the score
                            Normalize_Score = round_down(score/(math.sqrt(seq_q_align_score)*math.sqrt(seq_t_align_score)),2)
                            #print(Normalize_Score)
                            # get the no of matches in query and target seq
                            matches = fn_format_alignment(*align[0])[1]
                            #print(matches)
                            # getting the no of mismatches in qry and trg seq
                            mismatches = fn_format_alignment(*align[0])[2]
                            #print(mismatches)
                            mismatch_positive=fn_format_alignment(*align[0])[5]
                            #print(mismatch_positive)
                            # calcualting the percent identity %.
                            seq_len = fn_format_alignment(*align[0])[4]
                            #print(seq_len)
                            Identity = str(round_down(matches/seq_len * 100,2)) +'%' 
                            #print(Identity)

                            #Identity =str(round_down(matches/(matches + mismatches) * 100,2)) +'%' 
                            # calculating the similarity in %
                            #Similarity=str(round_down((matches + mismatches)/seq_len * 100,2)) + '%'
                            Similarity=str(round_down((matches + mismatch_positive)/seq_len * 100,2)) + '%'
                            #print(Similarity)
                            writer.writerow({'Query_Seq_ID': seq_q.id,                                      
                                             'Target_Seq_ID': seq_t.id,
                                             'Align_Score': score, 'Normalize_Score':Normalize_Score,
                                             'Identity':Identity, 'Similarity':Similarity })
                            count +=1
                            print(count)
                            
                            print(seq_q.id)
                            print('.............................')
                            
                            print(seq_t.id)
                            print('length: ' + str(min(len(seq_q.seq), len(seq_t.seq))))
                            print('length2: ' + str(fn_format_alignment(*align[0])[4]))
                            print('Score: ' + str(fn_format_alignment(*align[0])[0]))
                            print('matches: ' + str(fn_format_alignment(*align[0])[1]))
                            print('mismatches: ' + str(fn_format_alignment(*align[0])[2]))
                            print('mismatches+ : '+ str(mismatch_positive))
                            print('gaps: ' + str(fn_format_alignment(*align[0])[3]))

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







