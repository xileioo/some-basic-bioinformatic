# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:26:28 2022

@author: xileioo
"""

import pandas as pd
import re

fileOut = open('d:/work/reference/human/gencode.v39.gene.txt','w')
file = open('d:/work/reference/human/gencode.v39.annotation.gff3','r')
lines = file.readlines()[7:]

count = 0

for line in lines:
    line = line.rstrip('\n')
    arr = line.split('\t')
    if(len(arr) == 9 and arr[2] == 'gene'):
        #print(line)
        chr_num = arr[0]
        start = arr[3]
        end = arr[4]
        info = arr[8]
        count += 1
        search_res = re.search(';gene_id=(.*);gene_type=(.*);gene_name=(.*);level=',info)
        final = chr_num + '\t' + start + '\t' + end + '\t' + search_res.group(1) + '\t' + search_res.group(2) + '\t' + search_res.group(3) + '\n'
        #print(final)
        fileOut.writelines(final)
    #if(count > 10):
    #    break
