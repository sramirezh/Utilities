#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 18:08:10 2018
Script to display the qstat as I want.

@author: simon
"""

import pandas as pd
from lxml import etree
import re
from subprocess import Popen,PIPE

f=open("out.qstat","w+")
p=Popen(['qstat', '-x'],stdout=f)
p.wait() 
f.close()
data="out.qstat"

tree = etree.parse(data)

names = []
values = []
for p in tree.iter() :
    names.append(tree.getpath(p).replace("/",".")[1:])
    values.append(p.text)
    
#df = pd.DataFrame({'key' : lstKey, 'value' : lstValue})
#df.sort_values('key')


names=names[1::]
values=values[1::]

#parameter_names=[]
#for n in names:
#    string=n.split('.',2)[-1]
#    if 'Job[' in string: #To not include this patter
#        continue
#    parameter_names.append(string)
#    
#parameter_names=list(set(parameter_names)) #Deleting the repeated occurrences

"""
'find my jobs"
"""
job_numbers=[]
for i,val in enumerate(values,start=0):
    if val is None:
        continue
    if 'sr802@dexter.cm' in val:
        job_numbers.append(re.findall(r"[-+]?\d*\.?\d+", names[i])[0])
        
"""
#Building the table
"""

table_names=["Job_Id", "Job_Name","queue", "Resource_List.nodes","job_state", "resources_used.walltime", "Resource_List.walltime"]
print_names=["Job_Id", "Job_Name","queue", "Resources","job_state", "Running Time", "Walltime"]
jobs=[]
for job in job_numbers:
    string= "Data.Job[%s]" %job
    index=names.index(string) #Find where this data starts
    del values[0:index]
    del names[0:index]
    job_values=[]
    for parameter in table_names:
        string_parameter=string+"."+parameter
        try:
            index2=names.index(string_parameter)
            text=values[index2].replace(".dexter.cm.cluster","")
        except ValueError:
            text="--"
        job_values.append(text) #This just applies to Job_Id, but I am lazy
    jobs.append(job_values)

df=pd.DataFrame(jobs,columns=print_names)
print(df)
    








    


