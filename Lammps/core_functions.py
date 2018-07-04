from __future__ import division
from subprocess import Popen,PIPE
from shlex import split
import argparse
import os
import pandas as pd

def bash_command(cmd):
    """
    function that evaluates simple bash commands with pipes,

    Input:
        cmd is a string of the command just as you write it on the shell but inside 3" in case you have
        several quoted text inside
    Returns: two elements
        out The output
        err the possible errors
    """
    cmd=str(cmd)
    chain=cmd.split("|")
    n_pipes=len(chain)

    for i in xrange(n_pipes):
        if i==0:
            p=Popen(split(chain[0]),stdout=PIPE)
        else:
            p=Popen(split(chain[i]), stdin=p.stdout, stdout=PIPE)

    return p.communicate()


def is_valid_file(parser, arg):
    """
    Checks if a file exists and returns the proper output to the parser handler.
    """
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def parameter_finder(List, String):
    """
    Finds a string on a List and returns the position on the list
    It is case insensitive
    """
    List=map(lambda x:x.lower(),List)
    String=String.lower()
    cont=0
    indexes=[]
    for s in List:
        if String in s:
            indexes.append(cont)
        cont+=1
    length=len(indexes)
    if length>1: print "There were several ocurrences"
    if length==0: print "No ocurrences found"

    return indexes

def read_data_file(input_file):
    """
    Reads a data file either with a header or not.
    It assumes that the header is commented with "#" and that the it last line contains the name of the variables
    Args:
        input_file file name

    Returns: A panda data frame, the column names can be obtained by data.columns.values and the numeric parameters with  data.values
    """
    header_lines=0
    last_pound_pos=-1
    with open(input_file, 'r') as data_file:
        while(data_file.read(1)=='#'):
            last_pound_pos = data_file.tell()
            header=data_file.readline()
            header_lines+=1

        #Read the next lines
        data_1=data_file.readline().split()
        data_2=data_file.readline().split()
        data_file.seek(last_pound_pos+1) #Goes back to the last line of the header

        if header_lines==0:
            data=pd.read_csv(data_file,sep=" ",header=None).dropna(axis=1,how='all')

        else:
            if len(data_1)!=len(data_2): #If there is a line containing the number of particles,
                data_file.readline()
            data_file.readline()
            data=pd.read_csv(data_file,sep=" ",header=None).dropna(axis=1,how='all')
            data.columns=header.split()

    return data
