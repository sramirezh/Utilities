from __future__ import division
from subprocess import Popen,PIPE
from shlex import split
import argparse
import os

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
