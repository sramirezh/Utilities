from __future__ import division
from subprocess import Popen,PIPE
from shlex import split


def bash_command(cmd):
    """
    function that evaluates simple bash commands with pipes,
    
    Input:
        cmd is a string of the command just as you write it on the shell
    Returns: two elements
        out The output 
        err the possible errors
    """
    chain=cmd.split("|")
    n_pipes=len(chain)
    
    for i in xrange(n_pipes):
        if i==0:
            p=Popen(split(chain[0]),stdout=PIPE)
        else:
            p=Popen(split(chain[i]), stdin=p.stdout, stdout=PIPE)        
    
    return p.communicate()

    
    