from __future__ import division
from subprocess import Popen,PIPE
from shlex import split
import argparse
import os
import pandas as pd
import re
import numpy as np
import sys
from cycler import cycler
import warnings
warnings.filterwarnings("ignore")


try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
#    from matplotlib.backends.backend_pdf import PdfPages
except ImportError as err:
    print err

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


def parameter_finder(target_list, search_list, msgflag=False,exact=False):
    """
    Finds all the elements in the search list in the target list
    the search is non case sensitive

    Args:
        msg is True to tell that there were occurrences or not.
        exact True if the search requires the exact string, False only if it requires to contain the search list, example "vx" is in "vx_sol" 
    Returns:
        indexes with the positions of the search_list found in the target_list
    """
    target_list=map(lambda x:x.lower(),target_list)

    indexes=[]


    if isinstance(search_list, str):
        cont=0
        search_list=search_list.lower()
        for t in target_list:
            if exact==False and search_list in t:
                indexes.append(cont)
            elif exact==True and search_list==t:
                indexes.append(cont)
            cont+=1
    if isinstance(search_list,list):
        search_list=map(lambda x:x.lower(),search_list)

        for s in search_list:
            for cont,t in enumerate(target_list):
                if exact==False and s in t:
                    print "Wrdong"
                    print s,t
                    indexes.append(cont)
                elif exact==True and s==t:
                    print s,t
                    indexes.append(cont)

    
    if msgflag==True:
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

            try:
                data=pd.read_csv(data_file,sep=" ",header=None).dropna(axis=1,how='all')
                data.columns=header.split()
            except:
                raise Exception("The input file '%s' is corrupted, usually the problem is because "\
                                "there is an end of a line that has an additional space" %input_file)

    return data

def extract_digits(strings):
    """
    input:
        strings: Array or single string
    Returns:
         An array of all the digits in a string
         if it is an array of strings, extracts the numbers and returns them sorted

    """
    if isinstance(strings, str):
        numbers=re.findall(r"[-+]?\d*\.?\d+",strings)
        output=([num.strip('.') for num in numbers])


    if isinstance(strings, list):
        output=[]
        for element in strings:
            numbers=re.findall(r"[-+]?\d*\.?\d+",element)
            output.append([num.strip('.') for num in numbers])

        if np.shape(output)[1]==1:
            """
            If there are several parameters it is difficult to reshape
            """
            output=np.sort(np.array(output,dtype=float).reshape((len(output))))

    return np.array(output,dtype=float)


def blockPrint():
    """
    Redirecting stdout to nothing
    """
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def set_plot_appearance():

    """
    Defines the appearence of the plot

    use rcParams.keys() to see the available parameters
    """
    axis_font=24
    tick_font=20
    legend_font=18


    #plt.rcParams['lines.linewidth'] = 1.5
    #plt.rcParams['lines.markeredgewidth'] = 1.0
    #plt.rcParams['lines.markersize'] = 2.5

    # Colors
    plt.rcParams['axes.prop_cycle'] = cycler(color='rbkgcmy')


    # Fonts and symbols
    #plt.rcParams['font.family'] = 'serif'
    #plt.rcParams['font.serif'] = 'Times New Roman'
    #plt.rcParams['font.weight'] = 'normal'
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams['text.usetex'] = True
    #plt.rcParams['mathtext.rm'] = 'serif'
    #plt.rcParams['mathtext.it'] = 'serif:italic'
    #plt.rcParams['mathtext.fontset'] = 'stix'


    # Axes
    plt.rcParams['axes.edgecolor'] = (0.0, 0.0, 0.0)
    #plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['axes.spines.right'] = True
    plt.rcParams['axes.spines.top'] = True
    plt.rcParams['axes.labelsize'] = axis_font
    plt.rcParams['axes.grid'] = False

    # Ticks
    plt.rcParams['xtick.direction'] = "in"
    plt.rcParams['ytick.direction'] = "in"
    plt.rcParams['ytick.labelsize'] = tick_font
    plt.rcParams['xtick.labelsize'] = tick_font

    # Legend

    plt.rcParams['legend.fontsize'] = legend_font
    plt.rcParams['legend.loc'] ='upper left'
    plt.rcParams['legend.labelspacing'] = 0.5
    plt.rcParams['legend.borderpad'] =0.4
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.edgecolor'] = 'k'
    # Fonts and symbols


    # Errorbar plots
    plt.rcParams['errorbar.capsize'] = 4

def integrate(x,y,xmin,xmax):
    """
    Integrate the data in x and y from xmin to xmax
    Args:
        x are the independent variables
        y the integrand
        xmin the lower limit of the integral
        xmax the upper limit of the integral
    """
    MinIndex=np.min(np.where(x>=xmin))
    MaxIndex=np.max(np.where(x<=xmax))
    I=np.trapz(y[MinIndex:MaxIndex],x[MinIndex:MaxIndex])

    return I



def modify_file(file_name,key_word,modified_line,copy_name=None,n_ocurrence=0):
    """
    Replaces the entire line of the n-th occurrence of the key in the file
    Args:
        file_name Path to the file to modify
        key_word to find
        modified_line text of the entire line
        copy_name if the output file is going to be different from file_name
        n_ocurrence the ocurrence of the key that is of interest
        
    """
    
    
    if copy_name==None:
        copy_name=file_name
    f=open(file_name,'r')
    
    lines=f.readlines()
    line_number=parameter_finder(lines,key_word)[n_ocurrence]
    lines[line_number]=modified_line
    f.close()
    s=open(copy_name,'w')
    s.writelines(lines)
    s.close()
    
    print "Modified the file %s"%(file_name.split('/')[-1])