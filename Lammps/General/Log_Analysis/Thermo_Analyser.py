"""
This script reads the log file and gets the
"""
import numpy as np
import pandas as pd
import os
import sys
import linecache
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat
"""
Initialisation and controling inputs
"""


def read_chunk(file_name,i_line,f_line):
    data=[]
    for i in range(i_line,f_line):
        data.append(linecache.getline(file_name, i+1).strip('\n').split())
    data=np.array(data,dtype='double')
    return data

def save_file(data,header):
    """
    Writes the file
    """
    np.savetxt("Parameters.dat",data,header=header)
    print("Generated the file Parameters.dat")



def data_extract(file_name):
    """
    Need to change this For the time being use gawk in os and awk in ubuntu.
    This function extracts the data from the timesteps after discarding the defined amount
    
    """
    if sys.platform=='darwin':
        awk_cmd='gawk'
    else:
        awk_cmd='awk'
    out,err=cf.bash_command("""grep -n "Per MPI" %s| %s -F":" '{print $1}'"""%(file_name,awk_cmd) )
    initial_line=np.array(out.split(),dtype=int)+1

    out2,err=cf.bash_command("""grep -n "Loop time" %s| %s -F":" '{print $1}'"""%(file_name,awk_cmd))
    final_line=np.array(out2.split(),dtype=int)-1
    #Check if there was minimization
    out3,err=cf.bash_command("""grep -n "Minimization stats" %s"""%file_name)

    if out3:
        initial_line=np.delete(initial_line,0)
        final_line=np.delete(final_line,0)

    """Checking if the simulation did not finish"""
    if np.size(initial_line)>np.size(final_line):
        print("\nThe log file is incomplete!, lets analyse it in anycase \n")
        out3,err=cf.bash_command("""wc -l %s |awk '{print $1}'"""%file_name)
        line=int(out3.split()[0])
        final_line=np.append(final_line,line)


    header=linecache.getline(file_name, initial_line[0]).strip('\n').split()
    header_string=" ".join(header)
    number_chunks=np.size(initial_line)

    for i in range(number_chunks):
        if i==0:
            total=read_chunk(file_name, initial_line[i],final_line[i])
        else:
            data=read_chunk(file_name,initial_line[i],final_line[i])
            total=np.delete(total,-1,axis=0) #As there are repeated timesteps, I choose to keep the last version of the variables
            total=np.vstack([total,data])


    return total, header_string

def discard_data(data,nmin):
    """
    Function to discard in number or percentage
    Args:
        data is a numpy array containing the data
        nmin could be a percentage between 0-1 or an integer
    """

    if nmin<1:
        print("Discarding %d%% of the timesteps for the analysis"%(int(nmin*100)))
        discard=int(nmin*len(data))
        data=data[discard:]
    else:
        print("Discarding %d out of %d timesteps for the analysis" %(nmin,len(data)))
        data=data[int(nmin):]

    return data

def thermo_analyser(file_name,minimum = 0.3):
    """
    Args:
        file_name
        min Number or percentage (between 0-1) of samples to be discarded
    Returns:
        thermo_data a dataframe that has the columns as the average, etc and the rows for each property
        check the columns with thermo_data.index.values, thermo_data.column.values 
    """

    data,header=data_extract(file_name)
    save_file(data,header)
    stat.fast_averager("Parameters.dat",minimum,"thermo.dat")

    thermo_data = read_thermo_output("thermo.dat")

    return thermo_data


def read_thermo_output(file_thermo):
    """
    Creates a data frame from the file "file_thermo"

    Returns:
        df a dataframe that has the columns as the average, etc and the rows for each property
        check the columns with thermo_data.index.values, thermo_data.column.values
    """

    f=open(file_thermo,'r')
    g = open("thermo_mod.dat",'w')
    lines=f.readlines()
    separator = "="
    array = []
    for l in lines:
        l = l.replace("=","\t")
        l = l.split()
        array.append(l)

    df = pd.DataFrame(array[1::], columns = array[0])
    df = df.set_index('Property')

    return df



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
    parser.add_argument('FileName', metavar='InputFile',help='Input filename',type=lambda x: cf.is_valid_file(parser, x))
    parser.add_argument('--min', help='Number or percentage (between 0-1) of samples to be discarded', default=0, type=float)
    args = parser.parse_args()

    thermo_analyser(args.FileName,args.min)

    
    
    