"""
This script reads the log file and gets the
"""
import numpy as np
import pandas as pd
import os
import sys
import linecache
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat
"""
Initialisation and controling inputs
"""


def read_chunk(file_name,i_line,f_line):
    """
    Args:
    i_lines is the line where the name of the columns are, that is why the loop
    starts at i_lines+1, 
    f_line is the number of the line after the end of the 
    data, so the loop reads until one line before
    """
    data=[]
    
    for i in range(i_line+1,f_line):
        data.append(linecache.getline(file_name, i).strip('\n').split())
    linecache.clearcache()
    data = np.array(data, dtype = float)
    return data


#def read_chunk(file_name, i_line, f_line):
#    """
#    TODO: This function is 10x slower than the one with linecache but is more 
#    reliable, what is the problem of the linecache function above?
#    
#    i_line: initial line to read
#    f_line: last line to read+1
#    
#    """
#    data = []
##    print (i_line, f_line)
#    with open(file_name) as fp:
#        for i, line in enumerate(fp):
#            if i >= i_line and i<f_line-1:
#                data.append(line.strip('\n').split())
#            else:
#                pass
##                print (i)
##    print (data[0], data[-1])
#    data = np.array(data, dtype = float)
#    fp.close()
#    return data
    
            

def save_file(data,header):
    """
    Writes the file
    """
    np.savetxt("Parameters.dat",data,header=header)
    print("Generated the file Parameters.dat")



def data_extract(file_name):
    """
    Extracts the data from log.lammps type files and if the file did not finish,
    still builds an array with all the results. Notice that the structure from
    a lammps log file is as follows:
        
        Per MPI rank memory allocation (min/avg/max) = 2.743 | 2.782 | 3.106 Mbytes
        Elapsed Step c_tempall c_tempfluid Press Volume Lz
       0        0    1.0089874    1.0482703    1.2968607    10213.444      35.4657
       .
       .
       .
       10000000 10000000   0.99920245    1.0381044     1.358909    10213.444      35.4657
       Loop time of 8850.06 on 16 procs for 10000000 steps with 7713 atoms
        
       
    Args:
        file_name: Name of the log file, usually: log.lammps
    
    returns:
        total: Array with the data
        header_string: name of the variables in the data array

    TODO: At the moment it assummes that the header is the same, so the
    thermo data is the same for all the chunks, it could be generalised so
    if the printed data changed during the run, it could be taken into account.
    
    """
    if sys.platform=='darwin':
        awk_cmd='awk' #gwak was necessary before catalina OS
    else:
        awk_cmd='awk'
       
    # Getting the beginning of a chunk in the log
    out,err = cf.bash_command("""grep -n "Per MPI" %s| %s -F":" '{print $1}'"""%(file_name,awk_cmd) )
    initial_line = np.array(out.split(), dtype=int) + 1
    
    # Getting the end of a chunk
    out2,err=cf.bash_command("""grep -n "Loop time" %s| %s -F":" '{print $1}'"""%(file_name,awk_cmd))
    final_line=np.array(out2.split(), dtype=int)
    
    #Check if there was minimization
    out3,err=cf.bash_command("""grep -n "Minimization stats" %s"""%file_name)
    if out3:
        # Deletes the first chunk containing the minimisation
        initial_line = np.delete(initial_line,0)
        final_line = np.delete(final_line,0)

    # Checks if the simulation finished, if not, the last line of the file is 
    # the last line of the chunk
    if np.size(initial_line) > np.size(final_line):
        print("\nThe log file is incomplete!, lets analyse it in anycase \n")
        out3,err = cf.bash_command("""wc -l %s |awk '{print $1}'"""%file_name)
        line = int(out3.split()[0])
        final_line = np.append(final_line,line)


    header = linecache.getline(file_name, initial_line[0]).strip('\n').split()
    linecache.clearcache()
    header_string = " ".join(header)
    number_chunks = np.size(initial_line)

    for i in range(number_chunks):
        if i == 0:
            total = read_chunk(file_name, initial_line[i], final_line[i])
        else:
            data = read_chunk(file_name, initial_line[i], final_line[i])
            #The last step of the nth chunk is equal to the first step in 
            # the next one, so we need to delete it
            total=np.delete(total, -1, axis=0) 
            total=np.vstack([total, data])


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

def thermo_analyser(file_name, minimum = 0.3):
    """
    Args:
        file_name
        min Number or percentage (between 0-1) of samples to be discarded
    Returns:
        thermo_data a dataframe that has the columns as the average, etc and the rows for each property
        check the columns with thermo_data.index.values, thermo_data.column.values 
    """

    data, header = data_extract(file_name)
    save_file(data, header)
    stat.fast_averager("Parameters.dat", minimum, "thermo.dat")

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
    lines=f.readlines()
    separator = "="
    array = []
    for l in lines:
        l = l.replace(separator,"\t")
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

    
    
    