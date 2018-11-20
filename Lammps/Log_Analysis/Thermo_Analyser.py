"""
This script reads the log file and gets the
"""
import numpy as np
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


def read_chunk(i_line,f_line):
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



def data_extract(file_name,discard):
    """
    This function extracts the data from the timesteps after discarding the defined amount
    """

    out,err=cf.bash_command("""grep -n "Per MPI" %s| awk -F":" '{print $1}'"""%file_name )
    initial_line=np.array(out.split(),dtype=int)+1

    out2,err=cf.bash_command("""grep -n "Loop time" %s| awk -F":" '{print $1}'"""%file_name)
    final_line=np.array(out2.split(),dtype=int)-1
    #Check if there was minimization
    out3,err=cf.bash_command("""grep -n "Minimization stats" %s"""%file_name)

    if out3:
        initial_line=np.delete(initial_line,0)
        final_line=np.delete(final_line,0)

    """Checking if the simulation did not finish"""
    if np.size(initial_line)>np.size(final_line):
        print "\nThe log file is incomplete!, lets analyse it in anycase \n"
        out3,err=cf.bash_command("""wc -l %s |awk '{print $1}'"""%file_name)
        line=int(out3.split()[0])
        final_line=np.append(final_line,line)


    header=linecache.getline(file_name, initial_line[0]).strip('\n').split()
    header_string=" ".join(header)
    number_chunks=np.size(initial_line)

    for i in xrange(number_chunks):
        if i==0:
            total=read_chunk(initial_line[i],final_line[i])
        else:
            data=read_chunk(initial_line[i],final_line[i])
            total=np.delete(total,-1,axis=0) #As there are repeated timesteps, I choose to keep the last version of the variables
            total=np.vstack([total,data])

    total=total[discard:]

    return total, header_string





if __name__=="__main__":
    parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
    parser.add_argument('FileName', metavar='InputFile',help='Input filename',type=lambda x: cf.is_valid_file(parser, x))
    parser.add_argument('--discard', help='Number of initial steps to discard', default=0, type=int)
    args = parser.parse_args()
    file_name=args.FileName

    data,header=data_extract(file_name,args.discard)
    save_file(data,header)
    stat.fast_averager("Parameters.dat",0,"thermo.dat")
