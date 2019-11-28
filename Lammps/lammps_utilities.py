import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.Log_Analysis.Thermo_Analyser as ta



def read_box_limits(log_name):
    """
    TODO make a lammps utilities function Copied from first_n_analysis
    Reads the box limits from log.lammps
    ONLY required for .xyz not for .dump
    Args:
        None: log_name name of the log file
    returns:
        volume
        limits

    """
    if not os.path.exists(log_name):
        print ("The log file specified does not exist")
        
        sys.exit("The log file specified does not exist")
        
        
    out,err=cf.bash_command("""grep -n "orthogonal box" %s | awk -F":" '{print $1}' """%log_name)
    line=int(out.split()[0])
    limits=linecache.getline(log_name, line)
    limits=re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", limits)
    limits=np.array(np.reshape(limits,(2,3)),dtype=float) #To have the box as with columns [x_i_min,x_i_max]
    volume=(limits[1,0]-limits[0,0])*(limits[1,1]-limits[0,1])*(limits[1,2]-limits[0,2])

    return volume,limits

# TODO generalise to get all the dimension
def read_region_height(domain_name, geom_file = "in.geom"):
    """
    Reads the bulk limits from in.geom
    Args:
        domain_name: is the name of the region to be analysed
        geom_file: Is the name of the file that contains the region definition
    returns:
        height

    """
    if not os.path.exists(geom_file):
        print ("The geometry file specified does not exist")
        
        sys.exit("The geometry file specified does not exist")
        
        
    out,err=cf.bash_command("""grep -n "%s" %s | awk -F":" '{print $1}' """%(domain_name,geom_file))
    line=int(out.split()[0])
    limits=linecache.getline(geom_file, line)
    limits=np.array(re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", limits),dtype = float)
    height = limits[1]-limits[0]

    return height, limits
# Create a class region and then  simulation box that is a inherits the rest with properties xmax,xmin, lx, vol, etc


def read_value_from(input_file,key_phrase):
    """
    Extracts the numberical value from a line well specified by the key_phrase
    Args:
        input_file: name of the file to look for the value log_file, in.*, etc
        key_phrase has to be very specific such that is not repeated in the file, ex "atoms in group gSolv"
    returns:
        value requested

    """
    if not os.path.exists(input_file):
        print ("The file %s specified does not exist"%input_file)
        
        sys.exit("The file %s specified does not exist"%input_file)
        
        
    out,err=cf.bash_command("""grep "%s" %s"""%(key_phrase,input_file,))
    value = cf.extract_digits(str(out))

    return value