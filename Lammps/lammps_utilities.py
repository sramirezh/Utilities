import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.Thermo_Analyser as ta



def read_box_limits(log_name, entry_num = 0):
    """
    TODO make a lammps utilities function Copied from first_n_analysis
    Reads the box limits from log.lammps
    ONLY required for .xyz not for .dump
    Args:
        None: log_name name of the log file
        entry_num = 0, if the initial box, -1 to get the last box, for instance
        after barostating
    returns:
        volume
        limits

    TODO: delete this function once  Simulation class is implemented
    """
    if not os.path.exists(log_name):
        print ("The log file specified does not exist")
        
        sys.exit("The log file specified does not exist")
        
        
    out,err = cf.bash_command("""grep -n -a "orthogonal box" %s | awk -F":" '{print $1}' """%log_name)
    line = int(out.split()[entry_num])
    limits = linecache.getline(log_name, line)
    limits = re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", limits)
    limits = np.array(np.reshape(limits,(2,3)),dtype=float) #To have the box as with columns [x_i_min,x_i_max]
    volume = (limits[1,0]-limits[0,0])*(limits[1,1]-limits[0,1])*(limits[1,2]-limits[0,2])

    return volume, limits

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


def read_value_from(input_file, key_phrase):
    """
    Extracts the numerical value from a line well specified by the key_phrase
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


def solid_surface(data, atom_type):
    """
    Computes the limits of the solid surface
    Args:
        data with the atom_type, posx,posy,posz
        atom_type atom type in the trajectory file of the solid surface.
    Returns:
        Characteristics of the solid surface
        Writes a file Zshift.dat with the maximun height to be used by other codes.

    Writes the zshift
    """

    #Getting the maximum position of the surface.
    indexes=np.where(data[:,0]==atom_type)[0]
    Maxz=np.max(data[indexes,3])
    Minz=np.min(data[indexes,3])

    print("The maximum height of the solid surface is %lf" %Maxz)
    print("The minimum height of the solid surface is %lf" %Minz)
    print("The height of the solid surface is %lf" %(Maxz-Minz))

    #Writing the Zshift
    f=open("Zshift.dat",'w')
    f.writelines("%lf \n" %Maxz)
    f.close

def get_variables(file_name):
    """
    Creates a dictionary from all the numerical variables from the given input
    file
    Args:
        file_name: Name of the input file to extract values.

    Returns:
        var_dir: dictionary containing the variable name as the keys and the 
        values. 
        You can use dictionary.keys() and dictionary.values on the output
    TODO: extract also values from symbolic variables depending on other 
          variables.
    """
    f = open(file_name)
    names =[]
    values =[]
    for line in f.readlines():
        if "variable" in line:
            temp = line.split()
            var_val = temp[3]
            if '$' not in var_val:
                values.append(cf.extract_digits(var_val)[0])
                names.append(temp[1])
    f.close()
    var_dir = dict(zip(names, values))
        
    return var_dir


# Once Simulation gets all the variables, this class could be unified
class SimulationType(object):
    """Class to define the type of simulation analysed, for example N2 or
    Octane
    """
    def __init__(self, name):
        """
        Args:
            name: to identify the simulation type, eg. N2
        """
        self.name = name    
        
    def print_params(self, logger):
        """
        Args:
            logger is an object of the class cf.log
        """
        logger.info("\nUsing the parameters from %s\n"%self.name)
        dictionary = vars(self)
        for key in dictionary.keys():
            if key != 'name':
                logger.info("Using %s as %s"%(key, dictionary[key]))
        logger.info("\n")


# Todo could create a method to get all the variables and their values and then create a dictionary
class Simulation(object):
    """
    TODO
    There are two classes that could be useful and probably merged later
    -system in widom_pylammps
    -simulation from qsub simulation results

    Atributes:
        log_file = name of the log file
        thermo_ave = average of the parameters in the log file [DATAFRAME]
        thermo_data = data from the log file [DATAFRAME]
        volume = intial volume from the box
        limits = initial limits of the box

    TODO: Disbale the pritns when calling thermo_analyser

    """
    def __init__(self, log_file):
        self.log_file = log_file
        self._initial_check()
        self._get_thermo_data()
        self._read_box_limits()
    
    def _initial_check(self):
        if not os.path.exists(self.log_file):
            print ("The input file %s does not exist"%self.log_file)
            sys.exit("The input file %s does not exist"%self.log_file)

    def _get_thermo_data(self):
        self.thermo_ave = ta.thermo_analyser(self.log_file).astype('float')
        self.thermo_data = cf.read_data_file("Parameters.dat")

    def _read_box_limits(self):
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
        if not os.path.exists(self.log_file):
            print("The log file specified does not exist")
            
            sys.exit("The log file specified does not exist")
            
        out, err = cf.bash_command("""grep -n -a "orthogonal box" %s | awk -F":" '{print $1}' """%self.log_file)
        line = int(out.split()[0])
        limits = linecache.getline(self.log_file, line)
        limits = re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", limits)
        limits = np.array(np.reshape(limits, (2,3)),dtype=float) # To have the box as with columns [x_i_min,x_i_max]
        self.volume = (limits[1, 0] - limits[0, 0]) * (limits[1, 1] - limits[0, 1]) * (limits[1, 2] - limits[0, 2])
        self.limits = limits

    def is_2d(self):
        out, err = cf.bash_command("""grep -n "enforce2d" %s | awk -F":" '{print $1}' """%self.log_file)

        if out:
            return True
        else:
            return False