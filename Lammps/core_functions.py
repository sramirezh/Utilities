
from subprocess import Popen, PIPE
from shlex import split
import argparse
import os
import pandas as pd
import re
import numpy as np
import sys
from cycler import cycler
import warnings
import pickle as pickle
import matplotlib.pyplot as plt
import logging
from distutils.spawn import find_executable

warnings.filterwarnings("ignore")


class log(logging.Logger):
    """
    Class extending the logger object from the module log
    Class to print everything to a file and to screen.
    Instead of using print, use log.info("text") or read the documentation
    """
    def __init__(self, path, cwd):
        """
        file_name is the name of the python file currently running
        path: given by __file__
        """
        logging.Logger.__init__(self, "")  # run parent __init__ class\
        self.setLevel(logging.DEBUG)
        self.path = os.path.dirname(path)
        self.file_name = path.split('/')[-1]
        self.cwd = cwd
        log_name = self.file_name.split('.')[0]
        self.log_file = "%s.log"%log_name       
        
        self._set_console_handler()
        self._set_file_handler()

        self.print_basic_info()

    def _set_console_handler(self):
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)  # or any other level
        self.addHandler(ch)
        
    def _set_file_handler(self):
        
        fh = logging.FileHandler(self.log_file, 'w')
        fh.setLevel(logging.INFO)  # or any level you want
        self.addHandler(fh)

    def _get_git_hash(self):
        os.chdir(self.path)
        hash_pre, err = bash_command("""git rev-parse HEAD""")
        hash_post = hash_pre.strip().decode("utf-8")
        self.hash = hash_post
        os.chdir(self.cwd)
    
    def print_basic_info(self):
        import datetime
        self._get_git_hash()
        self.info("\n--------------------------------------------------------")
        self.info(datetime.datetime.now())
        self.info("Running %s"%self.file_name)
        self.info("The git hash is %s" % self.hash)
        self.info("Running inside %s" % self.cwd)
        self.info("--------------------------------------------------------\n")
        


def bash_command(cmd):
    """
    function that evaluates simple bash commands with pipes,

    Input:
        cmd is a string of the command just as you write it on the shell but 
        inside 3" in case you have
        several quoted text inside
    Returns: two elements
        out The output
        err the possible errors
    """
    cmd = str(cmd)
    chain = cmd.split("|")
    n_pipes = len(chain)

    for i in range(n_pipes):
        if i == 0:
            p = Popen(split(chain[0]), stdout=PIPE)
        else:
            p = Popen(split(chain[i]), stdin=p.stdout, stdout=PIPE)

    return p.communicate()


def is_valid_file(parser, arg):
    """
    Checks if a file exists and returns the proper output to the parser handler.
    """
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def parameter_finder(target_list, search_list, msgflag=False, exact=False):
    """
    Finds all the elements in the search list in the target list
    the search is non case sensitive

    Args:
        Target list is the list where I want to look for the items in the 
        search
         list
        Search list is the list of items I want to look for
        msg is True to tell that there were occurrences or not.
        exact True if the search requires the exact string, False only if it 
        requires to contain the search list, example "vx" is in "vx_sol" 
    Returns:
        indexes with the positions of the search_list found in the target_list
    """
    target_list = [x.lower() for x in target_list]

    indexes = []


    if isinstance(search_list, str):
        cont = 0
        search_list = search_list.lower()
        for t in target_list:
            if exact == False and search_list in t:
                indexes.append(cont)
            elif exact == True and search_list == t:
                indexes.append(cont)
            cont += 1
    if isinstance(search_list, list):
        search_list = [x.lower() for x in search_list]

        for s in search_list:
            s = str(s)
            for cont, t in enumerate(target_list):
                if exact == False and s in t:
                    print((s, t))
                    indexes.append(cont)
                elif exact == True and s == t:
                    print((s, t))
                    indexes.append(cont)

    
    if msgflag == True:
        length = len(indexes)
        if length > 1: print("There were several ocurrences")
        if length == 0: print("No ocurrences found")

    return indexes


def read_data_file(input_file):
    """
    Reads a data file either with a header or not.
    It assumes that the header is commented with "#" and that the it last line contains the name of the variables
    Args:
        input_file file name

    Returns: A panda data frame, the column names can be obtained by data.columns.values and the numeric parameters with  data.values
    """
    header_lines = 0
    last_pound_pos =- 1
    with open(input_file, 'r') as data_file:
        while(data_file.read(1) == '#'):
            last_pound_pos = data_file.tell()
            header = data_file.readline()
            header_lines += 1

        #Read the next lines
        data_1 = data_file.readline().split()
        data_2 = data_file.readline().split()
        data_file.seek(last_pound_pos+1) #Goes back to the last line of the header

        if header_lines == 0:
            data=pd.read_csv(data_file, sep=" ", header=None).dropna(axis=1, how='all')

        else:
            if len(data_1) != len(data_2): #If there is a line containing the number of particles,
                data_file.readline()
            data_file.readline()

            try:
                data = pd.read_csv(data_file, sep=" ", header=None).dropna(axis=1, how='all')
                data.columns = header.split()
            except:
                raise Exception("The input file '%s' is corrupted, usually the problem is because "\
                                "there is an end of a line that has an additional space" %input_file)

    return data

def extract_digits(strings, sort=True):
    """
    input:
        strings: Array or single string
    Returns:
         An array of all the digits in a string
         if it is an array of strings, extracts the numbers and returns them sorted
         indexes is contains the indexes to run the input strings such that they are organised

        later the strings can be organised as 
        a = cf.extract_digits(files)
        index_files = a[1]
        files = [files[i] for i in index_files]

    """
    if isinstance(strings, str):
        numbers = re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", strings)
        output = ([num.strip('.') for num in numbers])
        return np.array(output, dtype=float)


    if isinstance(strings, list):
        output = []
        for element in strings:
            numbers = re.findall(r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", element)
            output.append([num.strip('.') for num in numbers])

        if np.shape(output)[1] == 1:
            """
            If there are several parameters it is difficult to reshape
            """
            if sort == True:
                unsorted = np.array(output,dtype=float).reshape((len(output)))
                columns = np.column_stack(([unsorted,np.arange(0,len(unsorted))]))
                columns = columns[columns[:,0].argsort()]
                index_sorted = columns[:,1].astype(int)



                return np.array(unsorted[index_sorted] ,dtype=float),index_sorted
            else:
                output = np.array(output,dtype=float).reshape((len(output)))
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



def set_plot_appearance(presentation_type = False):

    """
    Defines the appearence of the plot

    use rcParams.keys() to see the available parameters
    """
    if presentation_type == False: 
        axis_font=24
        tick_font=20
        legend_font=12


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
        if find_executable('latex'):
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

    if presentation_type == True: 
        axis_font= 30
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
        if find_executable('latex'):
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
        plt.rcParams['legend.loc'] = 'upper left'
        plt.rcParams['legend.labelspacing'] = 0.5
        plt.rcParams['legend.borderpad'] = 0.4
        plt.rcParams['legend.frameon'] = True
        plt.rcParams['legend.fancybox'] = False
        plt.rcParams['legend.edgecolor'] = 'k'
        # Fonts and symbols


        # Errorbar plots
        plt.rcParams['errorbar.capsize'] = 4

def integrate(x, y, xmin, xmax):
    """
    Integrate the data in x and y from xmin to xmax
    
    Parameters
    ----------
    
        x : are the independent variables
        y : the integrand could be a unumpy object
        xmin : the lower limit of the integral
        xmax : the upper limit of the integral
    """
    MinIndex = np.min(np.where(x >= xmin))
    MaxIndex = np.max(np.where(x <= xmax))
    integral = np.trapz(y[MinIndex:MaxIndex], x[MinIndex:MaxIndex])

    return integral






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
    
    print(("Modified the file %s"%(file_name.split('/')[-1])))


def load_instance(file_name):
    """
    Loads the data structure to be used later
    """
    file1 = open(file_name, 'rb')
    instance = pickle.load(file1)
    file1.close()

    return instance

def save_instance(instance, file_name):
    """
    Saves the instance 
    args:
        instance: is the objet that is going to be saved, if it is for a plot, save the ax 
        file_name name without extension
    Returns: Nothing
    """
    afile = open(r'%s.pkl'%file_name, 'wb')
    pickle.dump(instance, afile)
    afile.close()   


def str2list(input):
    """
    Check if the input is a string, if so, converts it to a list

    Returns:
        a list
    """
    if isinstance(input,str):
        return [input]
    
    else:
        return input



def beware_msg(msg):
    """
    Prints a message that can be easily spotted
    """
    print ("\n\n\n************************************************************")
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n")
    print (msg)
    print ("\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print ("************************************************************\n\n\n")


def group_consecutive(data):
    """
    Groups consecutive indexes in an array an return a list with all the groups
    Args:
        data 1D data array
    Returns:
        results a list with all the consecutive indexes grouped
    """
    from itertools import groupby
    from operator import itemgetter
    results=[]
    for k, g in groupby(enumerate(data), lambda i_x: i_x[0]-i_x[1]):
        results.append(list(map(itemgetter(1), g)))
        
         
    return results
         

def plateau_finder(data,tol=0.0003):
    """
    Function that finds the plateaus of a distribution y, that has x spacing constant
    Args:
        data 1D data that has delta in the other dimension constant
        tol tolerance for the variance
    """
    from scipy.ndimage.filters import generic_filter
    tol=0.0003
    filt_data = generic_filter(data, np.std, size=3)
    plat_index=np.where(filt_data<(np.min(filt_data)+tol))[0]
    
    plateaus=group_consecutive(plat_index)
    
    return plateaus


def plot_zoom(ax, xlims):
    """
    Returns the zoomed axis such that it focuses on the limits given in the x
    axis
    
    Args:
        ax: Pyplot axis
        xlims: array with the limits in x, eg. [0,10]
        
    Returns:
        ax after the zooming
    """
    
    xmin, xmax, ymin, ymax = get_y_lims(ax, xlims)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    return ax


#   TODO generalise this
def get_y_lims(ax, xlims):
    """
    
    Returns the bounds for the given limits of the plot
    
    Assumes that all the plots share the same x positions
    
    Args:
        ax: axis with all the lines that want to be plotted
        xlims: the zooming region
    
    Returns:
        xmin: x coordinate in ax that is equal or more than the lower bound
        xmax: x coordinate in ax that is equal or less than the upper bound
        ymin: minimum value from all the lines in ax that is whithin the interval
        ymax: minimum value from all the lines in ax that is whithin the interval
        
    """
    # Assuming that all objects have the same x coordinates
    x = ax.lines[0].get_data()[0]
    
    min_index = np.min(np.where(x >= xlims[0]))
    max_index = np.max(np.where(x <= xlims[1]))
    xmax = x[max_index]
    xmin = x[min_index]
    
    ymax_array = []
    ymin_array = []
    
    for function in ax.lines:
        y = function.get_data()[1]
        
        ymin_array.append(np.min(y[min_index:max_index]))
        ymax_array.append(np.max(y[min_index:max_index]))
    
    ymax = max(ymax_array)
    ymin = min(ymin_array)

    return xmin, xmax, ymin, ymax