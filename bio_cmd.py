#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

# ************************************************************************ 
# DESCRIPTION:
# ........................................................................
# Command line version of the multiple charge BioFET-SIM program.
# Run by entering the below statement on the command line:
# python bio_run.py <INPUT.PY>
# ------------------------------------------------------------------------ 

# ************************************************************************ 
# PREPARATION
# ........................................................................
from bio_lib import *
import bio_com
import pickle
#import scipy
# ------------------------------------------------------------------------ 

# ************************************************************************ 
# CHECK ARGUMENTS AND SETUP
# ........................................................................
def get_info(input_file_name):
    dat = open(input_file_name, 'rb')
    val = pickle.load(dat)
    print val['comment']
    print 'Parameters:'
    params = val.keys()
    params.sort()
    info = ""
    for k in params:
        if k=='rho' or k=='comment':
            pass
        else:
            info += "%s%s\n" % (k.ljust(20), str(val[k]).rjust(10))
    return info

def check_args():
    # Check the '--info' and '--calc' arguments. Check the last argument is BFS input file.
    if len(sys.argv) == 3 and sys.argv[1]=='--calc':
        return True 
    if len(sys.argv) == 2 and (sys.argv[1]=='--set'):
        print "Adjust parameter using: <param> <val>"
        print "Example: $ python bio_run.py --set L_d 3.0 input.bfs"
        print "Adjustable parameters:"
        print "L_d, L_tf, eps_1, eps_2, eps_2, eps_3, lay_bf, lay_ox, mu, n_0"
        return False
    if len(sys.argv) == 5 and sys.argv[1]=='--set':
        dat = open(sys.argv[-1], 'rb')
        val = pickle.load(dat)
        target = os.path.splitext(sys.argv[-1])[0]
        set_par = sys.argv[2]
        new_val = float(sys.argv[3])
        val[set_par] = new_val
        adjust(target, val)
    # Generic error message.
    else:
        print "BioFET-SIM usage:"
        print "$ python bio_run.py --calc <input.bfs>"
        print "or"
        print "$ python bio_run.py --set <param> <val> <input.bfs>"
        return False
# ------------------------------------------------------------------------ 

# ************************************************************************ 
# SET: Adjust a parameter.
# ........................................................................
def adjust(target, val): 
    # Open, adjust and close the BFS input file.
    dat=open(target + '.bfs', 'wb')
    pickle.dump(val, dat)
    print "Parameter adjusted."
# ------------------------------------------------------------------------ 

# ************************************************************************ 
# COMPUTE
# - x_lbl: independent variable.
# ........................................................................
def compute(): 
    # --------------------------------------------------------------------
    # PARAMETERS
    # To adjust, remove first '#' and add value after '=' sign.
    # ....................................................................
    sensitivity = str(bio_com.compute(val['rho'], val['nw_len'], val['nw_rad'], val['lay_ox'], val['L_d'],
                                      val['L_tf'], val['lay_bf'], val['eps_1'], val['eps_2'],
                                      val['eps_3'], val['n_0'], val['nw_type'], val['num_prot']))
    return sensitivity
# ------------------------------------------------------------------------ 

# ************************************************************************ 
# LAUNCH:
# ........................................................................
if __name__ == '__main__': 
    if check_args():
        # Adjust a parameter.
        if sys.argv[1] == '--set':
            print "Adjust parameter."
        # Start calculation.
        if sys.argv[1] == '--calc':
            input_file_name = sys.argv[2]
            dat = open(input_file_name, 'rb')
            val = pickle.load(dat)
            dat.close() 
            print get_info(sys.argv[2])
            print "Base Conductance [nS]:", G0(val['nw_len'], val['nw_rad'], val['n_0'], val['mu'])
            print "Sensitivity:          ", compute()
# ------------------------------------------------------------------------ 
