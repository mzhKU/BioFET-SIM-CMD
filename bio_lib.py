#!/usr/bin/python

print "Content-type: text/plain\n"

# ************************************************************************
# BioFET-SIM Library Functions.
# ........................................................................  

import os
import sys
#os.environ['PYTHONPATH'] = '/Library/Frameworks/Python.framework/Versions/7.2/lib/python2.7/site-packages'
#sys.path.append('/Library/Frameworks/Python.framework/Versions/7.2/lib/python2.7/site-packages/')
os.environ['PYTHONPATH'] = '/usr/bin/python'

from scipy import log, sqrt, exp, pi, power
from scipy.special import i0, i1, k0, k1
from subprocess import Popen, PIPE
from string import Template
import itertools
import datetime
import os.path
#import bio_tst
import string
import urllib
#import numpy
import scipy
import cgi

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---------- Do not edit this module without appropriate care. -----------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

# ************************************************************************
# CONTENT:
# - constants.
# - special functions.
# - geometrical calculations; pKa calculation.
# - File generation; parsing; fixing; NW surface; PDB download.
# - HTML, Jmol and GNUPlot
# - BioFET-SIM Parameters
# - Show functions
# ------------------------------------------------------------------------ 

# ************************************************************************
# NOTES
# ........................................................................  
# - prepareMulti.html could be exported to a Django template.
# ------------------------------------------------------------------------ 

# ************************************************************************
# EDITS:
# ........................................................................  
# 16.02.2012: Distinction between nanowire and nanoribbon.
# 17.02.2012: Discarding nanoribbon in def Gamma_li
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: CONSTANTS.
# ........................................................................  
q_elem  = 1.602176487E-19 # Coulomb.
h_bar   = 1.054571628E-34 # Reduced Plank's constant.
eps_0   = 8.854187817e-12 # Vacuum permittivity.
el_mass = 9.10938215e-31  # Electron mass at rest.
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: GLOBAL VARIABLES IN APPLICATION SCOPE.
# ........................................................................  

# KU machine: Directories
#vmd_base_path = '/opt/vmd_package/Contents/vmd/'
#pdb_base_path = '/Users/mzhKU_work/Sites/ku_prototype/bfs_pdb/'
#results_path  = '/Users/mzhKU_work/Sites/ku_prototype/bfs_res/'
#jmol_path     = '/Users/mzhKU_work/Sites/ku_prototype/bfs_spt/'

# KU machine: Applications
#vmd_cmd_path              = '/Users/mzhKU_work/software/vmd_package/Contents/vmd/vmd_MACOSXX86'
#external_python_base_path = '/Library/Frameworks/Python.framework/Versions/Current/bin/python'
#propka_path               = '/Users/mzhKU_work/software/propka3/propka.py'
#pdb2pqr_base_path         = '/opt/pdb2pqr/pdb2pqr.py'

# PROPKA: Directories
vmd_base_path = '/var/www/vmd/bin_vmd/vmd'
pdb_base_path = '/var/www/propka/biofet-sim/bfs_pdb/'
results_path  = '/var/www/propka/biofet-sim/bfs_res/'
# PROPKA: Applications
external_python_base_path = '/usr/bin/python'
vmd_cmd_path              = '/var/www/vmd/bin_vmd/vmd'
pdb2pqr_base_path         = '/opt/pdb2pqr/pdb2pqr.py'
python3_path              = '/usr/local/bin/python3'
propka_path               = '/opt/propka30/propka.py'
gnuplot_exe               = '/usr/bin/gnuplot'
convert_exe               = '/usr/bin/convert'
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: SPECIAL FUNCTIONS.
# ........................................................................  
# Dimensionless function, quantifying actual sensitivity of the wire in
# the presence of Debye screening in the electrolyte and a finite
# Thomas-Fermi screening in the nanowire (Sorensen et al.).
def Gamma(nw_rad, lay_ox, L_d, L_tf, eps_1, eps_2, eps_3): 
    fact1 = (nw_rad + lay_ox)/L_d
    fact2 = nw_rad/L_tf
    fact3 = fact1**(-1)
    fact4 = (nw_rad + lay_ox)/nw_rad 
    num    = eps_1*k0(fact1)*(L_d/L_tf)*i1(fact2)
    denom1 = k0(fact1)*fact3 
    denom2 = log(fact4)*k1(fact1)*(eps_3/eps_2)
    denom3 = (denom1 + denom2)*eps_1*fact2*i1(fact2) 
    denom  = denom3 + eps_3*k1(fact1)*i0(fact2) 
    gamma = num/denom 
    return gamma

def G0(nw_len, nw_rad, n_0, mu):
    """Conductance in the absence of surface charge G0."""
    return pi*nw_rad**2*q_elem*n_0*mu/nw_len

# Dimensionless function to quantify the effect of sigma_bi, [nm].
def Gamma_li(nw_rad, li, L_d): #, mode):
    ## Eq. (2), Vacic et al., 2011, JACS 
    #if mode == 'Nanowire':
    #    dist_i    = nw_rad/(nw_rad + li) 
    #    gamma_li = 2*dist_i*(1 + sqrt(dist_i)*exp(li/L_d))**(-1)
    ## Eq. (3), Vacic et al., 2011, JACS 
    #else:
    #    gamma_li = 2*(1+exp(li/L_d))**(-1)
    # Eq. (2), Vacic et al., 2011, JACS 
    dist_i    = nw_rad/(nw_rad + li) 
    gamma_li = 2*dist_i*(1 + sqrt(dist_i)*exp(li/L_d))**(-1)
    return gamma_li

# Compute surface charge density.
def Sigma_bi(nw_rad, li_tot, nw_len, num_prot, q_i): 
    # Converts to meter and Coulomb.
    nw_rad = nw_rad * 1E-9 
    li_tot = li_tot * 1E-9
    nw_len = nw_len * 1E-9 
    q_i = float(q_i)*q_elem 
    sigma_bi = num_prot*q_i/(2*pi*(nw_rad + li_tot)*nw_len)
    return sigma_bi

# Get number of ionizable sites in the protein.
def get_m(of):
    return len(of)
# ------------------------------------------------------------------------

# ************************************************************************
# SECTION: GEOMETRICAL CALCULATIONS; pKa calculation:
# - Positioning of protein on surface
# - Averaging of amino acid side chains to locate charge
# - Reorientation to symmetry axes
# - Get protein bounding box dimensions
# ........................................................................  
# Get coordinate of the charge with lowest z-coordinate value.
# The offset is in the opposite direction. The returned value is in [Ang].
def get_z_offset(rho):
    z_values = [] 
    for ri in rho:
        z_values.append(float(ri[2]))
    # Generate an error to show mod.rho in cgi traceback.
    #z + "2"
    return min(z_values)*(-1)

def glu(atms): 
    """Return OE1/OE2 average.
    """
    oe1 =  [ float(atms[1][0]),
             float(atms[1][1]),
             float(atms[1][2]) ]
    oe2 =  [ float(atms[2][0]),
             float(atms[2][1]),
             float(atms[2][2]) ]
    return [ atms[0], (oe1[0]+oe2[0])/2.0,
                      (oe1[1]+oe2[1])/2.0, (oe1[2]+oe2[2])/2.0 ] 

def cys(atms):
    """Return SG location.
    """
    return [ atms[0],
             float(atms[1][0]),
             float(atms[1][1]),
             float(atms[1][2]) ]

def his(atms):
    """Return CG/ND1/CD2/CE1/NE2 average.
    """
    cg  =  [ float(atms[1][0]),
             float(atms[1][1]),
             float(atms[1][2]) ]
    nd1 =  [ float(atms[2][0]),
             float(atms[2][1]),
             float(atms[2][2]) ]
    cd2 =  [ float(atms[3][0]),
             float(atms[3][1]),
             float(atms[3][2]) ]
    ce1 =  [ float(atms[4][0]),
             float(atms[4][1]),
             float(atms[4][2]) ]
    ne2 =  [ float(atms[5][0]),
             float(atms[5][1]),
             float(atms[5][2]) ]
    return [ atms[0], (cg[0]+nd1[0]+cd2[0]+ce1[0]+ne2[0])/5.0,
                      (cg[1]+nd1[1]+cd2[1]+ce1[1]+ne2[1])/5.0,
                      (cg[2]+nd1[2]+cd2[2]+ce1[2]+ne2[2])/5.0 ] 

def asp(atms):
    """Return OD1/OD2 average.
    """
    od1 =  [ float(atms[1][0]),
             float(atms[1][1]),
             float(atms[1][2]) ]
    od2 =  [ float(atms[2][0]),
             float(atms[2][1]),
             float(atms[2][2]) ]
    return [ atms[0], (od1[0]+od2[0])/2.0,
                      (od1[1]+od2[1])/2.0, (od1[2]+od2[2])/2.0 ] 

def tyr(atms):
    """Return OH location.
    """
    return [ atms[0], float(atms[1][0]),
            float(atms[1][1]), float(atms[1][2]) ]

def arg(atms):
    """Return CZ location.
    """
    return [ atms[0], float(atms[1][0]),
            float(atms[1][1]), float(atms[1][2]) ]

def lys(atms):
    """Return NZ location.
    """
    return [ atms[0], float(atms[1][0]),
            float(atms[1][1]), float(atms[1][2]) ]

def oxt(atms):
    """Return OXT location.
    """ 
    return [ atms[0], float(atms[1][0]),
            float(atms[1][1]), float(atms[1][2]) ]

def np(atms):
    """Return N+ location.
    """ 
    return [ atms[0], float(atms[1][0]),
            float(atms[1][1]), float(atms[1][2]) ]

# Coordinate averages.
def av_res(ion_res):
    """Switch function providing the access to the function to call
    for the specific amino acid.
    """
    # 'ion_res' is "['LYS 2 A NZ 3.485 3.366 1.894']". From this, all
    # other properties can be evaluated.
    av_functions = {'ASP':asp, 'GLU':glu,
                    'HIS':his, 'CYS':cys,
                    'TYR':tyr, 'LYS':lys,
                    'ARG':arg}
    av_func = av_functions[ion_res[0].split()[0]]
    av_atms = [" ".join(ion_res[0].split()[0:3])]
    for i in ion_res:
        av_atms.append(i.split()[4:]) 
    return av_func(av_atms)

def av_trm(ion_trm):
    """Return a terminus in the same format as the ionizable residues.
    """
    trm_lbl = [" ".join(ion_trm.split()[:4])]
    coords = [float(i) for i in ion_trm.split()[4:]]
    # Combine two lists element wise.
    return list(itertools.chain(trm_lbl, coords))

# Return the VMD script for controlling the reorientation.
def get_reo_src(com_xyz, target):
    # <<PATH>>
    reo_src =   'mol load pdb %s-fix.pdb\n' % (pdb_base_path + target)\
              + 'set sel [atomselect 0 \"protein\"]\n'\
              + 'atomselect0 moveby {%.6f  %.6f  %.6f}\n' %\
                 (com_xyz[0], com_xyz[1], com_xyz[2])\
              + 'package require Orient\n'\
              + 'namespace import Orient::orient\n'\
              + 'set sel [atomselect top \"protein\"]\n'\
              + 'set I [draw principalaxes $sel]\n'\
              + 'set A [orient $sel [lindex $I 2] {0 0 1}]\n'\
              + '$sel move $A\n'\
              + 'set I [draw principalaxes $sel]\n'\
              + 'set A [orient $sel [lindex $I 1] {0 1 0}]\n'\
              + '$sel move $A\n'\
              + 'set I [draw principalaxes $sel]\n'\
              + '$sel writepdb %s-reo.pdb' % (pdb_base_path + target)
    return reo_src

def get_box_dimensions(pqr):
    """Return bounding box dimensions"""
    X = []
    Y = []
    Z = [] 
    for atm in pqr.split('\n'):
        if len(atm.split()) > 0:
            if atm.split()[0] == 'ATOM':
                x = atm[30:38].strip()
                y = atm[38:46].strip()
                z = atm[46:54].strip()
        X.append(float(x))
        Y.append(float(y))
        Z.append(float(z)) 
    #for atm in pqr[1:]:
    #    if len(atm.split()) > 0:
    #        x = atm.split()[1]
    #        y = atm.split()[2]
    #        z = atm.split()[3]
    #    X.append(float(x))
    #    Y.append(float(y))
    #    Z.append(float(z)) 
    # Scan all residues.
    x_min = min(X)
    x_max = max(X)
    y_min = min(Y)
    y_max = max(Y)
    z_min = min(Z)
    z_max = max(Z) 
    return x_max-x_min, y_max-y_min, z_max-z_min
# ........................................................................
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: FILE GENERATION; PARSING; FIXING; NW SURFACE
# - PDB file upload
# - PDB file parsing
# - pKa calculation
# - pKa file parsing
# - PQR file writing
# - Computation of charge on residue;
# - Computation of number of proteins on NW
# - Generation of terminals
# - Resetting chains
# - Generate NW surface
# - Download pdb
# - Get Q_tot(pH)
# - Generate offline BFS input (parameters and coordinates)
# ........................................................................
def save_uploaded_file(form):
    # Upload handle.
    #upload_dir = "./uploaded"
    file_dat  = form["filename"].file
    file_name = form["filename"].filename
    #file_val = open(os.path.join(upload_dir, file_name), 'w')
    file_val = open('./bfs_pdb/' + file_name, 'w')
    file_val.write('./bfs_pdb/' + file_dat.read()) 
    file_val.close() 
    # Calculation setup.
    target = os.path.splitext(file_name)[0]
    pH     = float(form["pH"].value)
    return target, pH

# <BFS_CMD_INP>
#def generate_bfs_input(target, params, rho, x_lbl, num_prot, x_min, x_max,
#                       comment, num_qi, file_name):
#def generate_bfs_input(target, params, rho, x_lbl, num_prot,
#                       comment, num_qi, file_name):
def generate_bfs_input(target, params, rho, num_prot, comment, num_qi, file_name):
    import pickle
    f = open(pdb_base_path + file_name, 'wb')
    data = {}
    data['target']   = target
    data['rho']      = rho
    data['num_prot'] = num_prot
    #data['x_lbl']    = x_lbl
    #data['x_min']    = x_min
    #data['x_max']    = x_max
    data['num_qi']   = num_qi
    data['comment']  = comment
    for k in params.keys():
        data[k] = params[k]
    pickle.dump(data, f)
    f.close()

def get_pdb(target):
    pdb_file = open(target + '-reo.pdb', 'r')
    pdb_data = pdb_file.readlines()
    pdb_file.close() 
    return pdb_data

# Return individual PDB labels using PyMOL identifiers
def get_labels(line):
    res_atm = line[12:16].strip()
    res_nam = line[16:20].strip()
    res_chn = line[20:22].strip()
    res_ind = line[22:26].strip()
    res_x   = line[30:38].strip()
    res_y   = line[38:46].strip()
    res_z   = line[46:54].strip()
    return [res_atm, res_nam, res_ind, res_chn, res_x, res_y, res_z]

def rewrite_pdb(target, tmp_pqr):
    # Rewrite PDB file with the coordinates after the move.
    # Use original PDB file as template for residue info details.
    orig = open(pdb_base_path + target + '-reo.pdb', 'r').readlines()
    new  = ''
    cnt  = 0 
    tmp_pqr = tmp_pqr.split('\n')
    for l in orig:
        if l[:4] == 'ATOM':
            new += l[:30] + '%8.3f%8.3f%8.3f'% (float(tmp_pqr[cnt].split()[0]),
                                                float(tmp_pqr[cnt].split()[1]),
                                                float(tmp_pqr[cnt].split()[2])) + l[54:]
            cnt+=1 
    return new

def calc_pKas(target): 
    # <<PATH>>
    #pkap = Popen(['/usr/local/bin/python3',
    #              '/Users/mzh/software/propka30/propka.py',
    #              '%s-reo.pdb' % target], stdout=PIPE,
    #              stderr=PIPE, shell=False)
    # 11.04.2012: Base directory of process is script location.
    os.chdir('./bfs_reo')
    pkap = Popen(['/usr/local/bin/python3',
                  #'/opt/propka30/propka.py',
                  # <<PATH>>
                  #'/home/mzhpropka/software/propka30/propka.py',
                  '/opt/propka30/propka.py',
                  #'--grid', '0.0 14.0 0.1',
                  './%s-reo.pdb' % target], stdout=PIPE,
                  stderr=PIPE, shell=False)
    pka_dat = open('./' + target + '-reo.pka', 'w')
    pka_dat.write(pkap.stdout.read())
    pka_dat.close()
    os.chdir('..')

# Calculate Q_tot(pH)
# 12.04.2012
# - Discarded use
def get_Q_tot(target, pH):
    pka_dat = open(pdb_base_path + target + '-reo.pka', 'r')
    pka_val = pka_dat.readlines() 
    Q_tot_start_line = 0
    for line in enumerate(pka_val):
        if 'Protein charge of folded and unfolded' in line[1]:
            Q_tot_start_line = line[0] + 2
        break
    for Q_tot_line in pka_val[Q_tot_start_line:]:
        if Q_tot_line.split()[0] == 'The':
            break
        if "%2.1f" % float(Q_tot_line.split()[0]) == "%2.1f" % pH: 
            Q_tot = Q_tot_line.split()[2]
    return Q_tot

def calc_Q_tot(pqr):
    Q_tot = 0.0
    for i in pqr.split('\n'):
        Q_tot += round(float(i.split()[-2]), 2)
    return Q_tot
    
# pKa values of the residues.
def get_pKas(pka_dat):
    """The PROPKA output is provided as a stdout file handle.
    This avoids writing a pka file."""
    pka_val = pka_dat.readlines()
    # Locate 'SUMMARY' of pKa values in PROPKA output.
    pka_start_line = 0
    for line in enumerate(pka_val):
        if len(line[1].split()) != 0 and line[1].split()[0] == 'SUMMARY':
            pka_start_line = line[0] + 2 
    # Populate pKa list of residues and terminals.
    pka_tmp = []
    for pka_line in pka_val[pka_start_line:]:
        # Defining residue or terminal identifier 'id'.
        # The 'summary' lines are 5 elements long.
        if len(pka_line.split()) == 5:
            pka_tmp.append(pka_line.split())
    return pka_tmp

def set_pqr(target, av_RQ, pH, pka_dat):
    """PQR file to load in Jmol.
    'target':   Structure label to identify the pKa file.
    'av_RQ[0]':  "['LYS 2 A', 3.484, 3.366, 1.893]"
    'av_RQ[-1]': "['ASP 13 A OXT', 40.159,  16.562, -0.142]"
    In PDB/PQR format data, terminal is labeled 'OXT', in PROPKA
    it is labeled 'C-'.
    FIX:
    - Parse N+ in pKa file.
    - Parse for LIG.
    """ 
    pqr = ""
    # Get pKa values for which coordinates are available.
    # Match the label to fit the PROPKA summary label style.
    # If matched, append.
    pKas = get_pKas(pka_dat)
    cnt = 0 
    # Define a generic label for the charge carrier site,
    # residue or terminus: 'q_i_lbl'.
    for av_rq_i in av_RQ:
        # Amino acid charges. The label is 3 units long.
        # The termini labels are 4 units long.
        if len(av_rq_i[0].split()) == 3:
            q_i_lbl = " ".join(av_rq_i[0].split())
        # Termini charges, adapting to PROPKA terminus format.
        else:
            if av_rq_i[0].split()[-1] == 'N': 
                q_i_lbl = 'N+' + ' ' + " ".join(av_rq_i[0].split()[1:3]) 
            else:
                q_i_lbl = 'C-' + ' ' + " ".join(av_rq_i[0].split()[1:3]) 
        for pka_i in pKas:
            pka_i_lbl = " ".join(pka_i[:3])
            if pka_i_lbl == q_i_lbl:
                pqr += 'ATOM %7d' % int(av_rq_i[0].split()[1])\
                     + '   C ' + av_rq_i[0].split()[0]\
                     + ' ' + av_rq_i[0].split()[2]\
                     + '%16.3f' % av_rq_i[1]\
                     + '%8.3f' % av_rq_i[2]\
                     + '%8.3f' % av_rq_i[3]
                q_i = get_q_i(pka_i_lbl.split()[0], float(pka_i[3]), pH)
                # Prevent empty line at the end of the PQR file.
                if cnt < len(av_RQ)-1:
                    pqr += "%6.3f".rjust(6) % q_i + ' 1.0\n'
                    cnt += 1
                else:
                    # Strictly cannot append '\n' character.
                    pqr += "%6.3f".rjust(6) % q_i + ' 1.0'
    return pqr

def get_q_i(res_nam, pKa, pH): 
    """Calculate the charge of a residue depending on its pKa
    and pH.
    """
    q_i = 0.0
    # CYS residues forming disulfide bond are neutral.
    if pKa == 99.99:
        return q_i
    #exponent = numpy.power(10, pKa - pH)
    exponent = power(10, pKa - pH)
    q_i = exponent/(1.0 + exponent) 
    # Boolean algebra requires (...) when using 'OR' operator.
    if res_nam in ['ASP', 'GLU', 'C-' , 'TYR', 'Oco', 'CYS']:
        q_i -= 1.0
    return q_i

def write_pqr(target, pH, pqr):
    pqr_file = open(pdb_base_path + target + '-%05.2f-reo.pqr' % pH, 'w')
    pqr_file.write(pqr)
    pqr_file.close() 

# <<EDIT>>
# 29.05.2012: sim.av_RQ -> sim.rho;
#             Adjusted coordinate parsing to official PDB documentation description.
#def get_num_prot(av_RQ, nw_len, nw_rad):
def get_num_prot(rho, nw_len, nw_rad):
    """Compute the bounding box dimension parallel to the NW
    surface, i.e. the area of the x-, y-plane.
    'av_RQ' is formatted as PQR data, but provides only a
    default charge of '1.0'.
    """ 
    X = []
    Y = []
    Z = []
    #for atm in av_RQ:
    #    if len(atm.split()) > 0:
    #        if atm.split()[0] == 'ATOM':
    #            #x = atm[30:38].strip()
    #            #y = atm[38:46].strip()
    #            #z = atm[46:54].strip()
    #            x = atm[31:39].strip()
    #            y = atm[39:47].strip()
    #            z = atm[47:55].strip()
    #            print x, y, z
    #            X.append(float(x)*0.1)
    #            Y.append(float(y)*0.1)
    #            Z.append(float(z)*0.1) 
    for atm in rho:
        x, y, z = atm[0], atm[1], atm[2]
        X.append(float(x)*0.1)
        Y.append(float(y)*0.1)
        Z.append(float(z)*0.1) 
    # Scan all residues.
    x_min = min(X)
    x_max = max(X)
    y_min = min(Y)
    y_max = max(Y)
    z_min = min(Z)
    z_max = max(Z) 
    prot_xy = (x_max - x_min)*(y_max - y_min)
    nw_surface = 2*pi*nw_rad*nw_len
    n_bio_molecules = nw_surface/prot_xy 
    return n_bio_molecules

def fix_pdb(target):
    """Discard non-'ATOM' and non-'^TER' lines, else PDB2PQR does
    not include multiple chains."""
    # <<PATH>>
    if not os.path.exists('./bfs_fix/' + target + '-fix.pdb'):
        awk_cmd = ['awk', '/(ATOM|^TER)/,//']
        #awkp = Popen(awk_cmd, stdin=open('./' + target + '.pdb', 'r'),
        awkp = Popen(awk_cmd, stdin=open('./bfs_rec/' + target + '-rec.pdb', 'r'),
                     stdout=PIPE, stderr=PIPE, shell=False)
        nat = open('./bfs_rec/' + target + '-nat.pdb', 'w')
        nat.write(awkp.stdout.read())
        awkp.stdout.close()
        nat.close()
        # <<PATH>>
        #pqr_cmd = ['python', '/Users/mzh/software/pdb2pqr/pdb2pqr.py',
        #           '-v', '--chain', '--ff=CHARMM',
        #           target + '-nat.pdb', target + '-fix.pdb']
        # Deployment: copy ~/software/pdb2pqr to /opt/
        #pqr_cmd = ['python', '/home/mzhpropka/software/pdb2pqr/pdb2pqr.py',
        pqr_cmd = ['python', '/opt/pdb2pqr/pdb2pqr.py',
                   '-v', '--chain', '--ff=CHARMM',
                   './bfs_rec/' + target + '-nat.pdb', './bfs_fix/' + target + '-fix.pdb']
        #pqrp = Popen(pqr_cmd, stdout=PIPE, stdin=PIPE,
        #             stderr=open('pdb2pqr_err.dat', 'w'), shell=False)
        pqrp = Popen(pqr_cmd, stdout=PIPE, stdin=PIPE, shell=False)
        #bio_tst.test_fixPDB(pqrp)
        pqrp.communicate()
        pqrp.stdout.close()

def rechain(target):
    """Rechaining:
    [MODEL1 (Chain A, Chain B, ...), MODEL2 (Chain A, Chain B, ...), ...]
    into 
    [MODEL1 (Chain A, Chain B, Chain C, ...)]
    Discarding **non**-'ATOM', -'TER', -'ENDMDL' and -'MODEL' lines from the
    PDB file. Piping the awk stream directly to the 'val' variable used
    for rechaining."""
    # <<PATH>>
    if not os.path.exists('./bfs_rec/' + target + '-rec.pdb'):
        awk_cmd = ['awk', '/(ATOM|^TER|ENDMDL|MODEL)/,//']
        awkp = Popen(awk_cmd, stdin=open('./bfs_pdb/' + target + '.pdb', 'r'),
                     stdout=PIPE, stderr=PIPE, shell=False)
        val = awkp.stdout.readlines()
        rec = open('./bfs_rec/' + target + '-rec.pdb', 'w')
        # For security, 'model_id' is string, so it can not be used in 
        # mathematical operations. Here maxium 17 chains.
        codes = zip(range(1, 18), string.uppercase)
        models = [] 
        # Determine how many 'MODEL's are in the structure.
        mdl_cnt = 1
        for line in val:
            if len(line.split()) > 0:
                if line.split()[0] == 'MODEL':
                    models.append([str(mdl_cnt)])
                    mdl_cnt+=1
        # Determine how many 'TER's are in each model.
        mdl_cnt = -1
        ter_cnt = 0
        for line in val:
            if len(line.split()) > 0:
                if line.split()[0] == ('MODEL'):
                    # Set 'TER' and 'MODEL' counter.
                    #ter_cnt = 0
                    mdl_cnt += 1
                if line.split()[0] == ('TER'):
                    ter_cnt += 1
                models[mdl_cnt].append(ter_cnt) 
        # Generate unique chain identifiers.  
        atm_tot = 0
        for mdl in models:
            atm_tmp = 0
            # mdl:                        ['1', 0, 0, ..., 1, 1]
            # mdl_tmp:                    [     0, 0, ..., 1, 1]
            # mdl_tmp[atm_cnt]:           0
            # codes:                      [(1, 'A'), (2, 'B'), (3, 'C'), ...]
            # codes[mdl_tmp[atm_cnt]]:    (1, 'A')
            # codes[mdl_tmp[atm_cnt]][1]: 'A' 
            # NOTICE:
            # - 'atm_tot': Counts over all lines in the file.
            # - 'atm_tmp': Counts over all lines of a MODEL.
            for atm in mdl[1:]:
                #<CHECK_LINE>
                line_length = len(val[atm_tot])
                if line_length < 20:
                    op = val[atm_tot]
                else:
                    op = val[atm_tot][:20]\
                            + ' ' + codes[mdl[1:][atm_tmp]][1]\
                            + ' ' + val[atm_tot][23:]
                #op = val[atm_tot][:20]\
                #        + ' ' + codes[mdl[1:][atm_tmp]][1]\
                #        + ' ' + val[atm_tot][23:]
                rec.write(op)
                atm_tmp+=1
                atm_tot+=1 
        rec.close()

def get_nw_surface(z_dim, target):
    offset = -z_dim/2.0
    #print "<pre>bio_lib.get_nw_surface.offset: %4.2f</pre>" % offset
    nw  = ''
    nw += '4050\n'
    nw += 'NW Surface\n'
    # To display an atom is required only every 4. Angstrom.
    cnt = 0
    #for i in numpy.arange(-180.0, 180.0, 4.0):
    #    for j in numpy.arange(-90.0, 90.0, 4.0):
    for i in scipy.arange(-180.0, 180.0, 4.0):
        for j in scipy.arange(-90.0, 90.0, 4.0):
            #nw += 'C ' + '%2.1f ' % i + '%2.1f ' % j + '%2.1f\n' % offset
            nw += 'C %2.1f %2.1f %2.1f\n' % (i,j,offset)
    nw_dat = open(pdb_base_path + 'nw_%s.xyz' % target, 'w')
    nw_dat.write(nw)
    nw_dat.close() 

def download_pdb(target):
    if not os.path.exists(pdb_base_path + target + '.pdb'):
        address='http://www.pdb.org/pdb/files/%s.pdb1' % target
        dat = urllib.urlopen(address)
        val = open(pdb_base_path + target + '.pdb', 'w')
        for i in dat.readlines():
            val.write(i) 
        dat.close()
        val.close()
# ........................................................................  
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: Status Page, HTML, Jmol and GNUPlot
# - Status page
# - HTML page
# - Jmol script
# - GNUPlot output
# ........................................................................  
def get_status(target, new_message):
    stat_ini  = '<pre>\n'
    stat_ini += new_message + "\n"
    stat_ini += '</pre>\n'
    return stat_ini

def prepare_setup(target, pqr, pH, Q_tot): #, comment): #, mode):
    # Populate comment area.
    timestamp = ("%s" % datetime.datetime.now()).split('.')[0]
    comment  = '# BioFET-SIM Calculation\n'
    comment += '# Date of calculation: %s\n' % timestamp
    comment += '# Calculation target: %s\n' % target
    comment += '# pH: %s\n' % pH
    comment += '# <Add comment here>'
    # Substitute parameters.
    params = get_default_parameters()
    params['target']  = target
    params['pqr']     = pqr
    params['pH']      = pH
    params['Q_tot']   = Q_tot
    params['comment'] = comment
    #page_dat = open('page_setup.html', 'r')
    page_dat = open('resp.html', 'r')
    page_val = page_dat.read()
    page_dat.close()
    #page_template = Template(page_val) 
    #return page_template.safe_substitute(params)
    return page_val

def prepare_Jmol(target):
    sub = dict(target=target)
    jmol_instr_dat = open(jmol_path + 'mod_template.spt', 'r')
    jmol_instr_val = jmol_instr_dat.read()
    jmol_instr_tem = Template(jmol_instr_val) 
    jmol_instr_dat.close() 
    jmolScript = jmol_instr_tem.safe_substitute(sub) 
    js = open(jmol_path + 'mod_%s.spt' % target, 'w')
    js.write(jmolScript)
    js.close()

def prepare_pH_response_plot(target, pH_resp): #, mode):
    plot = ''
    for i in range(14):
        plot += str(i+1) + ' %4.2f'%pH_resp[i] + '\n'
    res_val = open(results_path + target + '-pH-reo.dat', 'w')
    res_val.write(plot)
    res_val.close()
    # KU machine.
    #gnus  = "set terminal svg\n"
    #gnus += "set output \'" + results_path + "%s-pH-reo.svg\'\n" % target
    #gnus += 'set style line 1 lt 1 lw 2 pt 7 ps 1\n'
    #gnus += "unset title\n"
    #gnus += "set nokey\n"
    #gnus += "set grid\n"
    #gnus += "set xlabel 'pH'\n"
    #gnus += "set ylabel 'Sensitivity(pH)'\n" 
    #gnus += "plot \'" + results_path + "%s-pH-reo.dat\' u ($1):($2) w lp ls 1\n" % target
    #gnus += "set output ''\n"
    #gnup = Popen([gnuplot_exe], stdin=PIPE, stdout=open(results_path + '%s-pH-reo.svg' % target, 'w'), stderr=PIPE, shell=False)

    # PROPKA
    gnus  = 'set terminal postscript eps enhanced color "Times-Roman" 22\n'
    gnus += "set output \'" + results_path + "%s-pH-reo.eps\'\n" % target
    gnus += 'set style line 1 lt 1 lw 2 pt 7 ps 1\n'
    gnus += "unset title\n"
    gnus += "set nokey\n"
    gnus += "set grid\n"
    gnus += "set xlabel 'pH'\n"
    gnus += "set ylabel 'Sensitivity(pH)'\n" 
    gnus += "plot \'" + results_path + "%s-pH-reo.dat\' u ($1):($2) w lp ls 1\n" % target
    gnus += "set output ''\n"
    gnup = Popen([gnuplot_exe], stdin=PIPE, stdout=open(results_path + '%s-pH-reo.eps' % target, 'w'), stderr=PIPE, shell=False)
    gnup.communicate(gnus)
    os.system(convert_exe + ' -resample 200x200 -density 200x200 '\
              + results_path + '%s-pH-reo.eps '% target\
              + results_path + '%s-pH-reo.png' % target)

def prepare_results(target, results, x_val, x_lbl, num_prot, dG_G0, G0, bfs_file_name, t): #, mode):
    num_prot_s = "%2.0f"%num_prot
    sub = dict(target=target, num_prot=num_prot_s,
               x_val=x_val, x_lbl=x_lbl,
               dG_G0=dG_G0, G0=G0, bfs_file_name=bfs_file_name) #, mode=mode)
    labels = {'L_d'    : 'Debye Length [nm]',
              'L_tf'   : 'Thomas-Fermi Length [nm]',
              'lay_ox' : 'Oxide Layer Thickness [nm]',
              'nw_len' : 'NW Length [nm]',
              'nw_rad' : 'NW Radius [nm]'} 
    plot = ""
    for xy in results:
        plot += xy
    res_val = open(results_path + target + '-reo.dat', 'w')
    res_val.write(plot)
    res_val.close()
    # KU machine
    #gnus  = "set terminal svg\n"
    #gnus += "set output \'" + results_path + "%s-reo.svg\'\n" % target
    #gnup = Popen([gnuplot_exe],
    #              stdin=PIPE, stdout=open(results_path + '%s-reo.svg' % target, 'w'),
    #              stderr=PIPE, shell=False)

    # PROPKA
    gnus  = 'set terminal postscript eps enhanced color "Times-Roman" 22\n'
    gnus += "set output \'" + results_path + "%s-reo.eps\'\n" % target
    gnus += 'set style line 1 lt 1 lw 6 pt 7 ps 1\n'
    gnus += "unset title\n"
    gnus += "set nokey\n"
    gnus += "set grid\n"
    gnus += "set xlabel \'%s\'\n" % labels[x_lbl]
    gnus += "set ylabel 'Sensitivity'\n" 
    gnus += "plot \'" + results_path + "%s-reo.dat\' u ($1):($2) w p ls 1\n" % target
    gnus += "set output ''\n"
    gnup = Popen([gnuplot_exe],
                   stdin=PIPE, stdout=open(results_path + '%s-reo.eps' % target, 'w'),
                   stderr=PIPE, shell=False)
    gnup.communicate(gnus)
    os.system(convert_exe + ' -resample 200x200 -density 200x200 '\
              + results_path + '%s-reo.eps '% target\
              + results_path + '%s-reo.png' % target)
# ........................................................................
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: BioFET-SIM parameters
# Status: not used.
# ........................................................................  
def get_default_parameters():
    params = {} 
    # ..............................
    # NW Properties
    # ..............................
    # NW Length [nm]
    params['nw_len']       =           2000.00
    params['delta_nw_len'] =            500.00 
    params['nw_len_x_min'] = (params['nw_len'] - params['delta_nw_len'])
    params['nw_len_x_max'] = (params['nw_len'] + params['delta_nw_len']) 
    # NW Radius [nm] 
    params['nw_rad']       =             10.00
    params['delta_nw_rad'] =              2.00
    params['nw_rad_x_min'] = (params['nw_rad'] - params['delta_nw_rad'])
    params['nw_rad_x_max'] = (params['nw_rad'] + params['delta_nw_rad']) 
    # Thomas-Fermi Length, [nm]; L_t(n_0)
    params['L_tf']         =              2.04
    params['delta_L_tf']   =              0.50
    params['L_tf_x_min']   =    (params['L_tf'] - params['delta_L_tf'])
    params['L_tf_x_max']   =    (params['L_tf'] + params['delta_L_tf']) 
    # NW Permittivity [eps_0] 
    params['eps_1']        =             12.00
    # Charge carrier mobility, [m^2 V^-1 s^-1]
    params['mu']           =           1.00E-2
    # Charge carrier density [m^-3]; n_0(L_tf)
    params['n_0']          =           1.11E24
    # {'P';'N'} NW doping type [unit?] 
    params['nw_type']      =               'P'

    # ..............................
    # Other layer properties
    # ..............................
    # Oxide lyer thickness
    params['lay_ox']       =              2.00
    params['delta_lay_ox'] =              1.00
    params['lay_ox_x_min'] = (params['lay_ox'] - params['delta_lay_ox'])
    params['lay_ox_x_max'] = (params['lay_ox'] + params['delta_lay_ox'])
    # Oxide layer permittivity [eps_0]
    params['eps_2']        =              3.90
    # Biofunctionalization layer thickness [nm]
    params['lay_bf']       =              1.00

    # ..............................
    # Solvent Properties
    # ..............................
    # Debye length [nm]
    #params['L_d']      = numpy.arange(0.1, 8.0, 0.2) # Solvent Debye length [nm]
    params['L_d']     =                   2.00
    params['delta_L_d'] =                 1.00
    params['L_d_x_min']    = (params['L_d'] - params['delta_L_d'])
    params['L_d_x_max']    = (params['L_d'] + params['delta_L_d'])
    # Solvent permittivity
    params['eps_3']    =                 78.00

    # Protein Properties, computed internally or user defined
    params['num_prot'] =                  4000 # Total number of proteins on NW 
    return params

def get_parameters(form):
    params = {}
    for k in form.keys():
        params[k] = form[k]
    return params 
# ........................................................................
# ------------------------------------------------------------------------ 

# ************************************************************************
# SECTION: Show Functions
# ........................................................................  
def show_rho(rho):
    print "<pre>"
    for i in rho:
        print i
    print "</pre>"
# ........................................................................
# ------------------------------------------------------------------------ 

# ************************************************************************
# ........................................................................
# TEST EXECUTE
if __name__ == '__main__': 
    #print calc_pKas('kk8add').read()
    print calc_Q_tot('kk8add', 7.4)
    #nw_rad = 10.0
    #lay_ox = 2.0
    #L_d = 1000.0
    #L_tf = 2.0
    #eps_1 = 12.0
    #eps_2 = 3.9
    #eps_3 = 78.0
    #print Gamma(nw_rad, lay_ox, L_d, L_tf, eps_1, eps_2, eps_3) 
# ------------------------------------------------------------------------ 
