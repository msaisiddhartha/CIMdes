import numpy as np
import sys

f_in = "input.in"
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1

nlines = file_len(f_in)
data = []
dict = {}

f = open(f_in,'r')
for i in range(nlines):
    datastr = f.readline()
    data = datastr.split()
    dict[data[0].strip(':')] = [float(data[1])]
    if len(data)>2:
        for j in range(2,len(data)):
            dict[data[0].strip(':')].append(float(data[j]))            

bnumber = int(np.asarray(dict['Blades']))
Z = np.zeros(bnumber)
cl = np.zeros(bnumber)
ang = np.zeros(bnumber)

WorkRatio = np.zeros((bnumber//2+1,1))
dalpha = np.zeros((bnumber//2,1))
#----------------------------Input Parameters--------------------------------
casename    = "dlrcc"
N           = dict['N'][0]             # Speed of Impeller [rpm]
P01         = dict['PT_in'][0]         # Inlet Pressure [Pa]
T01         = dict['TT_in'][0]        # Inlet Temperature [K]
Vt1         = 0                     # Inlet Tangentail Velocity[m^2/s]
mdot        = dict['mdot'][0]         # Mass flow rate [kg/s]
delTT       = 176.46                # Ovearll Temperature rise [K]
Rgas        = 287                   # Gas constant of Air [J/kg-K]
Cp          = 1006                  # Specific constant at constant pressure [J/kg-K]
for i in range(bnumber):
    Z[i]    = dict['Z'][i]        # Number of Blades starting with rotor [-]
    cl[i]   = dict['cl'][i]         # Average Tip clearance per blade row
    ang[i]  = dict['angles'][i]     # Row type (0-axial; 1-radial; 2-mixed)
nsect       = 5                     # Number of streamlines
nu          = 1.8e-6                # Kinematic viscosity of fluid

#-----------------------------Stage Parameters---------------------------------
if bnumber//2-1>1:
    for i in range(bnumber//2-1):
        WorkRatio[i]   = np.asarray(dict['R1_WR'][i])     # Ratio of work done per stage (Total = 1)
        dalpha[i]      = np.asarray(dict['dalpha'][i])   # Stator turning angle [deg]
else:
    WorkRatio[0]   = np.asarray(dict['R1_WR'][0])     # Ratio of work done per stage (Total = 1)
    dalpha[0]      = np.asarray(dict['dalpha'][0])   # Stator turning angle [deg]
        
WorkRatio[-1]   = 1-np.sum(dict['R1_WR'])
Y           = np.asarray(dict['Y'])              # Total Pressure loss coefficient across stator [-]
Bckswp      = -45               # Backsweep angle [deg]

#-----------------------------Flowpath Parameters------------------------------
R_hub_le = 0.0449341            # Radius of compressor inlet at hub
R_tip_le = 0.11274362           # Radius of compressor inlet at casing
R_hub_te = 0.2                  # Radius of compressor outlet at hub
R_tip_te = 0.2                  # Radius of compressor outlet at casing
R_DiffuserExit  = 0.3           # Radius of diffuser outlet

# (x,r)-coordinates for center of circle for hub
[X01, R01] = [-0.065005893244729426, 0.21920467299175289]
R = 0.186  # Radius of the circle

[X04,R04]       = [0,0.209]            #(x,r)-coordinates for center of ellipse for shroud
[ae,be]         = [0.105761,R04-R_tip_le]

#------------------------------Tblade3 Input-----------------------------------
airfoiltype = 'sect1'
chrdr_nd = 1.165
gap = np.zeros((bnumber+1,2))
gapc = 0
for i in range(bnumber+1):
    for j in range(2):
        gap[i,j] = dict['gap'][gapc]
        gapc+=1
#===============================================================================
#------------------------Variable Declaration I---------------------------------
#===============================================================================
nrows = len(Z)                          #Number of blade rows
case = 'mr'
nstns = nrows+1                         #Number of stations
nstations = nrows*2

xm   = np.zeros((nstations,1))
rm   = np.zeros((nstations,1))
x_s = np.zeros((nstations,2))
r_s = np.zeros((nstations,2))


#Scalar properties dor 1D
U   = np.zeros((nstations,1))
V   = np.zeros((nstations,1))
Vm  = np.zeros((nstations,1))
Vt  = np.zeros((nstations,1))
Vz  = np.zeros((nstations,1))
Vr  = np.zeros((nstations,1))
W   = np.zeros((nstations,1))
Wt  = np.zeros((nstations,1))
Wm  = np.zeros((nstations,1))

#Intensive properties
T       = np.zeros((nstations,1))
T0      = np.zeros((nstations,1))
P       = np.zeros((nstations,1))
P0      = np.zeros((nstations,1))
alpham   = np.zeros((nstations,1))
betam    = np.zeros((nstations,1))
alphaz   = np.zeros((nstations,1))
betaz    = np.zeros((nstations,1))
M       = np.zeros((nstations,1))
Mrel    = np.zeros((nstations,1))
T0rel   = np.zeros((nstations,1))
P0rel   = np.zeros((nstations,1))
rho     = np.zeros((nstations,1))
sw      = np.zeros((nstations,1))
a       = np.zeros((nstations,1))
area    = np.zeros((nstations,1))


#===============================================================================
#------------------------Variable Declaration II--------------------------------
#===============================================================================
#Along meanline
sol = np.zeros(nrows)               #Solidity
pitch = np.zeros(nrows)             #Pitch = s/c
Rx = np.zeros(nrows)                #Degree of reaction
phi = np.zeros(nrows)               #Flow co-efficient
DH = np.zeros(nrows)                #DeHaller Number
Df = np.zeros(nrows)                #Diffusion Factor
Cf = np.zeros(nrows)                #Skin friction coefficient
Ib = np.zeros(nrows//2+1)           #Work input factot = Cp*delTT

# Enthalpy loss for 7 indiviudal losses row-wise and overall
dH_Loss = np.zeros((7, nrows + 1))
# First elemt is for classificatio aero, int and ext losses and second element
dH = np.zeros((3, nrows + 1))
TR = np.zeros(nrows + 1)

#Along span and diffrent Sections (_s represents along span)
chord = np.zeros((nsect, nrows))
Vt_s = np.zeros((nsect, nstations))
Wt_s = np.zeros((nsect, nstations))
Vm_s = np.zeros((nsect, nstations))
T_s = np.zeros((nsect, nstations))
sw_s = np.zeros((nsect, nstations))
Wt_s = np.zeros((nsect, nstations))
betam_s = np.zeros((nsect, nstations))
alpham_s = np.zeros((nsect, nstations))
V_s = np.zeros((nsect, nstations))
Vr_s = np.zeros((nsect, nstations))
Vz_s = np.zeros((nsect, nstations))
W_s = np.zeros((nsect, nstations))
a_s = np.zeros((nsect, nstations))
M_s = np.zeros((nsect, nstations))
Mrel_s = np.zeros((nsect, nstations))
ywall_s = np.zeros((nsect, nstations))
sol_s = np.zeros((nsect, nrows))
pitch_s = np.zeros((nsect, nrows))

Vm_s1 = np.zeros((nsect, nstations))
T_s1 = np.zeros((nsect, nstations))
Vt_s1 = np.zeros((nsect, nstations))
sw_s1 = np.zeros((nsect, nstations))
Wt_s1 = np.zeros((nsect, nstations))
betam_s1 = np.zeros((nsect, nstations))
alpham_s1 = np.zeros((nsect, nstations))
V_s1 = np.zeros((nsect, nstations))
W_s1 = np.zeros((nsect, nstations))
a_s1 = np.zeros((nsect, nstations))
M_s1 = np.zeros((nsect, nstations))
Mrel_s1 = np.zeros((nsect, nstations))
#Df = np.zeros((nsect, nrows))

if nsect % 2 == 0:
    mean_num = nsect // 2
else:
    mean_num = nsect // 2 + 1

fout = []
data = []
