import numpy as np

#----------------------------Input Parameters--------------------------------
casename   = "multistage"
N           = 22363                               #Speed of Impeller [rpm]
P01         = 101325                            #Inlet Pressure [Pa]
T01         = 300                               #Inlet Temperature [K]
Vt1         = 0                                 #Inlet Tangentail Velocity[m^2/s]
mdot        = 4                                #Mass flow rate [kg/s]
delTT       = 176.46                           #Ovearll Temperature rise [K]
Rgas        = 287                              #Gas constant of Air [J/kg-K]
Cp          = 1006                               #Specific constant at constant pressure [J/kg-K]
PR          = 4.8                                #Overall Pressure Ratio
Z           = [24,24,24]                          #Number of Blades [-]
nsect       = 5                               #Number of streamlines
Beta1_Blade = np.array([-36,-62])       #Inlet relative flow angle[deg]

#-----------------------------Stage Parameters---------------------------------
Eta_R1          = 0.97                           #Rotor 1 Efficiency
WorkRatio_R1    = 0.35                     #Ratio of work done by Rotor 1
dalpha          = 25                             #Stator turning angle [deg]
Y               = 0.03                                #Total Pressure loss coefficient across stator [-]
Beta6_Blade     = -30                       #Backsweep angle [deg]

#-----------------------------Flowpath Parameters------------------------------
R_hub_le = 0.0449341
R_tip_le = 0.11274362
R_hub_te = 0.2
R_tip_te = 0.2

[X01, R01] = [-0.065005893244729426, 0.21920467299175289]     #(x,r)-coordinates for center of circle for hub
R = 0.186                         #Radius of the circle

[X04,R04] = [0,0.209]            #(x,r)-coordinates for center of ellipse for shroud
[ae,be] = [0.105761,R04-R_tip_le]

R_DiffuserExit = 0.3

#------------------------------Tblade3 Input-----------------------------------
airfoiltype = 'sect1'
chrdr_nd = 1.165
gap_rs = np.array([0.0025,0.00125])
s1_len = np.array([0.04,0.02])
gap_sr = np.array([0.005,0.003])

#===============================================================================
#------------------------Variable Declaration I---------------------------------
#===============================================================================
nrows = len(Z)                          #Number of blade rows
case = 'mr'
nstns = nrows+1                         #Number of stations
nstations = nrows*2

xm   = np.zeros(nstations)
rm   = np.zeros(nstations)
x_s = np.zeros((nstations,2))
r_s = np.zeros((nstations,2))


#Scalar properties dor 1D
U   = np.zeros(nstations)
V   = np.zeros(nstations)
Vm  = np.zeros(nstations)
Vt  = np.zeros(nstations)
W   = np.zeros(nstations)
Wt  = np.zeros(nstations)
Wm  = np.zeros(nstations)

#Intensive properties
T       = np.zeros(nstations)
T0      = np.zeros(nstations)
P       = np.zeros(nstations)
P0      = np.zeros(nstations)
alpha   = np.zeros(nstations)
beta    = np.zeros(nstations)
M       = np.zeros(nstations)
Mrel    = np.zeros(nstations)
T0rel   = np.zeros(nstations)
P0rel   = np.zeros(nstations)
rho     = np.zeros(nstations)
sw      = np.zeros(nstations)
a       = np.zeros(nstations)
area    = np.zeros(nstations)

#===============================================================================
#------------------------Variable Declaration II--------------------------------
#===============================================================================
#Along meanline
sol = np.zeros(nrows)
pitch = np.zeros(nrows)
Rx = np.zeros(nrows)
phi = np.zeros(nrows)
DH = np.zeros(nrows)                #DeHaller Number
Df = np.zeros(nrows)              #Diffusion Factor

#Along span and diffrent Sections (_s represents along span)
chord = np.zeros((nsect, nrows))
Vt_s = np.zeros((nsect, nstations))
Vm_s = np.zeros((nsect, nstations))
T_s = np.zeros((nsect, nstations))
sw_s = np.zeros((nsect, nstations))
ywall_s = np.zeros((nsect,nstations))
sol_s = np.zeros((nsect,nrows))
pitch_s = np.zeros((nsect,nrows))
#Df = np.zeros((nsect, nrows))
