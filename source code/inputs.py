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
[a,b] = [0.105761,R04-R_tip_le]

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

#Scalar properties dor 1D
U   = np.zeros(nrows)
V   = np.zeros(nrows)
Vm  = np.zeros(nrows)
Vt  = np.zeros(nrows)
W   = np.zeros(nrows)
Wt  = np.zeros(nrows)

#Intensive properties
T       = np.zeros(nrows)
T0      = np.zeros(nrows)
P       = np.zeros(nrows)
P0      = np.zeros(nrows)
alpha   = np.zeros(nrows)
beta    = np.zeros(nrows)
M       = np.zeros(nrows)
Mrel    = np.zeros(nrows)
T0rel   = np.zeros(nrows)
P0rel   = np.zeros(nrows)
sw      = np.zeros(nrows)

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
area = np.zeros(nstations)



#Along span and diffrent Sections (_s represents along span)
Vt_s = np.zeros((nsect, nstns))
Vm_s = np.zeros((nsect, nstns))
T_s = np.zeros((nsect, nstns))
sw_s = np.zeros((nsect, nstns))
ywall_s = np.zeros((nsect,nstns))
sol_s = np.zeros((nsect,nrows))
pitch_s = np.zeros((nsect,nrows))
#Df = np.zeros((nsect, nrows))
