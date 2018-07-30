import numpy as np
from inputs import *

npts = 50                               #Number of pts on streamline
np_inlt = 5
np_inlt_eq = 2
stp = 0.1
np_outlt = 25
#------------------------------------------------------------------------------
#==============================================================================
#Generating streamlines at hub nad tip using the krain single-stage compressor
#==============================================================================
#Hub
r_hub = np.linspace(R_hub_le,R_hub_te,npts)
x_hub =  X01 + (R**2-(r_hub-R01)**2)**0.5

#Casing
r_tip = np.linspace(R_tip_le,R_tip_te,npts)
x_tip = X04 + (a/b)*(b**2-(r_tip-R04)**2)**0.5

#Non-dimesionalizing with leading edge hub radius as scaling factor
bsf = r_tip[0]
x_hub_nd = x_hub/bsf
r_hub_nd = r_hub/bsf
x_tip_nd = x_tip/bsf
r_tip_nd = r_tip/bsf



#------------------------------------------------------------------------------
#==============================================================================
#Splitting the blade into 3 parts for rotor and stator
#==============================================================================
#----------------------------Rotor 1-------------------------------------------
x = np.array((npts,nstations))
r = np.array((npts,nstations))

r1 = np.array([r_hub[0],r_tip[0]])
x1 = np.array([x_hub[0],x_tip[0]])

r2 = np.array([0.075,0.12])
x2 = np.array([X01 + (R**2-(r2[0]-R01)**2)**0.5,X04 + (a/b)*(b**2-
                  (r2[1]-R04)**2)**0.5])

#-----------------------------Stator-------------------------------------------
r3 = r2+gap_rs
x3 = np.array([X01 + (R**2-(r3[0]-R01)**2)**0.5,X04 + (a/b)*(b**2-
                  (r3[1]-R04)**2)**0.5])

r4 = r3+s1_len
x4 = np.array([X01 + (R**2-(r4[0]-R01)**2)**0.5,X04 + (a/b)*(b**2-
                  (r4[1]-R04)**2)**0.5])

#--------------------------------Rotor 2---------------------------------------
r5 = r4+gap_sr
x5 = np.array([X01 + (R**2-(r5[0]-R01)**2)**0.5,X04 + (a/b)*(b**2-
                  (r5[1]-R04)**2)**0.5])
r6 = np.array([r_hub[-1],r_tip[-1]])
x6 = np.array([x_hub[-1],x_tip[-1]])

r_nd = np.array([r1,r2,r3,r4,r5,r6])/bsf
x_nd = np.array([x1,x2,x3,x4,x5,x6])/bsf
#------------------------------------------------------------------------------
#----------------Area calculation at inlet and oultlet-------------------------
# Station 1(Inlet)
r1m = np.mean(r1)
A1 = np.pi * (r1[1]**2 - r1[0]**2)

# Staion 2(R1 exit)
r2m = np.mean(r2)
A2 = np.pi * (sum(r2)) * ((r2[0] - r2[1])**2 + (x2[0] - x2[1])**2)**0.5

# Staion 3(S1 inlet)
r3m = np.mean(r3)
A3 = np.pi * (sum(r3)) * ((r3[0] - r3[1])**2 + (x3[0] - x3[1])**2)**0.5

# Staion 4(S1 outlet)
r4m = np.mean(r4)
A4 = np.pi * (sum(r4)) * ((r4[0] - r4[1])**2 + (x4[0] - x4[1])**2)**0.5

# Staion 5(R2 inlet)
r5m = np.mean(r5)
A5 = np.pi * (sum(r5)) * ((r5[0] - r5[1])**2 + (x5[0] - x5[1])**2)**0.5

# Station 6(Outlet)
r6m = np.mean(r6)
b6 = x6[0] - x6[1]
A6 = 2 * np.pi * b6 * r6[0]

#-----------------------------Streamlines--------------------------------------

# Inlet Streamlines
x_inlt = np.zeros(np_inlt) + (-stp * np_inlt)
r_hub_inlt = np.zeros(np_inlt)
r_tip_inlt = np.ones(np_inlt) * r_tip_nd[0]
for i in range(0, np_inlt - 1):
    x_inlt[i + 1] = x_inlt[i] + stp

for j in range(np_inlt - 1, -1, -1):
    if j >= np_inlt - np_inlt_eq:
        r_hub_inlt[j] = R01 / bsf - \
            ((R / bsf)**2 - (x_inlt[j] - (X01 / bsf))**2)**0.5
    else:
        r_hub_inlt[j] = r_hub_inlt[j + 1]

# Outlet Streamlines
R_DiffuserExit_nd = R_DiffuserExit / bsf
r_outlt = np.linspace(r_hub_nd[-1], R_DiffuserExit_nd, np_outlt + 1)
r_outlt = r_outlt[1:]
x_hub_outlt = np.zeros(np_outlt) + x_hub_nd[-1]
x_tip_outlt = np.zeros(np_outlt)

# Diffuser
Aoutlt = 2 * np.pi * (x_hub_nd[-1] - x_tip_nd[-1]) * r_hub_nd[-1]
Aexit = 0.99 * (Aoutlt)

area = np.linspace(Aoutlt, Aexit, np_outlt)
x_tip_outlt = x_hub_outlt - area / (2 * np.pi * r_outlt)

x_hub_nd = np.hstack((x_inlt, x_hub_nd, x_hub_outlt))
r_hub_nd = np.hstack((r_hub_inlt, r_hub_nd, r_outlt))
x_tip_nd = np.hstack((x_inlt, x_tip_nd, x_tip_outlt))
r_tip_nd = np.hstack((r_tip_inlt, r_tip_nd, r_outlt))

xsl = np.zeros((nsect, len(x_hub_nd)))
rsl = np.zeros((nsect, len(x_hub_nd)))
for j in range(len(x_hub_nd)):
    xsl[:, j] = np.linspace(x_hub_nd[j], x_tip_nd[j], nsect)
    rsl[:, j] = np.linspace(r_hub_nd[j], r_tip_nd[j], nsect)
