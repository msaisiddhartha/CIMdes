"""Multi-rotor compressor. All units are in SI system."""
import matplotlib.pyplot as py
import numpy as np
import subprocess

#----------------------------Input Parameters--------------------------------
casename = "multistage"
N = 22363                               #Speed of Impeller [rpm]
mdot = 4                                #Mass flow rate [kg/s]
Rgas = 287                              #Gas constant of Air [J/kg-K]
Cp = 1006                               #Specific constant at constant pressure [J/kg-K]
PR = 4.8                                #Overall Pressure Ratio
Z = [24,24,24]                          #Number of Blades [-]
nsect = 5                               #Number of streamlines

#-----------------------------Stage Parameters---------------------------------
PR_SingleStage = 4.7                    #Overall Pressure ratio Single stage
Eta_SingleStage = 0.94                  #Overall efficiency single stage

Eta_R1 = 0.97                           #Rotor 1 Efficiency
Eta_R2 = 0.95                           #Rotor 2 Efficiency
WorkRatio_R1 = 0.35                     #Ratio of work done by Rotor 1
dalpha = 25                             #Stator turning angle [deg]
Y = 0.03                                #Total Pressure loss coefficient across stator [-]
Beta6_Blade = -30                       #Backsweep angle [deg]

#-----------------------------Inlet Conditions---------------------------------
P01 = 101325                            #Inlet Pressure [Pa]
T01 = 300                               #Inlet Temperature [K]
Vt1 = 0                                 #Inlet Tangentail Velocity[m^2/s]
Beta1_Blade = np.array([-36,-62] )      #Inlet relative flow angle[deg]

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
nrows = len(Z)                          #Number of blade rows
nstns = nrows+1                         #Number of stations
nstations = nrows*2
airfoiltype = 'sect1'
case = "mr"
chrdr_nd = 1.165

#--------------------------------Streamlines-----------------------------------
npts = 50                               #Number of pts on streamline
np_inlt = 5
np_inlt_eq = 2
stp = 0.1
np_outlt = 25
gap_rs = np.array([0.0025,0.00125])
s1_len = np.array([0.04,0.02])
gap_sr = np.array([0.005,0.003])
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

#==============================================================================
#Calculating Thermodynamic Properties at inlet and outlet
#==============================================================================
"""The flow assumed here is isentropic flow. Inlet relative angles at hub and
shroud are given 36deg and 62deg respectively same as input from krain paper.
Pressure ratio is 4.7 and efficiency is calculated based upon inputs. The air
considered is a perfect gas."""

#Evaluating based on known parameters
w = 2*np.pi*N/60                        #Angular velocity[rad/s]
g = 1/(1-Rgas/Cp)                       #Ratio of Specific heats[-]


#0-D Calculations Thremodynamic properties
TR_SingleStage = (PR_SingleStage**((g-1)/g)-1)/Eta_SingleStage+1
Work_SingleStage = Cp*T01*(TR_SingleStage-1)
Work_MultiStage = Work_SingleStage

TR = (Work_MultiStage/(Cp*T01))+1
Eta = (PR**((g-1)/g)-1)/(TR-1)
P06 = PR*P01                            #Outlet Total pressure[Pa]
T06 = TR*T01                            #Outlet Total Temperature[K]
Work_MultiStage = Cp*(T06-T01)
Work_R1 = WorkRatio_R1*Work_MultiStage
Work_R2 = Work_MultiStage-Work_R1
T05 = T06-Work_R2/Cp
T02 = T05
TR_R1 = T02/T01
TR_R2 = T06/T05

#----------------Area calculation at inlet and oultlet-------------------------
#Station 1(Inlet)
r1m = np.mean(r1)
A1 = np.pi*(r1[1]**2-r1[0]**2)

#Staion 2(R1 exit)
r2m = np.mean(r2)
A2 = np.pi*(sum(r2))*((r2[0]-r2[1])**2+(x2[0]-x2[1])**2)**0.5

#Staion 3(S1 inlet)
r3m = np.mean(r3)
A3 = np.pi*(sum(r3))*((r3[0]-r3[1])**2+(x3[0]-x3[1])**2)**0.5

#Staion 4(S1 outlet)
r4m = np.mean(r4)
A4 = np.pi*(sum(r4))*((r4[0]-r4[1])**2+(x4[0]-x4[1])**2)**0.5

#Staion 5(R2 inlet)
r5m = np.mean(r5)
A5 = np.pi*(sum(r5))*((r5[0]-r5[1])**2+(x5[0]-x5[1])**2)**0.5

#Station 6(Outlet)
r6m = np.mean(r6)
b6 = x6[0]-x6[1]
A6 = 2*np.pi*b6*r6[0]

#------------------------------------------------------------------------------
#==============================================================================
#Meanline and Velocty Triangles Calculations
#==============================================================================
#-----------------------------------Inlet--------------------------------------
U1 = w*r1m
sw1 = r1*Vt1
Wt1 = Vt1-U1
beta1 = np.average(Beta1_Blade)
Vm1 = Wt1/np.tan(np.radians(beta1))
rho1 = mdot/(Vm1*A1)
V1 = (Vt1**2+Vm1**2)**0.5
T1 = T01-V1**2/(2*Cp)
P1 = P01*((T1/T01)**((g/(g-1))))
W1 = Wt1/np.sin(np.radians(beta1))
a1 = (g*Rgas*T1)**0.5
M1 = V1/a1
alpha1 = np.degrees(np.arctan(Vt1/Vm1))
M1rel = W1/a1
T01rel = T1+(W1**2/(2*Cp))
P01rel = P1+(0.5*rho1*W1**2)

#------------------------------------Station 2---------------------------------
U2 = w*r2m
Vt2 = (Work_R1+U1*Vt1)/U2
PR_R1 = (Eta_R1*(TR_R1-1)+1)**(g/(g-1))
Wt2 = Vt2-U2
P02 = PR_R1*P01
rho02 = P02/(Rgas*T02)
rho2 = rho02
error = 5
while abs(error)>1e-6:
    Vm2 = mdot/(rho2*A2)
    V2 = (Vm2**2+Vt2**2)**0.5
    T2 = T02-V2**2/(2*Cp)
    P2 = P02*((T2/T02)**((g/(g-1))))
    rho2p = P2/(Rgas*T2)
    error = (1-rho2p/rho2)*100
    rho2 = rho2p
a2 = (g*Rgas*T2)**0.5
Wm2 = Vm2
W2 = (Wt2**2+Wm2**2)**0.5
M2 = V2/a2
M2rel = W2/a2
beta2 = np.degrees(np.arctan(Wt2/Vm2))
alpha2 = np.degrees(np.arctan(Vt2/Vm2))
T02rel = T2+(W2**2/(2*Cp))
P02rel = P2+(0.5*rho2*W2**2)
sw2 = r2*Vt2

#-----------------------------R2 Inlet(Station 2)------------------------------
U5 = w*r5m
T05 = T02
P05 = P02-Y*(P02-P2)

alpha5 = alpha2-dalpha
rho05 = P05/(Rgas*T05)
rho5 = rho05
error = 5
while abs(error)>1e-6:
    Vm5 = mdot/(rho5*A5)
    Vt5 = np.tan(np.deg2rad(alpha5))*Vm5
    V5 = (Vm5**2+Vt5**2)**0.5
    T5 = T05-V5**2/(2*Cp)
    P5 = P05*((T5/T05)**((g/(g-1))))
    rho5p = P5/(Rgas*T5)
    error = (1-rho5p/rho5)*100
    rho5 = rho5p
a5 = (g*Rgas*T5)**0.5
Wt5 = Vt5-U5
PR_R2 = P06/P05
TR_R2 = T06/T05
Eta_R2 = (PR_R2**((g-1)/g)-1)/(TR_R2-1)
Wm5 = Vm5
W5 = (Wt5**2+Wm5**2)**0.5
M5 = V5/a5
M5rel = W5/a5
beta5 = np.degrees(np.arctan(Wt5/Vm5))
alpha5 = np.degrees(np.arctan(Vt5/Vm5))
T05rel = T5+(W5**2/(2*Cp))
P05rel = P5+(0.5*rho5*W5**2)
sw5 = r5*Vt5

#-----------------------------R2 Outlet(Station 4)-----------------------------
U6 = w*r6m
Vt6 = (Cp*(T06-T05)+U5*Vt5)/U6
Wt6 = Vt6-U6
beta6 =  Beta6_Blade
Vm6 = Wt6/np.tan(np.deg2rad(beta6))
rho06 = P06/(Rgas*T06)
rho6 = mdot/(A6*Vm6)
alpha6 = np.rad2deg(np.arctan(Vt6/Vm6))
Wm6 = Vm6
V6 = (Vt6**2+Vm6**2)**0.5
W6 = (Wt6**2+Wm6**2)**0.5
T6 = T06-V6**2/(2*Cp)
P6 = P06*((T6/T06)**((g/(g-1))))
a6 = (g*Rgas*T6)**0.5
M6 = V6/a6
M6rel = W6/a6
T06rel = T6+(W6**2/(2*Cp))
P06rel = P6+(0.5*rho6*W6**2)
sw6 = r6*Vt6

#-------------------------------Flow properties--------------------------------
#At hub and shroud
Vt = np.zeros((nsect,nstns))
Vm = np.zeros((nsect,nstns))
T = np.zeros((nsect,nstns))
sw = np.zeros((nsect,nstns))

span = np.linspace(0,1,nsect)

r_1 = span*(r1[1]-r1[0])+r1[0]
r_2 = span*(r2[1]-r2[0])+r2[0]
r_5 = span*(r5[1]-r5[0])+r5[0]
r_6 = span*(r6[1]-r6[0])+r6[0]
r_blde = np.transpose(np.array((r_1,r_2,r_5,r_6)))

Ur = w*r_blde

Vm[:,0]=Vm1
Vm[:,1]=Vm2
Vm[:,2]=Vm5
Vm[:,3]=Vm6

Vt[:,0]=r1m*Vt1/r_1
Vt[:,1]=r2m*Vt2/r_2
Vt[:,2]=r5m*Vt5/r_5
Vt[:,3]=r6m*Vt6/r_6

sw = r_blde*Vt
Wt = Vt-Ur
beta = np.degrees(np.arctan(Wt/Vm))
alpha = np.degrees(np.arctan(Vt/Vm))
V = (Vt**2+Vm**2)**0.5
W = (Wt**2+Vm**2)**0.5



T[:,0] = T1
T[:,1] = T2
T[:,2] = T5
T[:,3] = T6

a = (g*T*Rgas)**0.5
M = V/a
Mrel = W/a

#------------------------------------------------------------------------------
#==============================================================================
# Creating T-Blade3 input file....
#==============================================================================

#-----------------------------Streamlines--------------------------------------
#Inlet Streamlines
x_inlt = np.zeros(np_inlt)+(-stp*np_inlt)
r_hub_inlt = np.zeros(np_inlt)
r_tip_inlt = np.ones(np_inlt)*r_tip_nd[0]
for i in range(0,np_inlt-1):
    x_inlt[i+1] = x_inlt[i]+stp

for j in range(np_inlt-1,-1,-1):
    if j>=np_inlt-np_inlt_eq:
        r_hub_inlt[j] = R01/bsf-((R/bsf)**2-(x_inlt[j]-(X01/bsf))**2)**0.5
    else:
        r_hub_inlt[j] = r_hub_inlt[j+1]

#Outlet Streamlines
R_DiffuserExit_nd = R_DiffuserExit/bsf
r_outlt = np.linspace(r_hub_nd[-1],R_DiffuserExit_nd,np_outlt+1)
r_outlt = r_outlt[1:]
x_hub_outlt = np.zeros(np_outlt)+x_hub_nd[-1]
x_tip_outlt = np.zeros(np_outlt)

#Diffuser
Aoutlt = 2*np.pi*(x_hub_nd[-1]-x_tip_nd[-1])*r_hub_nd[-1]
Aexit = 0.99*(Aoutlt)

area = np.linspace(Aoutlt,Aexit,np_outlt)
x_tip_outlt = x_hub_outlt-area/(2*np.pi*r_outlt)

x_hub_nd=np.hstack((x_inlt,x_hub_nd,x_hub_outlt))
r_hub_nd=np.hstack((r_hub_inlt,r_hub_nd,r_outlt))
x_tip_nd=np.hstack((x_inlt,x_tip_nd,x_tip_outlt))
r_tip_nd=np.hstack((r_tip_inlt,r_tip_nd,r_outlt))

xsl = np.zeros((nsect,len(x_hub_nd)))
rsl = np.zeros((nsect,len(x_hub_nd)))
for j in range(len(x_hub_nd)):
    xsl[:,j] = np.linspace(x_hub_nd[j],x_tip_nd[j],nsect)
    rsl[:,j] = np.linspace(r_hub_nd[j],r_tip_nd[j],nsect)

#-------------------------------Inputs to T-blade3-----------------------------
cntr=0
stage_num=1
fout = []
data = []
chord = np.zeros((nsect,nrows))
ctrlspan = np.linspace(0,1,nsect)
def rowname(row_num, stage_num):
    if row_num%2==0:
        row_name = "R"+str(stage_num)
    if row_num%2==1:
        row_name = "S"+str(stage_num)
        stage_num+=1
    return row_name, stage_num

#Flow angle switch
def angles(curr_row):
    if curr_row==1:
        return 0
    if curr_row==2:
        ang = 2
    else:
        ang = 1
    return ang

print()
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print()
print("Creating geomturbo file...")
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

for k in range(nrows):
    row_name = rowname(k, stage_num)
    stage_num = row_name[1]
    current_sec=0
    f = open("tblade3input." + str(k+1) + "." + str(casename) + ".dat", "w")
    f.write("Input parameters (version 1.1)" + '\n')
    f.write("    " + str(row_name[0]) + '\n')
    f.write(" Blade row #:" + '\n')
    f.write("    " + str(k+1) + '\n')
    f.write(" Number of blades in this row:" + '\n')
    f.write("    " + str(Z[k]) + '\n')
    f.write(" Blade Scaling factor (mm):" + '\n')
    f.write("    " + str(bsf*1000) + '\n')
    f.write(" Number of streamlines:" + '\n')
    f.write("    " + str(nsect) + '\n')
    f.write(" Angles in the input file (0=Beta_z (default),1=Beta_r):" + '\n')
    f.write("    " + str(angles(k+1)) + '\n')
    f.write(" Airfoil camber defined by curvature control (0=no,1=yes):" + '\n')
    f.write("    " + str(0) + '\n')
    f.write(" Airfoil Thickness distribution (0=Wennerstrom,1=Spline):" + '\n')
    f.write("    " + str(0) + '\n')
    f.write(" Airfoil Thickness multiplier (0=no,1=yes):" + '\n')
    f.write("    " + str(0) + '\n')
    f.write(" Airfoil LE defined by spline (0=no,1=yes):" + '\n')
    f.write("    " + str(0) + '\n')
    f.write(" Non-dimensional Actual chord (0=no,1=yes,2=spline):" + '\n')
    f.write("    " + str(0) + '\n')
    f.write(" Sectionwise properties:" + '\n')
    f.write(" J      in_Beta     out_Beta     mrel_in      chord      t/c_max     Incidence     Deviation    Sec. Flow Angle" + '\n')
    for i in range(nsect):
        # stagger is positive in T-Blade3 for Compressors (convention)
        #f.write(str(i + stagger[i] + stagger[i] + mrel[i] + chrdr_nd[i] + "0.2000" + "0.0000" + "0.0000" + "0.0000") + '\n')
        f.write(('%02d'%(i+1)) + "   ")
        if k%2==0:
            f.write(('%2.8f'%beta[i][k]) + "  ")
            f.write(('%2.8f'%beta[i][k+1]) + "  ")
            f.write(('%2.8f'%Mrel[i][k]) + "  ")
        if k%2==1:
            f.write(('%2.8f'%alpha[i][k]) + "  ")
            f.write(('%2.8f'%alpha[i][k+1]) + "  ")
            f.write(('%2.8f'%M[i][k]) + "  ")
        f.write(('%2.8f'%chrdr_nd) + "  ")
        f.write(('%2.8f'%0.0150) + "  ")
        f.write(('%2.8f'%0.0000) + "  ")
        f.write(('%2.8f'%0.0000) + "  ")
        f.write(('%2.8f'%0.0000) + "\n")
    f.write('\n')
    f.write(" LE / TE curve (x,r) definition :" + '\n')
    f.write(" Number of Curve points :" + '\n')
    f.write("    " + str('2') + '\n')
    f.write("   xLE          rLE           xTE          rTE" + '\n')
    for i in range(2):
        f.write("    " + ('%2.8f'%x_nd[cntr][i]) + "  ")
        f.write(('%2.8f'%r_nd[cntr][i]) + "  ")
        f.write(('%2.8f'%x_nd[cntr+1][i]) + "  ")
        f.write(('%2.8f'%r_nd[cntr+1][i]) + "\n")
    cntr+=2
    f.write('\n')
    f.write(" # Airfoil type and Variable Radial Stacking information.         #" + '\n')
    f.write(" # stack_u: % chord stack (0.00 to 100.00).                       #" + '\n')
    f.write(" # stack_v: % below or above meanline stack (-100.00 to +100.00). #" + '\n')
    f.write(" # Use +200 for stacking on airfoil area centroid.                #" + '\n')
    f.write(" Variable Radial stacking (0=no,1=yes):" + '\n')
    f.write("    " + ('%01d'%0) + '\n')
    f.write(" J   type |stk_u |stk_v |umxthk |lethk |tethk  |Jcells(Grid:4n+1) |eta_ofst(<=10){%thkc/Jmax}  |BGgrid(0=no,1=yes) |" + '\n')
    for i in range(nsect):
    	f.write('%02d'%(i+1) + "   ")
    	f.write(str(airfoiltype) + "  ")
    	f.write(('%2.2f'%25.000) + "  ")
    	f.write(('%2.2f'%00.000) + "  ")
    	f.write(('%2.2f'%00.300) + "  ")
    	f.write(('%2.2f'%00.010) + "  ")
    	f.write(('%2.2f'%00.010) + "  ")
    	f.write(('%02d'%15) + "  ")
    	f.write(('%02d'%10) + "  ")
    	f.write(('%01d'%0)  + '\n')
    f.write('\n')
    f.write(" Control table for blending section variable:" + '\n')
    f.write("           5           0           0" + '\n')
    f.write("       span                       bf1         bf2" + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "           ")
    	f.write('%01d'%1  + "           ")
    	f.write('%01d'%0  + '\n')
    f.write('\n')
    f.write(" Stacking axis location(200=centroid):" + '\n')
    f.write("   " + str('100000')  + '\n')
    f.write('\n')
    f.write(" Control points for delta_m:" + '\n')
    f.write("           " + str(5)  + '\n')
    f.write(str('        span                   delta_m')  + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
    	f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for delta_theta:" + '\n')
    f.write("           " + str(5)  + '\n')
    f.write(str('        span                   delta_theta')  + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
    	f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for in_beta*:" + '\n')
    f.write("           " + str(5)  + '\n')
    f.write(str('        span                   in_beta*')  + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
    	f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for out_beta*:" + '\n')
    f.write("           " + str(5)  + '\n')
    f.write(str('        span                   out_beta*')  + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
    	f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for chord:" + '\n')
    f.write("           " + str(5)  + '\n')
    f.write(str('        span                   chord')  + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
    	f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for tm/c:" + '\n')
    f.write("           " + str(5)  + '\n')
    f.write(str('        span                   tm/c')  + '\n')
    for i in range(5):
    	f.write("  " + ('%2.15f'%ctrlspan[i]) + "     ")
    	f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Hub offset" + '\n')
    f.write(" " + str('0.000000000000000') + '\n')
    f.write(" Tip offset" + '\n')
    f.write(" " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write("  Streamline Data" + '\n')
    f.write("  x_s      r_s" + '\n')
    for i in range(nsect):
    	for j in range(len(xsl[i])):
    		f.write(('%2.8f'%xsl[i][j]) + "  ")
    		f.write(('%2.8f'%rsl[i][j]) + '\n')
    	f.write('0 0' + '\n')
    # f.write('\n')
    # End of File
    f.close()

#------------------------------------------------------------------------------
#==============================================================================
#Creating blade file for autogrid
#==============================================================================
    fout.append(str('output_'+str(row_name[0])+'.txt'))
    with open(fout[k],'w') as output_f:
        output = subprocess.Popen(["tblade3",f.name], stdout=output_f, stderr=output_f).communicate()[0]
    subprocess.call("geomturbo "+f.name+" 241 ")
    for line in open(fout[k],'r'):
        c=0
        for word in line.lower().split():
            word = word.strip("'?,.;!-/\":")
            if "chord_actual(mm)" in line:
                if c==0:
                    data=line.lower().split()
                    chord[current_sec][k]=float(data[1])
                    c+=1
                    current_sec+=1

#------------------------------------------------------------------------------
#==============================================================================
#Calculating Design Parameters through 1-D Analysis
#==============================================================================

"At all Sections along span"
#First cell distance using Blasius solution
ywall = np.zeros((nsect,nstns))
ywall[:,0] = 6*((Vm1/0.0000157)**(-7/8))*((0.2)**(1/8));
ywall[:,1] = 6*((Vm2/0.0000157)**(-7/8))*((0.2)**(1/8));
ywall[:,2] = 6*((Vm5/0.0000157)**(-7/8))*((0.2)**(1/8));
ywall[:,3] = 6*((Vm6/0.0000157)**(-7/8))*((0.2)**(1/8));

#Solidity
sol = np.zeros((nsect,nrows))
pitch = np.zeros((nsect,nrows))
for i in range(nrows):
    pitch[:,i] = (2*np.pi*r_blde[:,i+1]*1000)/Z[i]
    sol[:,i]=chord[:,i]/pitch[:,i]

#Diffusion Factor
DiffusionFactor = np.zeros((nsect,nrows))
DiffusionFactor[:,0] = (1-W[:,1]/W[:,0])+(sw[:,1]-sw[:,0])/((r_blde[:,0]+r_blde[:,1])*sol[:,0]*W[:,0])
DiffusionFactor[:,1] = (1-V[:,2]/V[:,1])+(sw[:,2]-sw[:,1])/((r_blde[:,1]+r_blde[:,2])*sol[:,1]*V[:,1])
DiffusionFactor[:,2] = (1-W[:,3]/W[:,2])+(sw[:,3]-sw[:,2])/((r_blde[:,2]+r_blde[:,3])*sol[:,2]*W[:,2])

#----------------------------"At meanline"---------------------------------
if nsect%2==0:
    mean_num = nsect//2
else:
    mean_num = nsect//2+1

sol_m = np.zeros(nrows)
pitch_m = np.zeros(nrows)
Rx = np.zeros(nrows)
phi = np.zeros(nrows)
DeHallerNumber = np.zeros(nrows)
DiffusionFactor_m = np.zeros(nrows)
for i in range(nrows):
    pitch_m[i] = (2*np.pi*r_blde[mean_num,i+1]*1000)/Z[i]
    sol_m[i]=chord[mean_num,i]/pitch_m[i]
    Rx[i] = (np.tan(np.deg2rad(beta[mean_num,i+1]))+np.tan(np.deg2rad(beta[mean_num,i])))*Vm[mean_num,i]/(2*Ur[mean_num,i])
    phi[i] = Vm[mean_num,i+1]/Ur[mean_num,i]
    DeHallerNumber[i] = W[mean_num,i+1]/W[mean_num,i]
    DiffusionFactor_m[i] = (1-W[mean_num,i+1]/W[mean_num,i])+(sw[mean_num,i+1]-sw[mean_num,i])/((r_blde[mean_num,i]+r_blde[mean_num,i+1])*sol_m[i]*W[mean_num,i])

print("Estimate of first cell wall distance =",min(min(ywall[0]),min(ywall[1])))

fmean = open("meanlineflowproperties.dat",'w')
fmean.write("Row    Solidity    DF    DeHallerNumber      Rx      phi\n")
for i in range(nrows):
    fmean.write("%02d"%(i+1) + "    ")
    fmean.write("  "+'%2.4f'%sol_m[i]+'  '+'%2.4f'%DiffusionFactor_m[i])
    fmean.write("    " + "%2.4f"%DeHallerNumber[i])
    fmean.write("    " + "\t%2.4f"%-Rx[i])
    fmean.write("    " + "%2.4f"%phi[i])
    fmean.write('\n')
fmean.write('\n')
fmean.write("Station    Swirl   Vt[m/s]    Vm[m/s]    T[k]      Mach\n")
for i in range(nstns):
    fmean.write("%02d"%(i+1) + "  ")
    fmean.write("  "+'%2.4f'%sw[0,i]+'    '+'%2.4f'%Vt[0,i]+'    '+'%2.4f'%Vm[0,i]+'    '+'%2.4f'%T[0,i]+'    '+'%2.4f'%M[0,i]+'\n')
fmean.write('\n\n')
fmean.write("Overall Pressure ratio = %2.4f"%PR+'\n')
fmean.write("Overall Efficiency = %2.4f"%Eta+'\n')
fmean.write('\n')
fmean.write("Rotor 1 Pressure ratio = %2.4f"%PR_R1+'\n')
fmean.write("Rotor 1 Efficiency = %2.4f"%Eta_R1+'\n')
fmean.write('\n')
fmean.write("Rotor 2 Pressure ratio = %2.4f"%PR_R2+'\n')
fmean.write("Rotor 2 Efficiency = %2.4f"%Eta_R2+'\n')
fmean.write('\n')
fmean.close()
#------------------------------------------------------------------------------
#==============================================================================
#Plotting meridional view, velocity triangle and blade angles
#==============================================================================
print()
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print()
print("Generating Plots...")
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
#-----------------------------Meridional view=---------------------------------
fignum = 1
sz=15
py.figure(fignum, figsize=(16,9))
py.plot(xsl[0],rsl[0],color='g',marker='.',lw=1)
py.plot(xsl[-1],rsl[-1],color='g',marker='.',lw=1)
for i in range(0,len(r_nd)):
    lines=py.plot(x_nd[i],r_nd[i])
    py.setp(lines,color = 'k',lw=0.75)
for i in range(1,len(r_nd)-1):
    lines=py.plot(xsl[i],rsl[i])
    py.setp(lines,color = 'blue',linestyle='--',lw=0.75)
py.xlabel('x-coordinate')
py.ylabel('r-coordinate')
py.text(0.25,0.41,'R1',size=sz)
py.text(0.75,0.81,'S1',size=sz)
py.text(1,1.3,'R2',size=sz)
py.text(-0.05,1.01,'1',size=sz)
py.text(0.3,1.075,'2',size=sz)
py.text(0.37,1.1,'3',size=sz)
py.text(0.61,1.25,'4',size=sz)
py.text(0.67,1.3,'5',size=sz)
py.text(0.87,1.778,'6',size=sz)
py.text(0.92,2.65,'7',size=sz)
py.axis('equal')
py.grid(True)
py.tight_layout()
py.savefig("streamlines.png")

fignum+=1

#---------------------------Velocity triangles---------------------------------
fig = py.figure(figsize=(15, 10))
for i in range(nstns):
    ax = fig.add_subplot(2,2,i+1)
    soa = np.array([[0,0,Vm[0][i],0],[0,0,Vm[0][i],Vt[0][i]],[Vm[0][i],0,0,
                     Vt[0][i]],[0,0,Vm[0][i],Wt[0][i]],[Vm[0][i],0,0,Wt[0][i]]])
    X, Y, U, V = zip(*soa)
    ax = py.gca()
    ax.set_xlim([-10, Vm[0][i]+10])
    if Wt[0][i]<0:
        ax.set_ylim([Vt[0][i]+10, Wt[0][i]-10])
        ax.text(25,-30,r'$\beta$ = %.2f$\degree$'%(beta[0][i]),size='16')
        ax.text(25,10,r'$\alpha$ = %.2f$\degree$'%(alpha[0][i]),size='16')
    else:
        ax.set_ylim([Vt[0][i]+10, -10])
        ax.text(25,10,r'$\beta$ = %.2f$\degree$'%(beta[0][i]),size='16')
        ax.text(25,75,r'$\alpha$ = %.2f$\degree$'%(alpha[0][i]),size='16')
    ax.invert_yaxis()
    ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
    ax.set_title("Station "+str(i+1))
py.tight_layout()
fig.savefig("veltri.png")

fignum+=1

#----------------------------Blade Angles--------------------------------------
py.figure(fignum, figsize=(16,9))
py.plot(beta[:,0],span,'k',marker='.',label=r'$\beta_{in}$')
py.plot(beta[:,1],span,'r',marker='.',label=r'$\beta_{out}$')
py.legend()
py.ylabel('span',size='24')
py.xlabel(r'$\beta_z$',size='24')
leg = py.gca().get_legend()
ltext = leg.get_texts()
py.setp(ltext,fontsize='16')
py.tight_layout()
py.savefig("rotor1.png")

fignum+=1

py.figure(fignum, figsize=(16,9))
py.plot(alpha[:,1],span,'k',marker='.',label=r'$\alpha_{in}$')
py.plot(alpha[:,2],span,'r',marker='.',label=r'$\alpha_{out}$')
py.ylabel('span',size='24')
py.xlabel(r'$\alpha_z$',size='24')
py.legend()
leg = py.gca().get_legend()
ltext = leg.get_texts()
py.setp(ltext,fontsize='16')
py.tight_layout()
py.savefig("stator1.png")

fignum+=1
py.figure(fignum, figsize=(16,9))
py.plot(beta[:,2],span,'k',marker='.',label=r'$\beta_{in}$')
py.plot(beta[:,3],span,'r',marker='.',label=r'$\beta_{out}$')
py.ylabel('span',size='24')
py.xlabel(r'$\beta_r$',size='24')
py.legend()
leg = py.gca().get_legend()
ltext = leg.get_texts()
py.setp(ltext,fontsize='16')
py.tight_layout()
py.savefig("rotor2.png")
#py.show()
