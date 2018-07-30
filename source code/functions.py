import numpy as np
import subprocess, os
from  inputs import *
from design import *

#==============================================================================
# Calculating Thermodynamic Properties at inlet and outlet
#==============================================================================
"""The flow assumed here is isentropic flow. Inlet relative angles at hub and
Pressure ratio is 4.7 and efficiency is calculated based upon inputs. The air
considered is a perfect gas."""

# Evaluating based on known parameters
w = 2 * np.pi * N / 60  # Angular velocity[rad/s]
g = 1 / (1 - Rgas / Cp)  # Ratio of Specific heats[-]
Vt[0] = Vt1
T0[0] = T01
P0[0] = P01

# 0-D Calculations Thremodynamic properties
delTT_R1 = WorkRatio_R1 * delTT
delTT_R2 = delTT - delTT_R1

T0[5] = delTT + T0[0]  # Outlet Total Temperature[K]
P0[5] = PR * P0[0]  # Outlet Total pressure[Pa]
Eta = ((P0[5] / P0[0])**((g - 1) / g) - 1) / ((T0[5] / T0[0]) - 1)

T0[4] = T0[5] - delTT_R2
rm, area, r_s, bsf, xsl, rsl = streamlines()
#-----------------------------------Inlet--------------------------------------

U[0] = w * rm[0]
Wt[0] = Vt[0] - U[0]
beta[0] = np.average(Beta1_Blade)
Vm[0] = Wt[0] / np.tan(np.radians(beta[0]))
rho[0] = mdot / (Vm[0] * area[0])
V[0] = (Vt[0]**2 + Vm[0]**2)**0.5
T[0] = T0[0] - V[0]**2 / (2 * Cp)
P[0] = P0[0] * ((T[0] / T0[0])**((g / (g - 1))))
W[0] = Wt[0] / np.sin(np.radians(beta[0]))
a[0] = (g * Rgas * T[0])**0.5
M[0] = V[0] / a[0]
alpha[0] = np.degrees(np.arctan(Vt[0] / Vm[0]))
Mrel[0] = W[0] / a[0]
T0rel[0] = T[0] + (W[0]**2 / (2 * Cp))
P0rel[0] = P[0] + (0.5 * rho[0] * W[0]**2)
sw[0] = rm[0]*Vt[0]
#------------------------------------Station 2---------------------------------
U[1] = w * rm[1]
T0[1] = delTT_R1 + T0[0]
Vt[1] = (Cp * (T0[1] - T0[0]) + U[0] * Vt[0]) / U[1]
PR_R1 = (Eta_R1 * (T0[1] / T0[0] - 1) + 1)**(g / (g - 1))
Wt[1] = Vt[1] - U[1]
P0[1] = PR_R1 * P0[0]
rho02 = P0[1] / (Rgas * T0[1])
rho[1] = rho02
error = 5
while abs(error) > 1e-6:
    Vm[1] = mdot / (rho[1] * area[1])
    V[1] = (Vm[1]**2 + Vt[1]**2)**0.5
    T[1] = T0[1] - V[1]**2 / (2 * Cp)
    P[1] = P0[1] * ((T[1] / T0[1])**((g / (g - 1))))
    rho2p = P[1] / (Rgas * T[1])
    error = (1 - rho2p / rho[1]) * 100
    rho[1] = rho2p
a[1] = (g * Rgas * T[1])**0.5
Wm[1] = Vm[1]
W[1] = (Wt[1]**2 + Wm[1]**2)**0.5
M[1] = V[1] / a[1]
Mrel[1] = W[1] / a[1]
beta[1] = np.degrees(np.arctan(Wt[1] / Vm[1]))
alpha[1] = np.degrees(np.arctan(Vt[1] / Vm[1]))
T0rel[1] = T[1] + (W[1]**2 / (2 * Cp))
P0rel[1] = P[1] + (0.5 * rho[1] * W[1]**2)
sw[1] = rm[1]*Vt[1]

#---------------------------S1 inlet(Station 3)---------------------------------
U[2] = w*rm[2]
T0[2]= T0[1]
Vt[2] = Vt[1]
Wt[2] = Wt[1]
P0[2] = P0[1]
rho[2] = rho[1]
Vm[2] = Vm[1]
V[2] = V[1]
T[2] = T[1]
P[2] = P[1]
a[2] = a[1]
Wm[2] = Wm[1]
W[2] = W[1]
M[2] = M[1]
Mrel[2] = Mrel[1]
beta[2] = beta[1]
alpha[2] = alpha[1]
T0rel[2] = T0rel[1]
P0rel[2] = P0rel[1]
sw[2] = sw[1]


#-----------------------------R2 Inlet(Station 5)------------------------------
U[4] = w * rm[4]
T0[4] = T0[1]
P0[4] = P0[1] - Y * (P0[1] - P[1])

alpha[4] = alpha[1] - dalpha
rho05 = P0[4] / (Rgas * T0[4])
rho[4] = rho05
error = 5
while abs(error) > 1e-6:
    Vm[4] = mdot / (rho[4] * area[4])
    Vt[4] = np.tan(np.deg2rad(alpha[4])) * Vm[4]
    V[4] = (Vm[4]**2 + Vt[4]**2)**0.5
    T[4] = T0[4] - V[4]**2 / (2 * Cp)
    P[4] = P0[4] * ((T[4] / T0[4])**((g / (g - 1))))
    rho5p = P[4] / (Rgas * T[4])
    error = (1 - rho5p / rho[4]) * 100
    rho[4] = rho5p
a[4] = (g * Rgas * T[4])**0.5
Wt[4] = Vt[4] - U[4]
Eta_R2 = ((P0[5] / P0[4])**((g - 1) / g) - 1) / ((T0[5] / T0[4]) - 1)
Wm[4] = Vm[4]
W[4] = (Wt[4]**2 + Wm[4]**2)**0.5
M[4] = V[4] / a[4]
Mrel[4] = W[4] / a[4]
beta[4] = np.degrees(np.arctan(Wt[4] / Vm[4]))
alpha[4] = np.degrees(np.arctan(Vt[4] / Vm[4]))
T0rel[4] = T[4] + (W[4]**2 / (2 * Cp))
P0rel[4] = P[4] + (0.5 * rho[4] * W[4]**2)
sw[4] = rm[4]*Vt[4]

#---------------------------S1 inlet(Station 4)---------------------------------
U[3] = w*rm[3]
T0[3]= T0[4]
Vt[3] = Vt[4]
Wt[3] = Wt[4]
P0[3] = P0[4]
rho[3] = rho[4]
Vm[3] = Vm[4]
V[3] = V[4]
T[3] = T[4]
P[3] = P[4]
a[3] = a[4]
Wm[3] = Wm[4]
W[3] = W[4]
M[3] = M[4]
Mrel[3] = Mrel[4]
beta[3] = beta[4]
alpha[3] = alpha[4]
T0rel[3] = T0rel[4]
P0rel[3] = P0rel[4]
sw[3] = sw[4]
#-----------------------------R2 Outlet(Station 6)-----------------------------
U[5] = w * rm[5]
Vt[5] = (Cp * (T0[5] - T0[4]) + U[4] * Vt[4]) / U[5]
Wt[5] = Vt[5] - U[5]
beta[5] = Beta6_Blade
Vm[5] = Wt[5] / np.tan(np.deg2rad(beta[5]))
rho06 = P0[5] / (Rgas * T0[5])
rho[5] = mdot / (area[5] * Vm[5])
alpha[5] = np.rad2deg(np.arctan(Vt[5] / Vm[5]))
Wm[5] = Vm[5]
V[5] = (Vt[5]**2 + Vm[5]**2)**0.5
W[5] = (Wt[5]**2 + Wm[5]**2)**0.5
T[5] = T0[5] - V[5]**2 / (2 * Cp)
P[5] = P0[5] * ((T[5] / T0[5])**((g / (g - 1))))
a[5] = (g * Rgas * T[5])**0.5
M[5] = V[5] / a[5]
Mrel[5] = W[5] / a[5]
T0rel[5] = T[5] + (W[5]**2 / (2 * Cp))
P0rel[5] = P[5] + (0.5 * rho[5] * W[5]**2)
PR_R2 = P0[5] / P0[4]
sw[5] = rm[5]*Vt[5]
#-------------------------------Flow properties--------------------------------

r_span = np.zeros((nsect, nstations))
span = np.linspace(0, 1, nsect)
for i in range(nstations):
    for j in range(nsect):
        r_span[j, i] = span[j] * (r_s[i, 1] - r_s[i, 0]) + r_s[i, 0]

Ur = w * r_span

for i in range(nstations):
    Vm_s[:, i] = Vm[i]
    T_s[:, i] = T[i]
    for j in range(nsect):
        Vt_s[j, i] = rm[i] * Vt[i] / r_span[j,i]

sw_s = r_span * Vt_s
Wt_s = Vt_s - Ur
beta_s = np.degrees(np.arctan(Wt_s / Vm_s))
alpha_s = np.degrees(np.arctan(Vt_s / Vm_s))
V_s = (Vt_s**2 + Vm_s**2)**0.5
W_s = (Wt_s**2 + Vm_s**2)**0.5
a_s = (g * T_s * Rgas)**0.5
M_s = V_s / a_s
Mrel_s = W_s / a_s

#------------------------------------------------------------------------------
#==============================================================================
# Creating T-Blade3 input file....
#==============================================================================
#-------------------------------Inputs to T-blade3-----------------------------

print()
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print()
print("Creating geomturbo file...")
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

def rowname(row_num, stage_num):
    if row_num % 2 == 1:
        row_name = "R" + str(stage_num)
    if row_num % 2 == 0:
        row_name = "S" + str(stage_num)
        stage_num += 1
    return row_name, stage_num

# Flow angle switch


def angles(curr_row):
    if curr_row == 1:
        return 0
    if curr_row == 2:
        ang = 2
    else:
        ang = 1
    return ang

def inc_angle(k):
    IncAngR1LE = np.zeros(nsect)
    if k == 0:
        IncAngR1LE[0] = beta[0][0] - beta1b[0]
        IncAngR1LE[-1] = beta[-1][0] - beta1b[-1]

    # Assuming linear variation from hub to shroud
        IncAngR1LE = np.linspace(IncAngR1LE[0], IncAngR1LE[-1], nsect)

    return IncAngR1LE


def create_tblade3(k, cntr, row_name, stagenum, fout, data):
    current_sec = 0
    f = open("tblade3input." + str(k + 1) + "." + str(casename) + ".dat", "w")
    f.write("Input parameters (version 1.1)" + '\n')
    f.write("    " + str(row_name[0]) + '\n')
    f.write(" Blade row #:" + '\n')
    f.write("    " + str(k + 1) + '\n')
    f.write(" Number of blades in this row:" + '\n')
    f.write("    " + str(Z[k]) + '\n')
    f.write(" Blade Scaling factor (mm):" + '\n')
    f.write("    " + str(bsf * 1000) + '\n')
    f.write(" Number of streamlines:" + '\n')
    f.write("    " + str(nsect) + '\n')
    f.write(" Angles in the input file (0=Beta_z (default),1=Beta_r):" + '\n')
    f.write("    " + str(angles(k + 1)) + '\n')
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
        f.write(('%02d' % (i + 1)) + "   ")
        if k % 2 == 0:
            f.write(('%2.8f' % beta_s[i, cntr]) + "  ")
            f.write(('%2.8f' % beta_s[i, cntr + 1]) + "  ")
            f.write(('%2.8f' % Mrel_s[i, cntr]) + "  ")
        if k % 2 == 1:
            f.write(('%2.8f' % alpha_s[i, cntr]) + "  ")
            f.write(('%2.8f' % alpha_s[i, cntr + 1]) + "  ")
            f.write(('%2.8f' % M_s[i, cntr]) + "  ")
        f.write(('%2.8f' % chrdr_nd) + "  ")
        f.write(('%2.8f' % 0.0150) + "  ")
        f.write(('%2.8f' % 0.0000) + "  ")
        f.write(('%2.8f' % 0.0000) + "  ")
        f.write(('%2.8f' % 0.0000) + "\n")
    f.write('\n')
    f.write(" LE / TE curve (x,r) definition :" + '\n')
    f.write(" Number of Curve points :" + '\n')
    f.write("    " + str('2') + '\n')
    f.write("   xLE          rLE           xTE          rTE" + '\n')
    for i in range(2):
        f.write("    " + ('%2.8f' % (x_s[cntr, i]/bsf)) + "  ")
        f.write(('%2.8f' % (r_s[cntr, i]/bsf)) + "  ")
        f.write(('%2.8f' % (x_s[cntr + 1, i])) + "  ")
        f.write(('%2.8f' % (r_s[cntr + 1, i])) + "\n")
    f.write('\n')
    f.write(" # Airfoil type and Variable Radial Stacking information.         #" + '\n')
    f.write(" # stack_u: % chord stack (0.00 to 100.00).                       #" + '\n')
    f.write(" # stack_v: % below or above meanline stack (-100.00 to +100.00). #" + '\n')
    f.write(" # Use +200 for stacking on airfoil area centroid.                #" + '\n')
    f.write(" Variable Radial stacking (0=no,1=yes):" + '\n')
    f.write("    " + ('%01d' % 0) + '\n')
    f.write(
        " J   type |stk_u |stk_v |umxthk |lethk |tethk  |Jcells(Grid:4n+1) |eta_ofst(<=10){%thkc/Jmax}  |BGgrid(0=no,1=yes) |" + '\n')
    for i in range(nsect):
        f.write('%02d' % (i + 1) + "   ")
        f.write(str(airfoiltype) + "  ")
        f.write(('%2.2f' % 25.000) + "  ")
        f.write(('%2.2f' % 00.000) + "  ")
        f.write(('%2.2f' % 00.300) + "  ")
        f.write(('%2.2f' % 00.010) + "  ")
        f.write(('%2.2f' % 00.010) + "  ")
        f.write(('%02d' % 15) + "  ")
        f.write(('%02d' % 10) + "  ")
        f.write(('%01d' % 0) + '\n')
    f.write('\n')
    f.write(" Control table for blending section variable:" + '\n')
    f.write("           5           0           0" + '\n')
    f.write("       span                       bf1         bf2" + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "           ")
        f.write('%01d' % 1 + "           ")
        f.write('%01d' % 0 + '\n')
    f.write('\n')
    f.write(" Stacking axis location(200=centroid):" + '\n')
    f.write("   " + str('100000') + '\n')
    f.write('\n')
    f.write(" Control points for delta_m:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   delta_m') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for delta_theta:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   delta_theta') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for in_beta*:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   in_beta*') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for out_beta*:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   out_beta*') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for chord:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   chord') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for tm/c:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   tm/c') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % span[i]) + "     ")
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
            f.write(('%2.8f' % xsl[i][j]) + "  ")
            f.write(('%2.8f' % rsl[i][j]) + '\n')
        f.write('0 0' + '\n')
    # f.write('\n')
    # End of File
    f.close()

    fout.append(str('output_' + str(row_name) + '.txt'))
    with open(fout[k], 'w') as output_f:
        output = subprocess.Popen(
            ["tblade3", f.name], stdout=output_f, stderr=output_f).communicate()[0]
    subprocess.call("geomturbo " + f.name + " 241 ")
    for line in open(fout[k], 'r'):
        c = 0
        for word in line.lower().split():
            word = word.strip("'?,.;!-/\":")
            if "chord_actual(mm)" in line:
                if c == 0:
                    data = line.lower().split()
                    chord[current_sec,k] = float(data[1])
                    c += 1
                    current_sec += 1
        return chord[:,k]
