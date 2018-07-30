"""Multi-rotor compressor. All units are in SI system."""
import matplotlib.pyplot as py
import numpy as np
import subprocess, timeit, os
from inputs import *
from analysis_plots import *
from design import *

start = timeit.default_timer()

workdir = os.getcwd()
#==============================================================================
# Calculating Thermodynamic Properties at inlet and outlet
#==============================================================================
"""The flow assumed here is isentropic flow. Inlet relative angles at hub and
Pressure ratio is 4.7 and efficiency is calculated based upon inputs. The air
considered is a perfect gas."""

# Evaluating based on known parameters
w = 2 * np.pi * N / 60  # Angular velocity[rad/s]
g = 1 / (1 - Rgas / Cp)  # Ratio of Specific heats[-]


# 0-D Calculations Thremodynamic properties
delTT_R1 = WorkRatio_R1 * delTT
delTT_R2 = delTT - delTT_R1

T06 = delTT + T01  # Outlet Total Temperature[K]
P06 = PR * P01  # Outlet Total pressure[Pa]
Eta = ((P06 / P01)**((g - 1) / g) - 1) / ((T06 / T01) - 1)

T05 = T06 - delTT_R2
#------------------------------------------------------------------------------
#==============================================================================
# Meanline and Velocty Triangles Calculations
#==============================================================================
#-----------------------------------Inlet--------------------------------------
U1 = w * r1m
sw1 = r1 * Vt1
Wt1 = Vt1 - U1
beta1 = np.average(Beta1_Blade)
Vm1 = Wt1 / np.tan(np.radians(beta1))
rho1 = mdot / (Vm1 * A1)
V1 = (Vt1**2 + Vm1**2)**0.5
T1 = T01 - V1**2 / (2 * Cp)
P1 = P01 * ((T1 / T01)**((g / (g - 1))))
W1 = Wt1 / np.sin(np.radians(beta1))
a1 = (g * Rgas * T1)**0.5
M1 = V1 / a1
alpha1 = np.degrees(np.arctan(Vt1 / Vm1))
M1rel = W1 / a1
T01rel = T1 + (W1**2 / (2 * Cp))
P01rel = P1 + (0.5 * rho1 * W1**2)

#------------------------------------Station 2---------------------------------
U2 = w * r2m
T02 = delTT_R1 + T01
Vt2 = (Cp * (T02 - T01) + U1 * Vt1) / U2
PR_R1 = (Eta_R1 * (T02 / T01 - 1) + 1)**(g / (g - 1))
Wt2 = Vt2 - U2
P02 = PR_R1 * P01
rho02 = P02 / (Rgas * T02)
rho2 = rho02
error = 5
while abs(error) > 1e-6:
    Vm2 = mdot / (rho2 * A2)
    V2 = (Vm2**2 + Vt2**2)**0.5
    T2 = T02 - V2**2 / (2 * Cp)
    P2 = P02 * ((T2 / T02)**((g / (g - 1))))
    rho2p = P2 / (Rgas * T2)
    error = (1 - rho2p / rho2) * 100
    rho2 = rho2p
a2 = (g * Rgas * T2)**0.5
Wm2 = Vm2
W2 = (Wt2**2 + Wm2**2)**0.5
M2 = V2 / a2
M2rel = W2 / a2
beta2 = np.degrees(np.arctan(Wt2 / Vm2))
alpha2 = np.degrees(np.arctan(Vt2 / Vm2))
T02rel = T2 + (W2**2 / (2 * Cp))
P02rel = P2 + (0.5 * rho2 * W2**2)
sw2 = r2 * Vt2

#-----------------------------R2 Inlet(Station 2)------------------------------
U5 = w * r5m
T05 = T02
P05 = P02 - Y * (P02 - P2)

alpha5 = alpha2 - dalpha
rho05 = P05 / (Rgas * T05)
rho5 = rho05
error = 5
while abs(error) > 1e-6:
    Vm5 = mdot / (rho5 * A5)
    Vt5 = np.tan(np.deg2rad(alpha5)) * Vm5
    V5 = (Vm5**2 + Vt5**2)**0.5
    T5 = T05 - V5**2 / (2 * Cp)
    P5 = P05 * ((T5 / T05)**((g / (g - 1))))
    rho5p = P5 / (Rgas * T5)
    error = (1 - rho5p / rho5) * 100
    rho5 = rho5p
a5 = (g * Rgas * T5)**0.5
Wt5 = Vt5 - U5
Eta_R2 = ((P06 / P05)**((g - 1) / g) - 1) / ((T06 / T05) - 1)
Wm5 = Vm5
W5 = (Wt5**2 + Wm5**2)**0.5
M5 = V5 / a5
M5rel = W5 / a5
beta5 = np.degrees(np.arctan(Wt5 / Vm5))
alpha5 = np.degrees(np.arctan(Vt5 / Vm5))
T05rel = T5 + (W5**2 / (2 * Cp))
P05rel = P5 + (0.5 * rho5 * W5**2)
sw5 = r5 * Vt5

#-----------------------------R2 Outlet(Station 4)-----------------------------
U6 = w * r6m
Vt6 = (Cp * (T06 - T05) + U5 * Vt5) / U6
Wt6 = Vt6 - U6
beta6 = Beta6_Blade
Vm6 = Wt6 / np.tan(np.deg2rad(beta6))
rho06 = P06 / (Rgas * T06)
rho6 = mdot / (A6 * Vm6)
alpha6 = np.rad2deg(np.arctan(Vt6 / Vm6))
Wm6 = Vm6
V6 = (Vt6**2 + Vm6**2)**0.5
W6 = (Wt6**2 + Wm6**2)**0.5
T6 = T06 - V6**2 / (2 * Cp)
P6 = P06 * ((T6 / T06)**((g / (g - 1))))
a6 = (g * Rgas * T6)**0.5
M6 = V6 / a6
M6rel = W6 / a6
T06rel = T6 + (W6**2 / (2 * Cp))
P06rel = P6 + (0.5 * rho6 * W6**2)
sw6 = r6 * Vt6
PR_R2 = P06 / P05
#-------------------------------Flow properties--------------------------------

span = np.linspace(0, 1, nsect)

r_1 = span * (r1[1] - r1[0]) + r1[0]
r_2 = span * (r2[1] - r2[0]) + r2[0]
r_5 = span * (r5[1] - r5[0]) + r5[0]
r_6 = span * (r6[1] - r6[0]) + r6[0]
r_blde = np.transpose(np.array((r_1, r_2, r_5, r_6)))

Ur = w * r_blde

Vm_s[:, 0] = Vm1
Vm_s[:, 1] = Vm2
Vm_s[:, 2] = Vm5
Vm_s[:, 3] = Vm6

Vt_s[:, 0] = r1m * Vt1 / r_1
Vt_s[:, 1] = r2m * Vt2 / r_2
Vt_s[:, 2] = r5m * Vt5 / r_5
Vt_s[:, 3] = r6m * Vt6 / r_6

sw_s = r_blde * Vt_s
Wt_s = Vt_s - Ur
beta_s = np.degrees(np.arctan(Wt_s / Vm_s))
alpha_s = np.degrees(np.arctan(Vt_s / Vm_s))
V_s = (Vt_s**2 + Vm_s**2)**0.5
W_s = (Wt_s**2 + Vm_s**2)**0.5

T_s[:, 0] = T1
T_s[:, 1] = T2
T_s[:, 2] = T5
T_s[:, 3] = T6

a = (g * T_s * Rgas)**0.5
M = V_s / a
Mrel = W_s / a

#------------------------------------------------------------------------------
#==============================================================================
# Creating T-Blade3 input file....
#==============================================================================
#-------------------------------Inputs to T-blade3-----------------------------
cntr = 0
stage_num = 1
fout = []
data = []
chord = np.zeros((nsect, nrows))
ctrlspan = np.linspace(0, 1, nsect)


def rowname(row_num, stage_num):
    if row_num % 2 == 0:
        row_name = "R" + str(stage_num)
    if row_num % 2 == 1:
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
            f.write(('%2.8f' % beta_s[i][k]) + "  ")
            f.write(('%2.8f' % beta_s[i][k + 1]) + "  ")
            f.write(('%2.8f' % Mrel[i][k]) + "  ")
        if k % 2 == 1:
            f.write(('%2.8f' % alpha_s[i][k]) + "  ")
            f.write(('%2.8f' % alpha_s[i][k + 1]) + "  ")
            f.write(('%2.8f' % M[i][k]) + "  ")
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
        f.write("    " + ('%2.8f' % x_nd[cntr][i]) + "  ")
        f.write(('%2.8f' % r_nd[cntr][i]) + "  ")
        f.write(('%2.8f' % x_nd[cntr + 1][i]) + "  ")
        f.write(('%2.8f' % r_nd[cntr + 1][i]) + "\n")
    cntr += 2
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
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "           ")
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
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for delta_theta:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   delta_theta') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for in_beta*:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   in_beta*') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for out_beta*:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   out_beta*') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for chord:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   chord') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "     ")
        f.write("  " + str('0.000000000000000') + '\n')
    f.write('\n')
    f.write(" Control points for tm/c:" + '\n')
    f.write("           " + str(5) + '\n')
    f.write(str('        span                   tm/c') + '\n')
    for i in range(5):
        f.write("  " + ('%2.15f' % ctrlspan[i]) + "     ")
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

#------------------------------------------------------------------------------
#==============================================================================
# Creating blade file for autogrid
#==============================================================================
    fout.append(str('output_' + str(row_name[0]) + '.txt'))
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
                    chord[current_sec][k] = float(data[1])
                    c += 1
                    current_sec += 1

#------------------------------------------------------------------------------
#==============================================================================
# Calculating Design Parameters through 1-D Analysis
#==============================================================================

"At all Sections along span"
# First cell distance using Blasius solution
ywall_s[:, 0] = 6 * ((Vm1 / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
ywall_s[:, 1] = 6 * ((Vm2 / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
ywall_s[:, 2] = 6 * ((Vm5 / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
ywall_s[:, 3] = 6 * ((Vm6 / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));

# Solidity
for i in range(nrows):
    pitch_s[:, i] = (2 * np.pi * r_blde[:, i + 1] * 1000) / Z[i]
    sol_s[:, i] = chord[:, i] / pitch_s[:, i]

# Diffusion Factor
#Df[:, 0] = (1 - W_s[:, 1] / W_s[:, 0]) + (sw_s[:, 1] - sw_s[:, 0]) / \
#    ((r_blde[:, 0] + r_blde[:, 1]) * sol_s[:, 0] * W_s[:, 0])
#Df[:, 1] = (1 - V_s[:, 2] / V_s[:, 1]) + (sw_s[:, 2] - sw_s[:, 1]) / \
#    ((r_blde[:, 1] + r_blde[:, 2]) * sol_s[:, 1] * V_s[:, 1])
#Df[:, 2] = (1 - W_s[:, 3] / W_s[:, 2]) + (sw_s[:, 3] - sw_s[:, 2]) / \
#    ((r_blde[:, 2] + r_blde[:, 3]) * sol_s[:, 2] * W_s[:, 2])

#----------------------------"At meanline"---------------------------------
if nsect % 2 == 0:
    mean_num = nsect // 2
else:
    mean_num = nsect // 2 + 1


for i in range(nrows):
    pitch[i] = (2 * np.pi * r_blde[mean_num, i + 1] * 1000) / Z[i]
    sol[i] = chord[mean_num, i] / pitch[i]
    Rx[i] = (np.tan(np.deg2rad(beta_s[mean_num, i + 1])) +
             np.tan(np.deg2rad(beta_s[mean_num, i]))) * Vm_s[mean_num, i] / (2 * Ur[mean_num, i])
    phi[i] = Vm_s[mean_num, i + 1] / Ur[mean_num, i]
    DH[i] = W_s[mean_num, i + 1] / W_s[mean_num, i]
    Df[i] = (1 - W_s[mean_num, i + 1] / W_s[mean_num, i]) + (sw_s[mean_num, i + 1] -
                                                                        sw_s[mean_num, i]) / ((r_blde[mean_num, i] + r_blde[mean_num, i + 1]) * sol[i] * W_s[mean_num, i])

print("Estimate of first cell wall distance =",
      min(min(ywall_s[0]), min(ywall_s[1])))

fmean = open("meanlineflowproperties.dat", 'w')
fmean.write("Row    Solidity    DF    DeHallerNumber      Rx      phi\n")
for i in range(nrows):
    fmean.write("%02d" % (i + 1) + "    ")
    fmean.write("  " + '%2.4f' %
                sol[i] + '  ' + '%2.4f' % Df[i])
    fmean.write("    " + "%2.4f" % DH[i])
    fmean.write("    " + "\t%2.4f" % -Rx[i])
    fmean.write("    " + "%2.4f" % phi[i])
    fmean.write('\n')
fmean.write('\n')
fmean.write("Station    Swirl   Vt[m/s]    Vm[m/s]    T[k]      Mach\n")
for i in range(nstns):
    fmean.write("%02d" % (i + 1) + "  ")
    fmean.write("  " + '%2.4f' % sw_s[0, i] + '    ' + '%2.4f' % Vt_s[0, i] + '    ' + '%2.4f' %
                Vm_s[0, i] + '    ' + '%2.4f' % T_s[0, i] + '    ' + '%2.4f' % M[0, i] + '\n')
fmean.write('\n\n')
fmean.write("Overall Pressure ratio = %2.4f" % PR + '\n')
fmean.write("Overall Efficiency = %2.4f" % Eta + '\n')
fmean.write('\n')
fmean.write("Rotor 1 Pressure ratio = %2.4f" % PR_R1 + '\n')
fmean.write("Rotor 1 Efficiency = %2.4f" % Eta_R1 + '\n')
fmean.write('\n')
fmean.write("Rotor 2 Pressure ratio = %2.4f" % PR_R2 + '\n')
fmean.write("Rotor 2 Efficiency = %2.4f" % Eta_R2 + '\n')
fmean.write('\n')
fmean.close()


plots(xsl, rsl, x_nd, r_nd, Vm_s, Vt_s, W_s, Wt_s, alpha_s, beta_s, span, nstns)

stop = timeit.default_timer()
print(" Execution Time: ", '%1.3f'%(stop - start) , "seconds")
