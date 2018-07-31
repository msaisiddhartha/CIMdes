"""Multi-rotor compressor. All units are in SI system."""
import matplotlib.pyplot as py
import numpy as np
import subprocess
import timeit
import os
from design import *
from inputs import *
from loss import *
from functions import *


from analysis_plots import *


start = timeit.default_timer()

workdir = os.getcwd()

rm, area, r_s, bsf, xsl, rsl, bw = streamlines()

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
sw[0] = rm[0] * Vt[0]

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
sw[1] = rm[1] * Vt[1]

#---------------------------S1 inlet(Station 3)---------------------------------
U[2] = w * rm[2]
T0[2] = T0[1]
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
U[3] = w * rm[3]
T0[3] = T0[1]
P0[3] = P0[1] - Y * (P0[1] - P[1])

alpha[3] = alpha[1] - dalpha
rho05 = P0[3] / (Rgas * T0[3])
rho[3] = rho05
error = 5
while abs(error) > 1e-6:
    Vm[3] = mdot / (rho[3] * area[3])
    Vt[3] = np.tan(np.deg2rad(alpha[3])) * Vm[3]
    V[3] = (Vm[3]**2 + Vt[3]**2)**0.5
    T[3] = T0[3] - V[3]**2 / (2 * Cp)
    P[3] = P0[3] * ((T[3] / T0[3])**((g / (g - 1))))
    rho5p = P[3] / (Rgas * T[3])
    error = (1 - rho5p / rho[3]) * 100
    rho[3] = rho5p
a[3] = (g * Rgas * T[3])**0.5
Wt[3] = Vt[3] - U[3]
Eta_R2 = ((P0[5] / P0[3])**((g - 1) / g) - 1) / ((T0[5] / T0[3]) - 1)
Wm[3] = Vm[3]
W[3] = (Wt[3]**2 + Wm[3]**2)**0.5
M[3] = V[3] / a[3]
Mrel[3] = W[3] / a[3]
beta[3] = np.degrees(np.arctan(Wt[3] / Vm[3]))
alpha[3] = np.degrees(np.arctan(Vt[3] / Vm[3]))
T0rel[3] = T[3] + (W[3]**2 / (2 * Cp))
P0rel[3] = P[3] + (0.5 * rho[3] * W[3]**2)
sw[3] = rm[3] * Vt[3]

#---------------------------S1 inlet(Station 3)---------------------------------
U[4] = w * rm[4]
T0[4] = T0[3]
Vt[4] = Vt[3]
Wt[4] = Wt[3]
P0[4] = P0[3]
rho[4] = rho[3]
Vm[4] = Vm[3]
V[4] = V[3]
T[4] = T[3]
P[4] = P[3]
a[4] = a[3]
Wm[4] = Wm[3]
W[4] = W[3]
M[4] = M[3]
Mrel[4] = Mrel[3]
beta[4] = beta[3]
alpha[4] = alpha[3]
T0rel[4] = T0rel[3]
P0rel[4] = P0rel[3]
sw[4] = sw[3]

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
sw[5] = rm[5] * Vt[5]
#-------------------------------Flow properties--------------------------------

r_span = np.zeros((nsect, nstations))
span = np.linspace(0, 1, nsect)
for i in range(nstations):
    for j in range(nsect):
        r_span[j, i] = span[j] * (r_s[i, 1] - r_s[i, 0]) + r_s[i, 0]

Ur = w * r_span

#free-vortex calculation along span
for i in range(nstations):
    Vm_s[:, i] = Vm[i]
    T_s[:, i] = T[i]
    for j in range(nsect):
        Vt_s[j, i] = rm[i] * Vt[i] / r_span[j, i]

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
# Calculating Design Parameters through 1-D Analysis
#==============================================================================

"At all Sections along span"
# First cell distance using Blasius solution
ywall_s = 6 * ((Vm / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8))
print("Estimate of first cell wall distance =", np.amin(ywall_s))

#----------------------------"At meanline"---------------------------------
if nsect % 2 == 0:
    mean_num = nsect // 2
else:
    mean_num = nsect // 2 + 1
stagenum = 1

keyword = "chord_actual(mm):"
c = 0
for i in range(nrows):
    stagenum = (i // 2) + 1
    if (i + 1) % 2 == 1:
        row_name = "R"
    if (i + 1) % 2 == 0:
        row_name = "S"
    create_tblade3(i, c, row_name, stagenum, data, nsect, bsf, beta_s, alpha_s, Mrel_s, M_s, x_s, r_s, span ,xsl, rsl)
    flog = open('output_' + row_name + str(stagenum) + '.txt', 'r')
    chord[:,i] = chord_lookup(flog, row_name, stagenum, keyword)
    pitch[i] = 2 * np.pi * rm[c + 1] / Z[i]
    sol[i] = chord[mean_num, i] / pitch[i]
    phi[i] = FlowCoeff(Vm[c+1], U[c+1])
#Checking for axial or not
    if rm[c]/rm[c+1] >=0.7:
        axial = True
    if rm[c]/rm[c+1] <=0.7:
        axial = False

#Check for rotor or stator
    if row_name == "R":
        #Rx[i] = DegofReac(P[c + 2], P[c + 1], P[c])
        DH[i] = dehaller(W[c+1], W[c])
        Re = Recalc(chord[mean_num, i], rho[c], W[c])
        Df[i] = DiffusionFact(Cp, sw[i], DH[i], Wt[c+1], Wt[c], W[c], sol[i], W[c+1], rm[c], W_s[-1, c], U[c], T0[c+1], T0[c],  r_s[c, -1], axial)
    if row_name == "S":
        DH[i] = dehaller(V[c+1], V[c])
        Re = Recalc(chord[mean_num, i], rho[c], V[c])
        Df[i] = DiffusionFact(Cp, sw[i], DH[i], Vt[c+1], Vt[c], V[c], sol[i], V[c+1], rm[c], V_s[-1, c], U[c], T0[c+1], T0[c],  r_s[c, -1], axial)

    c += 2

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
                Vm_s[0, i] + '    ' + '%2.4f' % T_s[0, i] + '    ' + '%2.4f' % M_s[0, i] + '\n')
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


plots(xsl, rsl, Vm_s, Vt_s, W_s, Wt_s, alpha_s, beta_s, span, nstns, bsf)

stop = timeit.default_timer()
print(" Execution Time: ", '%1.3f' % (stop - start), "seconds")
