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

rm, area, r_s, bsf, xsl, rsl, bw, gamma = streamlines()

# Evaluating based on known parameters
w = 2 * np.pi * N / 60  # Angular velocity[rad/s]
g = 1 / (1 - Rgas / Cp)  # Ratio of Specific heats[-]
U = w*rm
Vt[0] = Vt1
T0[0] = T01
P0[0] = P01

delTT_row = WorkRatio * delTT
T0[-1] = T0[0] + delTT
# 0-D Calculations Thremodynamic properties

Wt[0] = Vt[0] - U[0]
rho0 = P0[0] / (Rgas * T0[0])
rho[0], Vm[0], V[0], T[0], P[0] = ThermoPropRotor(rho0, area[0], Vt[0], T0[0], P0[0], g)
for i in range(0, nstations-2):
    if i%4==0:
        T0[i+1] = T0[i] + delTT_row[i//4]
        P0[i+1] = P0[i] * (T0[i+1]/T0[i])**(g / (g - 1))
        Vt[i+1] = (Cp * (T0[i+1] - T0[i]) + U[i] * Vt[i]) / U[i+1]
        Wt[i+1] = Vt[i+1] - U[i+1]
        rho0 = P0[i+1] / (Rgas * T0[i+1])
        rho[i+1], Vm[i+1], V[i+1], T[i+1], P[i+1] = ThermoPropRotor(rho0, area[i+1], Vt[i+1], T0[i+1], P0[i+1], g)

    if i%4 == 1 or i%4 ==3:
        T0[i+1] = T0[i]
        P0[i+1] = P0[i]
        Vt[i+1] = rm[i] * Vt[i]/rm[i+1]
        rho[i+1] = rho[i]
        Vm[i+1] = Vm[i]
        V[i+1] = V[i]
        T[i+1] = T[i]
        P[i+1]=P[i]
        Wt[i+1] = Vt[i+1] - U[i+1]
        alpham[i+1] = np.degrees(np.arctan(Vt[i+1] / Vm[i+1]))

    if i%4 == 2:
        T0[i+1] = T0[i]
        P0[i+1] = P0[i] - Y * (P0[i] - P0[i])
        alpham[i+1] = alpham[i] - dalpha[i//4]
        rho0 = P0[i+1] / (Rgas * T0[i+1])
        rho[i+1], Vm[i+1], V[i+1], T[i+1], P[i+1], Vt[i+1] = ThermoPropStator(rho0, area[i+1], alpham[i+1], T0[i+1], P0[i+1], g)
        Wt[i+1] = Vt[i+1] - U[i+1]

#-----------------------------R2 Outlet(Station 6)-----------------------------
T0[4] = T0[5] - delTT_row[1]
P0[5] = P0[4] * (T0[5]/T0[4])**(g / (g - 1))
Vt[5] = (Cp * (T0[5] - T0[4]) + U[4] * Vt[4]) / U[5]
Wt[5] = Vt[5] - U[5]
betam[5] = Beta6_Blade
Vm[5] = Wt[5] / np.tan(np.deg2rad(betam[5]))
rho06 = P0[5] / (Rgas * T0[5])
rho[5] = mdot / (area[5] * Vm[5])
alpham[5] = np.rad2deg(np.arctan(Vt[5] / Vm[5]))
V[5] = (Vt[5]**2 + Vm[5]**2)**0.5
T[5] = T0[5] - V[5]**2 / (2 * Cp)
P[5] = P0[5] * ((T[5] / T0[5])**((g / (g - 1))))

a = (g * Rgas * T)**0.5
Wm = Vm
W = (Wt**2 + Wm**2)**0.5
M = V / a
Mrel = W / a
betam = np.degrees(np.arctan(Wt / Vm))
alpham = np.rad2deg(np.arctan(Vt / Vm))
T0rel = T + (W**2 / (2 * Cp))
P0rel = P + (0.5 * rho * W**2)
sw = rm * Vt
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
betam_s = np.degrees(np.arctan(Wt_s / Vm_s))
alpha_s = np.degrees(np.arctan(Vt_s / Vm_s))
V_s = (Vt_s**2 + Vm_s**2)**0.5
W_s = (Wt_s**2 + Vm_s**2)**0.5
a_s = (g * T_s * Rgas)**0.5
M_s = V_s / a_s
Mrel_s = W_s / a_s

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
    create_tblade3(i, c, row_name, stagenum, data, nsect, bsf, betam_s, alpha_s, Mrel_s, M_s, x_s, r_s, span ,xsl, rsl)
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

#Check for rotor or stator and calcualte properties and losses
    if row_name == "R":
        DH[i] = dehaller(W[c+1], W[c])
        Re = Recalc(chord[mean_num, i], rho[c], W[c])
        Df[i] = DiffusionFact(Cp, sw[i], DH[i], Wt[c+1], Wt[c], W[c], sol[i], W[c+1], rm[c], W_s[-1, c], U[c], T0[c+1], T0[c],  r_s[c, -1], axial)
        dH_Loss[0, i] = IncLoss(Vm[c], Wt[c], betam[c])
        dH_Loss[1,i] = SkinFricLoss(W[c+1], W[c], r_s[c,0], r_s[c,1], betam[c], Z[i], Re, chord[mean_num, i])
        dH_Loss[4,i] = RecirculationLoss(betam[c+1], Df[i], U[c+1])

    if row_name == "S":
        DH[i] = dehaller(V[c+1], V[c])
        Re = Recalc(chord[mean_num, i], rho[c], V[c])
        Df[i] = DiffusionFact(Cp, sw[i], DH[i], Vt[c+1], Vt[c], V[c], sol[i], V[c+1], rm[c], V_s[-1, c], U[c], T0[c+1], T0[c],  r_s[c, -1], axial)
        dH_Loss[0, i] = IncLoss(Vm[c], Wt[c], betam[c])
        dH_Loss[1,i] = SkinFricLoss(V[c+1], V[c], r_s[c,0], r_s[c,1], alpham[c], Z[i], Re, chord[mean_num, i])
        dH_Loss[4,i] = RecirculationLoss(alpham[c+1], Df[i], U[c+1])
    #Rx[i] = DegofReac(P[c + 2], P[c + 1], P[c])
    dH_Loss[2,i] = BladeLoadLoss(Df[i], U[c+1])
    dH_Loss[3,i] = ClearanceLoss(r_s[c,0], r_s[c,1], rm[c+1], cl[i], bw[c+1], U[c+1])
    dH_Loss[5,i] = np.fabs(LeakageLoss(bw[c+1], bw[c], rm[c+1], rm[c], Vt[c+1], Vt[c], Z[i], cl[i], rho[c+1], U[c+1], chord[mean_num, i]))
    dH_Loss[6,i] = DiskFricLoss(U[c+1], rm[c+1], rho[c+1], rho[c], Re)

    c += 2

for i in range(7):
    dH_Loss[i,3] = np.sum(dH_Loss[i,0:3])

#------------------------------------------------------------------------------
#==============================================================================
# Calculating Design Parameters through 1-D Analysis
#==============================================================================

"At all Sections along span"
# First cell distance using Blasius solution
ywall_s = 6 * ((Vm / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8))
print("Estimate of first cell wall distance =", np.amin(ywall_s))

PR_r = np.zeros(3)
TR_r = np.zeros(3)
Eta_r = np.zeros(3)

dH[0,0] = Cp*(T0[1] - T0[0])
dH[0,1] = Cp*(T0[3] - T0[2])
dH[0,2] = Cp*(T0[5] - T0[4])
dH[0,3] = Cp*(T0[5] - T0[0])

dH[1,0] = np.sum(dH_Loss[0:4,0])
dH[1,1] = np.sum(dH_Loss[0:4,1])
dH[1,2] = np.sum(dH_Loss[0:4,2])
dH[1,3] = np.sum(dH_Loss[0:4,3])

dH[2,0] = np.sum(dH_Loss[4:7,0])
dH[2,1] = np.sum(dH_Loss[4:7,1])
dH[2,2] = np.sum(dH_Loss[4:7,2])
dH[2,3] = np.sum(dH_Loss[4:7,3])

Eta_r = (dH[0,:]-dH[1,:])/(dH[0,:]+dH[2,:])

PR_r[0] = P0[1]/P0[0]
TR_r[0] = T0[1]/T0[0]

PR_r[1] = P0[5]/P0[4]
TR_r[1] = T0[5]/T0[4]


TR = T0[-1]/T0[0]
PR = (1+Eta_r[-1]*(TR-1))**(g/(g-1))

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
fmean.write("Overall Efficiency = %2.4f" % Eta_r[3] + '\n')
fmean.write('\n')
fmean.write("Rotor 1 Pressure ratio = %2.4f" % PR_r[0] + '\n')
fmean.write("Rotor 1 Efficiency = %2.4f" % Eta_r[0] + '\n')
fmean.write('\n')
fmean.write("Rotor 2 Pressure ratio = %2.4f" % PR_r[2] + '\n')
fmean.write("Rotor 2 Efficiency = %2.4f" % Eta_r[2] + '\n')
fmean.write('\n')
fmean.close()


plots(xsl, rsl, Vm_s, Vt_s, W_s, Wt_s, alpha_s, betam_s, span, nstns, bsf)

stop = timeit.default_timer()
print(" Execution Time: ", '%1.3f' % (stop - start), "seconds")
