"""Multi-rotor compressor. All units are in SI system."""
#import pylab as py
import numpy as np
import subprocess
import timeit
import os
import pandas as pd
from tabulate import tabulate

from inputs import *
from design import *
from loss import *
from functions import *
from analysis_plots import *

np.seterr(divide='ignore', invalid='ignore')

start = timeit.default_timer()

workdir = os.getcwd()


xm, rm, area, x_s, r_s, bsf, xsl, rsl, bw, gamma, minimum, maximum = streamlines()
stagenum = 1
rownum = 0
rowindex = []
Eta = np.ones((nrows + 1, 1))

# Evaluating based on known parameters
w = 2 * np.pi * N / 60  # Angular velocity[rad/s]
g = 1 / (1 - Rgas / Cp)  # Ratio of Specific heats[-]
U = w * rm
Vt[0] = Vt1
T0[0] = T01
P0[0] = P01

delTT_row = WorkRatio * delTT
T0[-1] = T0[0] + delTT
dH[0, -1] = Cp * (T0[-1] - T0[0])

# Flow properties along span
r_span = np.zeros((nsect, nstations))
span = np.linspace(0, 1, nsect)
for i in range(nstations):
    for j in range(nsect):
        r_span[j, i] = span[j] * (r_s[i, 1] - r_s[i, 0]) + r_s[i, 0]

Ur = w * r_span

# Rotor-1 inlet calculations
Wt[0] = Vt[0] - U[0]
rho0 = P0[0] / (Rgas * T0[0])
M[0] = Mach(area[0], T0[0], P0[0], g)
rho[0], Vm[0], V[0], T[0], P[0] = ThermoPropRotor(
    rho0, M[0, 0], area[0], Vt[0], T0[0], P0[0], g)
a[0], Wm[0], W[0], mach2,  Mrel[0], T0rel[0], P0rel[0], betam[0], alpham[0], Vz[0], Vr[0], betaz[0], alphaz[0] = Properties_add(
    g, T[0], Vm[0], Wt[0], V[0], P[0], rho[0], Vt[0], gamma[0], Vz[0])
sw[0] = rm[0] * Vt[0]
h[0]= Cp*T[0]
h0rel[0] = h[0] + 0.5*(W[0]**2)
I[0] = h[0] + 0.5*(W[0]**2-U[0]**2)
V_s[:, 0], W_s[:, 0], betam_s[:, 0], alpham_s[:, 0], M_s[:, 0], Mrel_s[:,
                                                                       0], Wt_s[:, 0], Vt_s[:, 0] = free_vortex(0, rm, Vt, r_span, Vm, T, g, Ur)
Vr_s[:, 0] = Vr[0]
Vz_s[:, 0] = Vz[0]
Vm_s[:, 0] = Vm[0]


keyword = "chord_actual(mm):"

# Calculations from rotor-1 exit onwards
for i in range(0, nstations - 1):

    error = 5

#----------------------------Detect rotor exit----------------------------------
    if i % 4 == 0:
        # Check axial or not
        if rm[i] / rm[i + 1] >= 0.7:
            axial = True
        if rm[i] / rm[i + 1] <= 0.7:
            axial = False

        row_name = "R"
        T0[i + 1] = T0[i] + delTT_row[i // 4]
        print("Checking for R" + str(stagenum) + " efficiency convergence")
        cntr = 0
        while np.fabs(error) > 1e-6:
            # for j in range(0,1):
            P0[i + 1] = P0[i] * (1 + Eta[rownum] *
                                 ((T0[i + 1] / T0[i]) - 1))**(g / (g - 1))
            # M[i + 1] = Mach(area[i + 1], T0[0], P0[0], g)
            Vt[i + 1] = (Cp * (T0[i + 1] - T0[i]) + U[i] * Vt[i]) / U[i + 1]
            Wt[i + 1] = Vt[i + 1] - U[i + 1]
            rho0 = P0[i + 1] / (Rgas * T0[i + 1])
            rho[i + 1], Vm[i + 1], V[i + 1], T[i + 1], P[i +
                                                         1] = ThermoPropRotor(rho0, M[i + 1, 0], area[i + 1], Vt[i + 1], T0[i + 1], P0[i + 1], g)
            a[i + 1], Wm[i + 1], W[i + 1], M[i+1], Mrel[i + 1], T0rel[i + 1], P0rel[i + 1], betam[i + 1], alpham[i + 1], Vz[i + 1], Vr[i + 1], betaz[i +
                                                                                                                                             1], alphaz[i + 1] = Properties_add(g, T[i + 1], Vm[i + 1], Wt[i + 1], V[i + 1], P[i + 1], rho[i + 1], Vt[i + 1], gamma[i + 1], Vz[i + 1])
            sw[i + 1] = rm[i + 1] * Vt[i + 1]
            h[i+1] = Cp*T[i+1]
            h0rel[i+1] = h[i+1] + 0.5*(W[i+1]**2)
            I[i+1] = h[i+1] + 0.5*(W[i+1]**2-U[i+1]**2)
        # Along span
            V_s[:, i + 1], W_s[:, i + 1], betam_s[:, i + 1], alpham_s[:, i + 1], M_s[:, i + 1], Mrel_s[:,
                                                                                                       i + 1], Wt_s[:, i + 1], Vt_s[:, i + 1] = free_vortex(i + 1, rm, Vt, r_span, Vm, T, g, Ur)
            Vr_s[:, i + 1] = Vr[i + 1]
            Vz_s[:, i + 1] = Vz[i + 1]
            Vm_s[:, i + 1] = Vm[i + 1]
            betaz_s = np.degrees(np.arctan(Wt_s / Vz_s))
            alphaz_s = np.degrees(np.arctan(Vt_s / Vz_s))
            betar_s = np.degrees(np.arctan(Wt_s / Vr_s))
            alphar_s = np.degrees(np.arctan(Vt_s / Vr_s))

            beta_in, beta_out = angle_def(rownum, betaz_s, betar_s, betam_s)

        # Create t-blade3 and calculate properties
            f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf,
                               beta_in, beta_out, Mrel_s, x_s, r_s, span, xsl, rsl)

            flog = open('output_' + row_name + str(stagenum) + '.txt', 'r')

        # Flow properties calculation
            chord[:, rownum] = chord_lookup(flog, row_name, stagenum, keyword)
            pitch[rownum] = 2 * np.pi * rm[i + 1] / Z[rownum]
            sol[rownum] = chord[mean_num, rownum] / pitch[rownum]
            phi[rownum] = FlowCoeff(Vm[i + 1], U[i + 1], bw[i], rm[i], rm[i+1], axial)
            DH[rownum] = dehaller(W[i + 1], W[i])
            Re = Recalc(chord[mean_num, rownum], rho[i], W[i])
            Df[rownum] = DiffusionFact(Cp, sw[rownum], DH[rownum], Wt[i + 1], Wt[i], W[i + 1], sol[rownum],
                                       W[i + 1], rm[i + 1], W_s[-1, i], U[i], T0[i + 1], T0[i],  r_s[i, -1], axial)

        # losses calculation
            dH_Loss[0, rownum] = IncLoss(Vm[i], Wt[i], betam[i])
            dH_Loss[1, rownum], Cf[rownum] = SkinFricLoss(rho[i], rho[i + 1], W[i + 1], W_s[0, i], W_s[-1, i], xm[i], xm[i + 1], r_s[i, 0],
                                                          r_s[i, 1], rm[i], rm[i + 1], betam_s[0, i], betam_s[-1, i], betam[i], Z[rownum], chord[mean_num, rownum], bw[i + 1], phi[rownum])
            dH_Loss[4, rownum] = RecirculationLoss(
                betam[i + 1], Df[rownum], U[i + 1])
            dH_Loss[2, rownum] = BladeLoadLoss(Df[rownum], U[i + 1])
            dH_Loss[3, rownum] = ClearanceLoss(
                r_s[i, 0], r_s[i, 1], rm[i + 1],  rho[i], rho[i + 1], cl[rownum], bw[i + 1], Z[rownum], Vt[i + 1], Vm[i], U[i + 1])
            dH_Loss[5, rownum] = np.fabs(LeakageLoss(bw[i + 1], bw[i], rm[i + 1], rm[i], Vt[i + 1],
                                                     Vt[i], Z[rownum], cl[rownum], rho[i + 1], U[i + 1], chord[mean_num, rownum]))
            dH_Loss[6, rownum] = DiskFricLoss(
                U[i + 1], rm[i + 1], rho[i + 1], rho[i])
            dH[0, rownum] = Cp * (T0[i + 1] - T0[i])
            TR[rownum] = T0[i + 1] / T0[i]
            # Entalpy loss due to internal losses for each blade row
            dH[1, rownum] = np.sum(dH_Loss[0:4, rownum], axis=0)
            # Entalpy loss due to external losses for each blade row
            dH[2, rownum] = np.sum(dH_Loss[4:7, rownum], axis=0)

        # Efficiency iteration
            Etap = (dH[0, rownum] - dH[1, rownum]) / \
                (dH[0, rownum] + dH[2, rownum])
            error = np.fabs((1 - Etap / Eta[rownum]) * 100)
            print("iter= " + str(cntr) + "\t" + "Efficiency error = " + "%.8f"  % float(error))
            #print("iter = " + str(cntr) + "\t" + "Efficiency error = %.8f" + %float(error))
            Eta[rownum] = Etap
            cntr += 1

    #Thickness and stagger retrieval from t-blade3
        print()
        a0 = (g * Rgas * T0[i])**0.5
        rho0 = P0[i] / (Rgas * T0[i])
        a0rel = (g * Rgas * T0rel[i])**0.5
        rho0rel = P0rel[i] / (Rgas * T0rel[i])
        thk_max = thk_lookup(stagenum, mean_num) * chord[mean_num, rownum]
        stagger = stagger_def(stagenum, mean_num, Bckswp)

    # Check for choking
        choked, Ath, Ast = choke(
            mdot, rho0, a0, g, U[i], r_s[i, 1], r_s[i, 0], beta_in[mean_num, i], 0, thk_max, Z[rownum])
        # if choked:
        #     print("Flow is choked.....")
        #     print("Analysis terminated!!!")
        #     print("Throat area = %0.6f" %float(Ath))
        #     print("Choke area = %0.6f" %float(Ast))
        #     break
        # if choked:
        #     print("Rotor " + str(rownum) + " is choked")
        #     print("Throat area = " + "%.8f" % float(Ath) +
        #           " Flow area = " + "%.8f" % float(Ast) + '\n')
        #     print("Removing choke flow........" + '\n')
        #     inc_max = 40
        #     inc = np.linspace(0, inc_max, 1001)
        #     l = 0
        #     while choked:
        #         choked, Ath, Ast = choke(
        #             mdot, rho0, a0, g, U[i], r_s[i, 1], r_s[i, 0], beta_in[mean_num, i], inc[l], thk_max, Z[rownum])

        #         l += 1
        #     print("Successfully removed choking by incidence angle adjustment")
        #     print("Incidence angle  = " + str(-inc[l]))
        #     choked, Ath, Ast = choke(
        #         mdot, rho0, a0, g, U[i], r_s[i, 1], r_s[i, 0], betam[i], 0, thk_max, Z[rownum])
        # else:
        #     print("Rotor " + str(rownum) + " is not choked")
        #     print("Throat area = %0.8f" %float(Ath) +
        #           " Flow area = %0.8f" %float(Ast) + '\n')

    # Slip Factor
        dbetadm = (betam[i + 1] - betam[i]) / \
            (chord[mean_num, rownum] * np.cos(np.radians(stagger)))
        slip_model, shape_factor, slip_rad, slip_turn, slip_pass = SlipFactor(
            betam[i+1], gamma[i + 1], Z[rownum], pitch[rownum], phi[rownum], dbetadm, rho[i + 1], bw[i + 1], thk_max / 1000)
        T0_ac = T0[i] + (U[i + 1] * (Vt[i + 1] - U[i + 1] *
                                         (1 - slip_model)) - U[i] * Vt[i]) / Cp
        V_slip = U[i+1]*(1-slip_model)
        Vt_ac =  Vt[i+1] - V_slip
        Wt_ac = Wt[i+1] + V_slip
        error_t0 = (T0_ac-T0[i+1])*100/T0[i+1]
        print(error_t0)
        rowindex.append(row_name + str(stagenum))
        rownum += 1

#----------------Check for R-S interface and calculated at inlets---------------
    if i % 4 == 1 or i % 4 == 3:
        T0[i + 1] = T0[i]
        P0[i + 1] = P0[i]
        # M[i + 1] = Mach(area[i + 1], T0[0], P0[0], g)
        Vt[i + 1] = rm[i] * Vt[i] / rm[i + 1]  # Angular momentum conservation
        rho[i + 1] = rho[i]
        Vm[i + 1] = Vm[i]
        V[i + 1] = V[i]
        T[i + 1] = T[i]
        P[i + 1] = P[i]
        Wt[i + 1] = Vt[i + 1] - U[i + 1]
        alpham[i + 1] = np.degrees(np.arctan(Vt[i + 1] / Vm[i + 1]))
        a[i + 1], Wm[i + 1], W[i + 1], M[i+1], Mrel[i + 1], T0rel[i + 1], P0rel[i + 1], betam[i + 1], alpham[i + 1], Vz[i + 1], Vr[i + 1], betaz[i +
                                                                                                                                         1], alphaz[i + 1] = Properties_add(g, T[i + 1], Vm[i + 1], Wt[i + 1], V[i + 1], P[i + 1], rho[i + 1], Vt[i + 1], gamma[i + 1], Vz[i + 1])
        sw[i + 1] = rm[i + 1] * Vt[i + 1]
        h[i+1] = Cp*T[i+1]
        h0rel[i+1] = h[i+1] + 0.5*(W[i+1]**2)
        I[i+1] = h[i+1] + 0.5*(W[i+1]**2-U[i+1]**2)

        # Along span
        V_s[:, i + 1], W_s[:, i + 1], betam_s[:, i + 1], alpham_s[:, i + 1], M_s[:, i + 1], Mrel_s[:,
                                                                                                   i + 1], Wt_s[:, i + 1], Vt_s[:, i + 1] = free_vortex(i + 1, rm, Vt, r_span, Vm, T, g, Ur)
        Vr_s[:, i + 1] = Vr[i + 1]
        Vz_s[:, i + 1] = Vz[i + 1]
        Vm_s[:, i + 1] = Vm[i + 1]
        betaz_s = np.degrees(np.arctan(Wt_s / Vz_s))
        alphaz_s = np.degrees(np.arctan(Vt_s / Vz_s))
        betar_s = np.degrees(np.arctan(Wt_s / Vr_s))
        alphar_s = np.degrees(np.arctan(Vt_s / Vr_s))

#-------------------Check for stator exit/rotor inlet---------------------------
    if i % 4 == 2:
        # Check axial or not
        if rm[i] / rm[i + 1] >= 0.7:
            axial = True
        if rm[i] / rm[i + 1] <= 0.7:
            axial = False
        
        T0[i + 1] = T0[i]
        alpham[i + 1] = alpham[i] - dalpha[i // 4]
        P0[i + 1] = P0[i] - Y * (P0[i] - P[i])
        # M[i + 1] = Mach(area[i + 1], T0[0], P0[0], g)
        rho0 = P0[i + 1] / (Rgas * T0[i + 1])
        rho[i + 1], Vm[i + 1], V[i + 1], T[i + 1], P[i + 1], Vt[i +
                                                                1] = ThermoPropStator(rho0, M[i + 1, 0], area[i + 1], alpham[i + 1], T0[i + 1], P0[i + 1], g)
        Wt[i + 1] = Vt[i + 1] - U[i + 1]
        # s
        a[i + 1], Wm[i + 1], W[i + 1], M[i+1], Mrel[i + 1], T0rel[i + 1], P0rel[i + 1], betam[i + 1], alpham[i + 1], Vz[i + 1], Vr[i + 1], betaz[i +
                                                                                                                                         1], alphaz[i + 1] = Properties_add(g, T[i + 1], Vm[i + 1], Wt[i + 1], V[i + 1], P[i + 1], rho[i + 1], Vt[i + 1], gamma[i + 1], Vz[i + 1])
        sw[i + 1] = rm[i + 1] * Vt[i + 1]
        h[i+1] = Cp*T[i+1]
        h0rel[i+1] = h[i+1] + 0.5*(W[i+1]**2)
        I[i+1] = h[i+1] + 0.5*(W[i+1]**2-U[i+1]**2)

    # Along span
        V_s[:, i + 1], W_s[:, i + 1], betam_s[:, i + 1], alpham_s[:, i + 1], M_s[:, i + 1], Mrel_s[:,
                                                                                                   i + 1], Wt_s[:, i + 1], Vt_s[:, i + 1] = free_vortex(i + 1, rm, Vt, r_span, Vm, T, g, Ur)
        Vr_s[:, i + 1] = Vr[i + 1]
        Vz_s[:, i + 1] = Vz[i + 1]
        Vm_s[:, i + 1] = Vm[i + 1]
        betaz_s = np.degrees(np.arctan(Wt_s / Vz_s))
        alphaz_s = np.degrees(np.arctan(Vt_s / Vz_s))
        betar_s = np.degrees(np.arctan(Wt_s / Vr_s))
        alphar_s = np.degrees(np.arctan(Vt_s / Vr_s))


        row_name = "S"

        alpha_in, alpha_out = angle_def(rownum, alphaz_s, alphar_s, alpham_s)
    
    # Create t-blade3 and calculate properties
        f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf,
                           alpha_in, alpha_out, M_s, x_s, r_s, span, xsl, rsl)

        flog = open('output_' + row_name + str(stagenum) + '.txt', 'r')

    # Flow properties calculation
        chord[:, rownum] = chord_lookup(flog, row_name, stagenum, keyword)
        pitch[rownum] = 2 * np.pi * rm[i + 1] / Z[rownum]
        sol[rownum] = chord[mean_num, rownum] / pitch[rownum]
        phi[rownum] = FlowCoeff(Vm[i + 1], U[i + 1], bw[i], rm[i], rm[i+1], axial)
        DH[rownum] = dehaller(V[i + 1], V[i])
        Re = Recalc(chord[mean_num, rownum], rho[i], V[i])
        Df[rownum] = DiffusionFact(Cp, sw[rownum], DH[rownum],  Vt[i + 1], Vt[i], V[i],
                                   sol[rownum], V[i + 1], rm[i], V_s[-1, i], U[i], T0[i + 1], T0[i],  r_s[i, -1], axial)
        Rx[rownum] = DegofReac(P[i], P[i - 2], P[i - 3])

    # losses calculation
        # dH_Loss[0, rownum] = IncLoss(Vm[i], Wt[i], betam[i])
        # dH_Loss[1, rownum], Cf[rownum] = SkinFricLoss(rho[i], rho[i + 1], V[i + 1], V_s[0, i], V_s[-1, i], xm[i], xm[i + 1], r_s[i, 0],
        #                                               r_s[i, 1], rm[i], rm[i + 1], alpham_s[0, i], alpham_s[-1, i], alpham[i], Z[rownum], chord[mean_num, rownum], bw[i + 1], phi[rownum])
        # dH_Loss[4, rownum] = RecirculationLoss(
        #     alpham[i + 1], Df[rownum], U[i + 1])
        # dH_Loss[2, rownum] = BladeLoadLoss(Df[rownum], U[i + 1])
        # dH_Loss[3, rownum] = ClearanceLoss(r_s[i, 0], r_s[i, 1], rm[i + 1],  rho[i],
        #                                    rho[i + 1], cl[rownum], bw[i + 1], Z[rownum], Vt[i + 1], Vm[i], U[i + 1])
        # dH_Loss[5, rownum] = np.fabs(LeakageLoss(bw[i + 1], bw[i], rm[i + 1], rm[i], Vt[i + 1],
        #                                          Vt[i], Z[rownum], cl[rownum], rho[i + 1], U[i + 1], chord[mean_num, rownum]))
        # dH_Loss[6, rownum] = DiskFricLoss(
        #     U[i + 1], rm[i + 1], rho[i + 1], rho[i])
        rowindex.append(row_name + str(stagenum))
        rownum += 1
        stagenum += 1

    subprocess.call("geomturbo " + f.name + " 241 ")

rowindex.append("Overall")
# subprocess.call("combine_geomturbo.sh")

#------------------------------------------------------------------------------
#==============================================================================
# Calculating Design Parameters through 1-D Analysis
#==============================================================================
TR[-1] = T0[-1] / T0[0]
dH_Loss[:, -1] = np.sum(dH_Loss[:, 0:3], axis=1)  # Overall enthaly loss
dH[1, :] = np.sum(dH_Loss[0:4, :], axis=0)        # Entalpy loss due to internal losses for each blade row
dH[2, :] = np.sum(dH_Loss[4:7, :], axis=0)      # Entalpy loss due to external losses for each blade row
Eta = (dH[0, :] - dH[1, :]) / (dH[0, :] + dH[2, :])  # Efficiency

for i in range(nrows):
    if i % 2 == 1:
        Eta[i] = 0

PR = (1 + Eta * (TR - 1))**(g / (g - 1))  # Pressure Ratio
Eta_poly = ((g - 1) / g) * np.log(PR) / np.log(TR)

#Flow coefficient of compressor
if r_span[-1,-1]/r_span[-1,0]>=1 and r_span[-1,-1]/r_span[-1,0]<=1.1:
    axial = True
else:
    axial = False

phi_ov = FlowCoeff(Vm[0], U[0], bw[0], rm[0], rm[-1], axial)

# First cell distance using Blasius solution
ywall_s = 6 * ((Vm_s / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8))
print("Estimate of first cell wall distance =", np.amin(ywall_s))

#------------------------------------------------------------------------------
#==============================================================================
# Exporting design point data
#==============================================================================

fmean = open("meanline_data_despt.dat", 'w')

#----------------------Stage Quantities-------------------------
stage_qty = ["Row", "Solidity", "DF", "Cf",
             "DeHallerNumber", "Rx", "phi", "PR", "Efficiency"]

row_list = [[] for i in range(nrows + 1)]

for i in range(nrows + 1):
    if i == 0:
        for j in range(len(stage_qty)):
            row_list[i].append(stage_qty[j])
    else:
        row_list[i].append("%02d" % (i))
        row_list[i].append("%.4f" % float(sol[i - 1]))
        row_list[i].append("%.4f" % float(Df[i - 1]))
        row_list[i].append("%.4f" % float(Cf[i - 1]))
        row_list[i].append("%.4f" % float(DH[i - 1]))
        row_list[i].append("%.4f" % float(Rx[i - 1]))
        row_list[i].append("%.4f" % float(phi[i - 1]))
        row_list[i].append("%.4f" % float(PR[i - 1]))
        row_list[i].append("%.4f" % float(Eta[i - 1]))

fmean.write('\n' + tabulate(row_list, headers="firstrow") + '\n')

#----------------------Station Quantities-------------------------
station_qty = ["J", "Swirl", "Vt[m/s]", "Vm[m/s]", "Vz[m/s]",  " Vr[m/s]",
               " V[m/s]",  " W[m/s]",  " T[k]",  "Mach", "Rel.Mach", "P0[Pa]",  "T0[k]"]
station_list = [[] for i in range(nstations + 1)]

for i in range(nstations + 1):
    if i == 0:
        for j in range(len(station_qty)):
            station_list[i].append(station_qty[j])
    else:
        station_list[i].append("%02d" % (i))
        station_list[i].append("%.4f" % float(sw[i - 1]))
        station_list[i].append("%.4f" % float(Vt[i - 1]))
        station_list[i].append("%.4f" % float(Vm[i - 1]))
        station_list[i].append("%.4f" % float(Vz[i - 1]))
        station_list[i].append("%.4f" % float(Vr[i - 1]))
        station_list[i].append("%.4f" % float(V[i - 1]))
        station_list[i].append("%.4f" % float(W[i - 1]))
        station_list[i].append("%.4f" % float(T[i - 1]))
        station_list[i].append("%.4f" % float(M[i - 1, 0]))
        station_list[i].append("%.4f" % float(Mrel[i - 1]))
        station_list[i].append("%.4f" % float(P0[i - 1]))
        station_list[i].append("%.4f" % float(T0[i - 1]))

fmean.write('\n' + tabulate(station_list, headers="firstrow") + '\n')

#------------------------------------------------------------------------------
loss_qty_rows = ["Inc. Loss", "Skin Friction Loss", "Blade Loading Loss",
                 "Clearance Loss", "Recirculation Loss", "Leakage Loss", "Disk Friction Loss"]
loss = pd.DataFrame(data=dH_Loss.T, index=rowindex,
                    columns=loss_qty_rows).astype('float')
fmean.write('\n' + tabulate(loss, headers=loss_qty_rows, numalign="left") + '\n')

#----------------------Overall Quantities-------------------------
fmean.write('\n\n')
fmean.write("Overall Pressure ratio = %2.4f" % PR[-1] + '\n')
fmean.write("Overall Isentropic Efficiency = %2.4f" % (Eta[-1]*100) + '\n')
fmean.write("Overall Polytropic Efficiency = %2.4f" % (Eta_poly[-1]*100) + '\n')
fmean.write("Exit rotor Slip Factor = %.4f" % slip_model + "\n")
fmean.write("Exit rotor Shape Factor = %.4f" % shape_factor + "\n")
fmean.write("FLow Coefficient = %.4f" % phi_ov + "\n")
fmean.write('\n')
fmean.close()

#==============================================================================

plots(xsl, rsl, Vm_s, Vt_s, W_s, Wt_s, alpham_s, betam_s, span, nstns, bsf, rowindex)
stop = timeit.default_timer()
print("Execution Time: ", '%1.3f' % (stop - start), "seconds")
