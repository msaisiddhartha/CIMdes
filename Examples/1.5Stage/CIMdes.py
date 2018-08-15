"""Multi-rotor compressor. All units are in SI system."""
import pylab as py
from matplotlib.pyplot import *
import numpy as np
import subprocess
import timeit
import os

from inputs import *
from design import *
from loss import *
from functions import *


from analysis_plots import *
np.seterr(divide='ignore', invalid='ignore')

start = timeit.default_timer()

workdir = os.getcwd()

rm, area, r_s, bsf, xsl, rsl, bw, gamma, phi_angle = streamlines()




# Evaluating based on known parameters
w = 2 * np.pi * N / 60  # Angular velocity[rad/s]
g = 1 / (1 - Rgas / Cp)  # Ratio of Specific heats[-]
U = w*rm
Vt[0] = Vt1
T0[0] = T01
P0[0] = P01

delTT_row = WorkRatio * delTT
T0[-1] = T0[0] + delTT
dH[0,-1] = Cp*(T0[-1] - T0[0])

#Flow properties along span
r_span = np.zeros((nsect, nstations))
span = np.linspace(0, 1, nsect)
for i in range(nstations):
    for j in range(nsect):
        r_span[j, i] = span[j] * (r_s[i, 1] - r_s[i, 0]) + r_s[i, 0]

Ur = w * r_span

#Rotor-1 inlet calculations
Wt[0] = Vt[0] - U[0]
rho0 = P0[0] / (Rgas * T0[0])
rho[0], Vm[0], V[0], T[0], P[0] = ThermoPropRotor(rho0, area[0], Vt[0], T0[0], P0[0], g)
a[0], Wm[0], W[0], M[0], Mrel[0], T0rel[0] ,P0rel[0] ,betam[0] ,alpham[0],Vz[0],Vr[0], betaz[0], alphaz[0] = Properties_add(g, T[0], Vm[0], Wt[0], V[0], W[0], P[0], rho[0], Vt[0], gamma[0], Vz[0])
sw[0] = rm[0] * Vt[0]
V_s[:,0], W_s[:,0], betam_s[:,0], alpham_s[:,0], M_s[:,0], Mrel_s[:,0], Wt_s[:,0], Vt_s[:,0] = free_vortex(0, rm, Vt, r_span, Vm, T, g, Ur)
Vr_s[:,0] = Vr[0]
Vz_s[:,0] = Vz[0]
Vm_s[:,0] = Vm[0]

stagenum = 1
rownum=0
Eta=np.ones((nrows+1,1))

keyword = "chord_actual(mm):"
#Calcualtions from rotor-1 exit onwards
for i in range(0, nstations-1):

    error = 5
    #Detect rotor exit
    if i%4==0:
        #Along meanline
        T0[i+1] = T0[i] + delTT_row[i//4]
        print("Chekcing for R" + str(stagenum) + " efficiency convergence")
        cntr=0
        while np.fabs(error)>1e-6:
            P0[i+1] = P0[i] * (1+Eta[rownum]*((T0[i+1]/T0[i])-1))**(g / (g - 1))
            Vt[i+1] = (Cp * (T0[i+1] - T0[i]) + U[i] * Vt[i]) / U[i+1]
            Wt[i+1] = Vt[i+1] - U[i+1]
            rho0 = P0[i+1] / (Rgas * T0[i+1])
            rho[i+1], Vm[i+1], V[i+1], T[i+1], P[i+1] = ThermoPropRotor(rho0, area[i+1], Vt[i+1], T0[i+1], P0[i+1], g)
            a[i+1], Wm[i+1], W[i+1], M[i+1], Mrel[i+1], T0rel[i+1] ,P0rel[i+1] ,betam[i+1] ,alpham[i+1],Vz[i+1],Vr[i+1], betaz[i+1], alphaz[i+1] = Properties_add(g, T[i+1], Vm[i+1], Wt[i+1], V[i+1], W[i+1], P[i+1], rho[i+1], Vt[i+1], gamma[i+1], Vz[i+1])
            sw[i+1] = rm[i+1] * Vt[i+1]

            #Along span
            V_s[:,i+1], W_s[:,i+1], betam_s[:,i+1], alpham_s[:,i+1], M_s[:,i+1], Mrel_s[:,i+1], Wt_s[:,i+1], Vt_s[:,i+1] = free_vortex(i+1, rm, Vt, r_span, Vm, T, g, Ur)
            Vr_s[:,i+1] = Vr[i+1]
            Vz_s[:,i+1] = Vz[i+1]
            Vm_s[:,i+1] = Vm[i+1]
            betaz_s = np.degrees(np.arctan(Wt_s / Vz_s))
            alphaz_s = np.degrees(np.arctan(Vt_s / Vz_s))
            betar_s = np.degrees(np.arctan(Wt_s / Vr_s))
            alphar_s = np.degrees(np.arctan(Vt_s / Vr_s))
            #Check axial or not
            if rm[i]/rm[i+1] >=0.7:
                axial = True
            if rm[i]/rm[i+1] <=0.7:
                axial = False

            row_name = "R"

            #Create t-blade3 and calculate properties
            if ang[rownum] == 0:

                f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf, betaz_s, betaz_s, Mrel_s, x_s, r_s, span ,xsl, rsl)
            elif ang[rownum] == 1:

                f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf, betar_s, betar_s, Mrel_s, x_s, r_s, span ,xsl, rsl)
            elif ang[rownum] == 2:

                f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf, betaz_s, betar_s, Mrel_s, x_s, r_s, span ,xsl, rsl)
            flog = open('output_' + row_name + str(stagenum) + '.txt', 'r')


            chord[:,rownum] = chord_lookup(flog, row_name, stagenum, keyword)
            pitch[rownum] = 2 * np.pi * rm[i + 1] / Z[rownum]
            sol[rownum] = chord[mean_num, rownum] / pitch[rownum]
            phi[rownum] = FlowCoeff(Vm[i+1], U[i+1])

            #losses calculation
            DH[rownum] = dehaller(W[i+1], W[i])
            Re = Recalc(chord[mean_num, rownum], rho[i], W[i])
            Df[rownum] = DiffusionFact(Cp, sw[rownum], DH[rownum], Wt[i+1], Wt[i], W[i], sol[rownum], W[i+1], rm[i], W_s[-1, i], U[i], T0[i+1], T0[i],  r_s[i, -1], axial)
            dH_Loss[0, rownum] = IncLoss(Vm[i], Wt[i], betam[i])
            dH_Loss[1,rownum] = SkinFricLoss(W[i+1], W[i], r_s[i,0], r_s[i,1], betam[i], Z[rownum], Re, chord[mean_num, rownum])
            dH_Loss[4,rownum] = RecirculationLoss(betam[i+1], Df[rownum], U[i+1])
            dH_Loss[2,rownum] = BladeLoadLoss(Df[rownum], U[i+1])
            dH_Loss[3,rownum] = ClearanceLoss(r_s[i,0], r_s[i,1], rm[i+1], cl[rownum], bw[i+1], U[i+1])
            dH_Loss[5,rownum] = np.fabs(LeakageLoss(bw[i+1], bw[i], rm[i+1], rm[i], Vt[i+1], Vt[i], Z[rownum], cl[rownum], rho[i+1], U[i+1], chord[mean_num, rownum]))
            dH_Loss[6,rownum] = DiskFricLoss(U[i+1], rm[i+1], rho[i+1], rho[i], Re)
            dH[0,rownum] = Cp*(T0[i+1] - T0[i])
            TR[rownum] = T0[i+1]/T0[i]

            dH[1,rownum] = np.sum(dH_Loss[0:4,rownum], axis =0)   #Entalpy loss due to internal losses for each blade row
            dH[2,rownum] = np.sum(dH_Loss[4:7,rownum], axis =0)   #Entalpy loss due to external losses for each blade row
            Etap = (dH[0,rownum]-dH[1,rownum])/(dH[0,rownum]+dH[2,rownum])
            error = np.fabs((1-Etap/Eta[rownum])*100)
            print("iter = "+ str(cntr) + "\t" + "Efficiency error = " + str(error))
            Eta[rownum] = Etap
            cntr+=1
        print()
        #Slip Factor
        if stagenum == 2:
            stagger  = stagger_def(stagenum, mean_num, Bckswp)
            thk_max  = thk_lookup(stagenum, mean_num) * chord[mean_num, rownum]
            dbetadm = (Bckswp-betam[i])/(chord[mean_num, rownum] * np.cos(np.radians(stagger)))
            slip_model, shape_factor, slip_rad, slip_turn, slip_pass = SlipFactor(Bckswp, gamma[i+1], Z[rownum], pitch[rownum], phi[rownum], dbetadm, rho[i+1], bw[i+1], thk_max)
            Vslip = U[i+1]*(1+phi[rownum]*np.tan(np.radians(Bckswp))) - Vt[i+1]

            slip_calc = 1-Vslip/U[i+1]
        rownum+=1

    #Check for R-S interface and keep properties constant along interface
    if i%4 == 1 or i%4 ==3:
        T0[i+1] = T0[i]
        P0[i+1] = P0[i]
        Vt[i+1] = rm[i] * Vt[i]/rm[i+1]         #Angular momentum conservation
        rho[i+1] = rho[i]
        Vm[i+1] = Vm[i]
        V[i+1] = V[i]
        T[i+1] = T[i]
        P[i+1]=P[i]
        Wt[i+1] = Vt[i+1] - U[i+1]
        alpham[i+1] = np.degrees(np.arctan(Vt[i+1] / Vm[i+1]))
        a[i+1], Wm[i+1], W[i+1], M[i+1], Mrel[i+1], T0rel[i+1] ,P0rel[i+1] ,betam[i+1] ,alpham[i+1],Vz[i+1],Vr[i+1], betaz[i+1], alphaz[i+1] = Properties_add(g, T[i+1], Vm[i+1], Wt[i+1], V[i+1], W[i+1], P[i+1], rho[i+1], Vt[i+1], gamma[i+1], Vz[i+1])
        sw[i+1] = rm[i+1] * Vt[i+1]

        #Along span
        V_s[:,i+1], W_s[:,i+1], betam_s[:,i+1], alpham_s[:,i+1], M_s[:,i+1], Mrel_s[:,i+1], Wt_s[:,i+1], Vt_s[:,i+1] = free_vortex(i+1, rm, Vt, r_span, Vm, T, g, Ur)
        Vr_s[:,i+1] = Vr[i+1]
        Vz_s[:,i+1] = Vz[i+1]
        Vm_s[:,i+1] = Vm[i+1]
        betaz_s = np.degrees(np.arctan(Wt_s / Vz_s))
        alphaz_s = np.degrees(np.arctan(Vt_s / Vz_s))
        betar_s = np.degrees(np.arctan(Wt_s / Vr_s))
        alphar_s = np.degrees(np.arctan(Vt_s / Vr_s))

    #Check for stator exit and calculate properties
    if i%4 == 2:
        T0[i+1] = T0[i]
        P0[i+1] = P0[i] - Y * (P0[i] - P0[i])
        alpham[i+1] = alpham[i] - dalpha[i//4]
        rho0 = P0[i+1] / (Rgas * T0[i+1])
        rho[i+1], Vm[i+1], V[i+1], T[i+1], P[i+1], Vt[i+1] = ThermoPropStator(rho0, area[i+1], alpham[i+1], T0[i+1], P0[i+1], g)
        Wt[i+1] = Vt[i+1] - U[i+1]
        #s
        a[i+1], Wm[i+1], W[i+1], M[i+1], Mrel[i+1], T0rel[i+1] ,P0rel[i+1] ,betam[i+1] ,alpham[i+1],Vz[i+1],Vr[i+1], betaz[i+1], alphaz[i+1] = Properties_add(g, T[i+1], Vm[i+1], Wt[i+1], V[i+1], W[i+1], P[i+1], rho[i+1], Vt[i+1], gamma[i+1], Vz[i+1])
        sw[i+1] = rm[i+1] * Vt[i+1]

        #Along span
        V_s[:,i+1], W_s[:,i+1], betam_s[:,i+1], alpham_s[:,i+1], M_s[:,i+1], Mrel_s[:,i+1], Wt_s[:,i+1], Vt_s[:,i+1] = free_vortex(i+1, rm, Vt, r_span, Vm, T, g, Ur)
        Vr_s[:,i+1] = Vr[i+1]
        Vz_s[:,i+1] = Vz[i+1]
        Vm_s[:,i+1] = Vm[i+1]
        betaz_s = np.degrees(np.arctan(Wt_s / Vz_s))
        alphaz_s = np.degrees(np.arctan(Vt_s / Vz_s))
        betar_s = np.degrees(np.arctan(Wt_s / Vr_s))
        alphar_s = np.degrees(np.arctan(Vt_s / Vr_s))
        #Check axial or not
        if rm[i]/rm[i+1] >=0.7:
            axial = True
        if rm[i]/rm[i+1] <=0.7:
            axial = False

        row_name = "S"

        #Create t-blade3 and calculate properties
        if ang[rownum] == 0:

            f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf, alphaz_s, alphaz_s, M_s, x_s, r_s, span ,xsl, rsl)
        if ang[rownum] == 1:

            f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf, alphar_s, alphar_s, M_s, x_s, r_s, span ,xsl, rsl)
        if ang[rownum] == 2:

            f = create_tblade3(rownum, i, row_name, stagenum, data, nsect, bsf, alphaz_s, alphar_s, M_s, x_s, r_s, span ,xsl, rsl)
        flog = open('output_' + row_name + str(stagenum) + '.txt', 'r')
        chord[:,rownum] = chord_lookup(flog, row_name, stagenum, keyword)
        pitch[rownum] = 2 * np.pi * rm[i + 1] / Z[rownum]
        sol[rownum] = chord[mean_num, rownum] / pitch[rownum]
        phi[rownum] = FlowCoeff(Vm[i+1], U[i+1])

        #losses calculation
        DH[rownum] = dehaller(V[i+1], V[i])
        Re = Recalc(chord[mean_num, rownum], rho[i], V[i])
        Df[rownum] = DiffusionFact(Cp, sw[rownum], DH[rownum], Vt[i+1], Vt[i], V[i], sol[rownum], V[i+1], rm[i], V_s[-1, i], U[i], T0[i+1], T0[i],  r_s[i, -1], axial)
        dH_Loss[0, rownum] = IncLoss(Vm[i], Wt[i], betam[i])
        dH_Loss[1,rownum] = SkinFricLoss(V[i+1], V[i], r_s[i,0], r_s[i,1], alpham[i], Z[rownum], Re, chord[mean_num, rownum])
        dH_Loss[4,rownum] = RecirculationLoss(alpham[i+1], Df[rownum], U[i+1])
        dH_Loss[2,rownum] = BladeLoadLoss(Df[rownum], U[i+1])
        dH_Loss[3,rownum] = ClearanceLoss(r_s[i,0], r_s[i,1], rm[i+1], cl[rownum], bw[i+1], U[i+1])
        dH_Loss[5,rownum] = np.fabs(LeakageLoss(bw[i+1], bw[i], rm[i+1], rm[i], Vt[i+1], Vt[i], Z[rownum], cl[rownum], rho[i+1], U[i+1], chord[mean_num, rownum]))
        dH_Loss[6,rownum] = DiskFricLoss(U[i+1], rm[i+1], rho[i+1], rho[i], Re)

        rownum+=1
        stagenum+=1

    subprocess.call("geomturbo " + f.name + " 241 ")


TR[-1] = T0[-1]/T0[0]
dH_Loss[:,3] = np.sum(dH_Loss[:,0:3], axis =1)  #Overall enthaly loss
dH[1,:] = np.sum(dH_Loss[0:4,:], axis =0)   #Entalpy loss due to internal losses for each blade row
dH[2,:] = np.sum(dH_Loss[4:7,:], axis =0)   #Entalpy loss due to external losses for each blade row

Eta = (dH[0,:]-dH[1,:])/(dH[0,:]+dH[2,:])       #Efficiency

PR = (1+Eta*(TR-1))**(g/(g-1))              #Pressure Ratio
#------------------------------------------------------------------------------
#==============================================================================
# Calculating Design Parameters through 1-D Analysis
#==============================================================================

#At all Sections along span
ywall_s = 6 * ((Vm_s / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8))   # First cell distance using Blasius solution
print("Estimate of first cell wall distance =", np.amin(ywall_s))

fmean = open("meanlineflowproperties.dat", 'w')
fmean.write("Row    Solidity    DF    DeHallerNumber      Rx      phi      PR       Efficiency\n")
for i in range(nrows):
    fmean.write("%02d" % (i + 1) + "    ")
    fmean.write("  " + '%2.4f' %
                sol[i] + '  ' + '%2.4f' % Df[i])
    fmean.write("    " + "%2.4f" % DH[i])
    fmean.write("    " + "\t%2.4f" % -Rx[i])
    fmean.write("    " + "%2.4f" % phi[i])
    fmean.write("    " + "%2.4f" % PR[i])
    fmean.write("    " + "%2.4f" % Eta[i])
    fmean.write('\n')
fmean.write('\n')
fmean.write("Station    Swirl   Vt[m/s]    Vm[m/s]    T[k]      Mach\n")
for i in range(nstns):
    fmean.write("%02d" % (i + 1) + "  ")
    fmean.write("  " + '%2.4f' % sw_s[mean_num, i] + '    ' + '%2.4f' % Vt_s[mean_num, i] + '    ' + '%2.4f' %
                Vm_s[mean_num, i] + '    ' + '%2.4f' % T_s[mean_num, i] + '    ' + '%2.4f' % M_s[mean_num, i] + '\n')
fmean.write('\n\n')
fmean.write("Overall Pressure ratio = %2.4f" % PR[-1] + '\n')
fmean.write("Overall Efficiency = %2.4f" % Eta[-1] + '\n')
fmean.write('\n')
fmean.close()


plots(xsl, rsl, Vm_s, Vt_s, W_s, Wt_s, alpham_s, betam_s, span, nstns, bsf)

stop = timeit.default_timer()
print(" Execution Time: ", '%1.3f' % (stop - start), "seconds")
