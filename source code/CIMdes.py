"""Multi-rotor compressor. All units are in SI system."""
import matplotlib.pyplot as py
import numpy as np
import subprocess, timeit, os
from design import *
from inputs import *
from functions import *


from analysis_plots import *


start = timeit.default_timer()

workdir = os.getcwd()

#------------------------------------------------------------------------------
#==============================================================================
# Calculating Design Parameters through 1-D Analysis
#==============================================================================

"At all Sections along span"
# First cell distance using Blasius solution
ywall_s[:, 0] = 6 * ((Vm[0] / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
ywall_s[:, 1] = 6 * ((Vm[1] / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
ywall_s[:, 2] = 6 * ((Vm[4] / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
ywall_s[:, 3] = 6 * ((Vm[5] / 0.0000157)**(-7 / 8)) * ((0.2)**(1 / 8));
print("Estimate of first cell wall distance =",
      min(min(ywall_s[0]), min(ywall_s[1])))

# Solidity
#for i in range(nrows):
    #pitch_s[:, i] = (2 * np.pi * r_blde[:, i + 1] * 1000) / Z[i]
#    sol_s[:, i] = chord[:, i] / pitch_s[:, i]

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

stage_num = 1
fout = []
data = []

c = 0
for i in range(nrows):
    row_name, stagenum = rowname(i+1, stage_num)
    chord[:,i] = create_tblade3(i, c, row_name, stagenum, fout, data)
    c += 2
c = 0
for i in range(nrows):
    pitch[i] = (2 * np.pi * rm[c + 1] * 1000) / Z[i]
    sol[i] = chord[mean_num, i] / pitch[i]
    Rx[i] = (np.tan(np.deg2rad(beta[c + 1])) + np.tan(np.deg2rad(beta[c]))) * Vm[c] / (2 * U[c])
    phi[i] = Vm[c + 1] / U[c]
    DH[i] = W[c+1] / W[c]
    Df[i] = (1 - W[c+1] / W[c]) + (sw[c+ 1] - sw[c]) / ((rm[c] + rm[c + 1]) * sol[i] * W[c])
    c+=2

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
print(" Execution Time: ", '%1.3f'%(stop - start) , "seconds")
