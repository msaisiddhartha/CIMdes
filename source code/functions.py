import numpy as np
import subprocess, os
from  inputs import *
from design import *

def DegofReac(P2, P1, P0):
    Rx = P1 - P0 / (P2 - P0)
    return Rx

def Recalc(L, rho, V):
    Re = rho*V*L/nu
    return Re

def dehaller(W1, W0):
    DH = W1 / W0
    return DH

def FlowCoeff(Vm1, U0):
    phi = Vm1 / U0
    return phi

def chord_lookup(flog, row_name, stagenum, keyword):
    ch = np.zeros(nsect)
    current_sec = 0
    lineread = flog.readlines()
    for num, line in enumerate(lineread, 1):
        linelist = line.split()
        if keyword in linelist:
            ch[current_sec] = float(linelist[1])
            current_sec += 1
    flog.close()
    return ch / 1000


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

def DiffusionFact(Cp, sw, DH, Wt1, Wt0, W0, sol, W1, rmean, W_tip0, U1, T01, T00, rtip, axial):
    if axial:
        Df = 1 - DH + (Wt1 - Wt0) / (2 * sol * W0)
    if not axial:
        Wratio = W1 / W_tip0
        num = 0.75 * Cp * (T01 - T00) / U1**2
        rad_ratio = rtip / rmean
        Df = 1 - Wratio + (Wratio * num) / \
            (sol * (1 - rad_ratio) + 2 * rad_ratio)
    return Df


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

# Flow angle switch
def create_tblade3(k, cntr, row_name, stagenum, data, nsect, bsf, beta_s, alpha_s, Mrel_s, M_s, x_s, r_s, span, xsl, rsl):
    f = open("tblade3input." + str(k + 1) + "." + str(casename) + ".dat", "w")
    f.write("Input parameters (version 1.1)" + '\n')
    f.write("    " + row_name + str(stagenum) + '\n')
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
        f.write("    " + ('%2.8f' % (x_s[cntr, i] / bsf)) + "  ")
        f.write(('%2.8f' % (r_s[cntr, i] / bsf)) + "  ")
        f.write(('%2.8f' % (x_s[cntr + 1, i] / bsf)) + "  ")
        f.write(('%2.8f' % (r_s[cntr + 1, i] / bsf)) + "\n")
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
    f.close()

#Executing T-blade3 and generate .geomturbo file
    fout = open('output_' + row_name + str(stagenum) + '.txt', 'w')
    subprocess.Popen(["tblade3", f.name], stdout=fout,
                     stderr=fout).communicate()[0]
    subprocess.call("geomturbo " + f.name + " 241 ")
    fout.close()
