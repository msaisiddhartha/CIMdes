import numpy as np
import subprocess, os, errno
from  inputs import *
from design import *


def ThermoPropRotor(rho0, area, Vt, T0, P0, g):
    density = rho0
    error = 5
    while abs(error) > 1e-6:
        Vmerid = mdot / (density * area)
        Vabs = (Vmerid**2 + Vt**2)**0.5
        Ts = T0 - Vabs**2 / (2 * Cp)
        Ps = P0 * ((Ts / T0)**((g / (g - 1))))
        rho2p = Ps / (Rgas * Ts)
        error = (1 - rho2p / density) * 100
        density = rho2p
        #print(density)
    return density, Vmerid, Vabs, Ts, Ps

def ThermoPropStator(rho0, area, alpha_out, T0, P0, g):
    density = rho0
    error = 5
    while abs(error) > 1e-6:
        Vmerid = mdot / (density * area)
        Vtang = np.tan(np.deg2rad(alpha_out)) * Vmerid
        Vabs = (Vmerid**2 + Vtang**2)**0.5
        Ts = T0 - Vabs**2 / (2 * Cp)
        Ps = P0 * ((Ts / T0)**((g / (g - 1))))
        rho2p = Ps / (Rgas * Ts)
        error = (1 - rho2p / density) * 100
        density = rho2p

    return density, Vmerid, Vabs, Ts, Ps, Vtang

def Properties_add(g, T1, Vm1, Wt1, V1, W1, P1, rho1, Vt1, phi_angle1,Vz1):
    a1      = (g * Rgas * T1)**0.5
    Wm1     = Vm1
    W1      = (Wt1**2 + Wm1**2)**0.5
    M1      = V1 / a1
    Mrel1   = W1 / a1
    T0rel1  = T1 + (W1**2 / (2 * Cp))
    P0rel1  = P1 + (0.5 * rho1 * W1**2)
    betam1  = np.degrees(np.arctan(Wt1 / Vm1))
    alpham1 = np.degrees(np.arctan(Vt1 / Vm1))
    Vz1     = Vm1*np.cos(np.radians(phi_angle1))
    Vr1     = Vm1*np.sin(np.radians(phi_angle1))
    betaz1  = np.degrees(np.arctan(Wt1 / Vz1))
    alphaz1 = np.degrees(np.arctan(Vt1 / Vz1))
    return a1, Wm1, W1, M1, Mrel1, T0rel1 ,P0rel1 ,betam1 ,alpham1,Vz1,Vr1, betaz1, alphaz1

def free_vortex(k, rm1, Vt1, r_span1, Vm1, T1, g, Ur1):
    Vm_s1[:,k] = Vm1[k]
    T_s1[:,k] = T1[k]
    for j in range(nsect):
        Vt_s1[j, k] = rm1[k] * Vt1[k] / r_span1[j, k]
    sw_s1[:,k]      = r_span1[:,k]*Vt_s1[:,k]
    Wt_s1[:,k]      = Vt_s1[:,k] - Ur1[:,k]
    betam_s1[:,k]   = np.degrees(np.arctan(Wt_s1[:,k] / Vm_s1[:,k]))
    alpham_s1[:,k]  = np.degrees(np.arctan(Vt_s1[:,k] / Vm_s1[:,k]))
    V_s1[:,k]       = (Vt_s1[:,k]**2 + Vm_s1[:,k]**2)**0.5
    W_s1[:,k]       = (Wt_s1[:,k]**2 + Vm_s1[:,k]**2)**0.5
    a_s1[:,k]       = (g * T_s1[:,k] * Rgas)**0.5
    M_s1[:,k]       = V_s1[:,k] / a_s1[:,k]
    Mrel_s1[:,k]    = W_s1[:,k] / a_s1[:,k]

    return V_s1[:,k], W_s1[:,k], betam_s1[:,k], alpham_s1[:,k], M_s1[:,k], Mrel_s1[:,k], Wt_s1[:,k], Vt_s1[:,k]

def span_var(k, rm1, Vt1, r_span1, Vm1, T1, g, Ur1):
    #exponential law
    Vm_s1[:,k] = Vm1[k]
    T_s1[:,k] = T1[k]
    a=100
    b=40
    if k%2==0:
        for j in rage(nsect):
            Vt_s1[j,k] = a-(b*r_span1[j,k]/rm1[k])
    else:
        for j in rage(nsect):
            Vt_s1[j,k] = a+(b*r_span1[j,k]/rm1[k])

    sw_s1[:,k]      = r_span1[:,k]*Vt_s1[:,k]
    Wt_s1[:,k]      = Vt_s1[:,k] - Ur1[:,k]
    betam_s1[:,k]   = np.degrees(np.arctan(Wt_s1[:,k] / Vm_s1[:,k]))
    alpham_s1[:,k]  = np.degrees(np.arctan(Vt_s1[:,k] / Vm_s1[:,k]))
    V_s1[:,k]       = (Vt_s1[:,k]**2 + Vm_s1[:,k]**2)**0.5
    W_s1[:,k]       = (Wt_s1[:,k]**2 + Vm_s1[:,k]**2)**0.5
    a_s1[:,k]       = (g * T_s1[:,k] * Rgas)**0.5
    M_s1[:,k]       = V_s1[:,k] / a_s1[:,k]
    Mrel_s1[:,k]    = W_s1[:,k] / a_s1[:,k]

def DegofReac(P2, P1, P0):
    Rx = (P1 - P0) / (P2 - P0)
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

def stagger_def(rn, mn, beta0):
    flog = open("splinedata_section." + str(mn) + ".R" + str(rn) + ".dat", 'r')
    nskip = 4
    npt = 3
    for i in range(nskip):
        datastr = flog.readline()
    line = np.zeros((npt,1))
    for i in range(npt):
        datastr = flog.readline()
        datam = datastr.split()
        line[i] = datam[3]
    camber_ang = np.degrees(np.arctan(line[2]))
    stagger_ang = beta0 - camber_ang
    flog.close()
    return stagger_ang

def thk_lookup(rn, mn):
    flog = open("blade_Section_data.R" + str(rn) + ".dat", 'r')
    nskip = 4+ mean_num
    for i in range(nskip):
        datastr = flog.readline()
    datastr = flog.readline()
    datam = datastr.split()
    thk = np.float(datam[11])
    return thk

def angle_def(rn, betaz_s, betar_s, betam_s):
    if ang[rn] == 0:
        beta_in_blade = betaz_s
        beta_out_blade = betaz_s
    elif ang[rn] == 1:
        beta_in_blade = betar_s
        beta_out_blade = betar_s
    elif ang[rn] == 2:
        beta_in_blade = betaz_s
        beta_out_blade = betar_s
    return beta_in_blade, beta_out_blade

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
def create_tblade3(k, cntr, row_name, stagenum, data, nsect, bsf, beta_in, beta_out, Mrel_s, x_s, r_s, span, xsl, rsl):
    x_cp_le = np.linspace(x_s[cntr,0], x_s[cntr,1], nsect)
    r_cp_le = np.linspace(r_s[cntr,0], r_s[cntr,1], nsect)
    x_cp_te = np.linspace(x_s[cntr+1,0], x_s[cntr+1,1], nsect)
    r_cp_te = np.linspace(r_s[cntr+1,0], r_s[cntr+1,1], nsect)
    f = open(str(casename) + "." + str(k + 1) + ".dat", "w")
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
    f.write("    " + str(ang[k]) + " inci_dev_spline" + '\n')
    f.write(" Airfoil camber defined by curvature control (0=no,1=yes):" + '\n')
    f.write("    " + str(1) + '\t'+ 'spanwise_spline' +  '\n')
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
        f.write(('%2.8f' % beta_in[i, cntr]) + "  ")
        f.write(('%2.8f' % beta_out[i, cntr + 1]) + "  ")
        f.write(('%2.8f' % Mrel_s[i, cntr]) + "  ")
        f.write(('%2.8f' % chrdr_nd) + "  ")
        f.write(('%2.8f' % 0.0150) + "  ")
        f.write(('%2.8f' % 0.0000) + "  ")
        f.write(('%2.8f' % 0.0000) + "  ")
        f.write(('%2.8f' % 0.0000) + "\n")
    f.write('\n')
    f.write(" LE / TE curve (x,r) definition :" + '\n')
    f.write(" Number of Curve points :" + '\n')
    f.write("    " + str(nsect) + '\n')
    f.write("   xLE          rLE           xTE          rTE" + '\n')
    for i in range(nsect):
        f.write("    " + ('%2.8f' % (x_cp_le[i] / bsf)) + "  ")
        f.write(('%2.8f' % (r_cp_le[i] / bsf)) + "  ")
        f.write(('%2.8f' % (x_cp_te[i] / bsf)) + "  ")
        f.write(('%2.8f' % (r_cp_te[i] / bsf)) + "\n")
    f.write('\n')
    f.write(" # Airfoil type and Variable Radial Stacking information.         #" + '\n')
    f.write(" # stack_u: % chord stack (0.00 to 100.00).                       #" + '\n')
    f.write(" # stack_v: % below or above meanline stack (-100.00 to +100.00). #" + '\n')
    f.write(" # Use +200 for stacking on airfoil area centroid.                #" + '\n')
    f.write(" Variable Radial stacking (0=no,1=yes):" + '\n')
    f.write("    " + ('%01d' % 0) + '\n')
    f.write(
        " J   type |stk_u |stk_v |umxthk |lethk |tethk  |Jcells(Grid:4n+1) |eta_ofst(<=10){thkc/Jmax}  |BGgrid(0=no,1=yes) |" + '\n')
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
    #subprocess.call("geomturbo " + f.name + " 241 ")
    fout.close()
    return f

def make_sure_path_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
