import numpy as np
from inputs import *
#from functions import *
import math

# Suffix nomenclature: 0-inlet/LE; 1-exit/TE


def SlipFactor(beta_b1, gamma1, Z, s, phi1, dBetadm1, rho1, blade_width1, thk): #Qiu
    F = 1 - 2 * np.sin(np.pi / Z) * np.sin((np.pi / Z) +
                                           np.radians(-beta_b1)) * np.cos(np.radians(-beta_b1)) * np.sin(np.radians(gamma1)) - thk / (s * np.cos(np.radians(beta_b1)))  # Shape Factor
    dSlip_rad = F * np.pi * \
        np.cos(np.radians(beta_b1)) * np.sin(np.radians(gamma1)) / Z
    dSlip_turn = F * s * phi1 * dBetadm1 * \
        0.1 / (4 * np.cos(np.radians(beta_b1)))
    #dSlip_pass = F * phi1 * s * np.sin(np.radians(beta_b1)) / (4 * rho1 * blade_width1)
    dSlip_pass = 0
    slip = 1 - dSlip_rad - dSlip_turn - dSlip_pass

    return slip, F, dSlip_rad, dSlip_turn, dSlip_pass


def choke(m, rho01, a01, g, Utip, r0tip, r0hub, beta1b, inc, thk, Zblade):
    Astar = m / (rho01 * a01 * ((2 + (g - 1) * (Utip / a01)**2) /
                                (g + 1))**((g + 1) / (2 * (g - 1))))
    Athroat = np.pi * (r0tip**2 - r0hub**2) * np.cos(np.radians(beta1b - inc)
                                                     ) - (thk * Zblade * (r0tip - r0hub)) / 1000
    if Astar < Athroat:
        return True, Athroat, Astar
    else:
        return False, Athroat, Astar


def IncLoss(Vm0, Wt0, betab0):  # Conrad
    f_inc = 0.4
    dH_inc = f_inc * (Wt0 - Vm0 * np.tan(np.radians(betab0)))**2
    return dH_inc


def BladeLoadLoss(Df, U):  # Coppage
    dH_bld = 0.05 * (Df * U)**2
    return dH_bld


def SkinFricLoss(rho0, rho1, W1, W_h0, W_t0, x0, x1, r_h0, r_t0, r1, r2, beta_h0, beta_t0, beta1, Z, Lb, bw1, floc):  # Jansen
    W_avg = (2 * W1 + W_h0 + W_t0) / 4
    r_hub_rat = r_h0 / r_t0
    d_hyd_rat = np.cos(np.radians(beta1)) / ((Z / np.pi) + (2 * r1 * np.cos(np.radians(beta1)) / bw1)) + \
        (0.5 * ((r_h0 + r_t0) / r1) * 0.5 *
         (np.cos(np.radians(beta_h0)) + np.cos(np.radians(beta_t0)))) / (Z / np.pi + (r_h0 + r_t0) * 0.5 * (np.cos(np.radians(beta_h0)) + np.cos(np.radians(beta_t0))) / (r_t0 - r_h0))
    d_hyd = d_hyd_rat * 2 * r1
    #p = Z / (np.pi * np.cos(np.radians(betaf)))
    Lz = 2*r2*(0.014 + (0.023*r2/r_h0) + 1.58 * floc)
    Lb = (np.pi / 4) * (2 * r1 - (r_h0 + r_t0) - bw1 + (2 * Lz)) * (1 /
                                                                    (np.cos(np.radians(beta_h0)) + np.cos(np.radians(beta_t0)) + np.cos(np.radians(beta1))))
    rho_avg = 0.5*(rho0+rho1)
    Re = rho_avg*d_hyd*W_avg/nu
    cf = 0.0412 * Re**(-0.1925)
    dH_sf = 2 * cf * Lb * W_avg**2 / d_hyd
    return dH_sf, cf


def DiskFricLoss(U1, rmean1, rho1, rho0):  # Daily and Nece
    Re_df = U1 * rho1 * rmean1 / nu
    if Re_df >= 3e5:
        f_df = 0.0622 / Re_df**0.2
    if Re_df < 3e5:
        f_df = 2.67 / Re_df**0.5
    rho_avg = 0.5 * (rho0 + rho1)
    dH_df = f_df * (rho_avg * rmean1**2 * U1**3) / (4 * mdot)
    #dH_df = 0.01356 * rho1 * U1**3 * (2 * rmean1)**2 / (mdot * Re)
    return dH_df


def ClearanceLoss(r_h0, r_t0, r1, rho0, rho1, clearance, bw1, Z, Vt1, Vm0, U1):  # Jansen
    num1 = 4*np.pi*np.fabs(Vt1)*Vm0/ (bw1*Z)
    num2 = np.fabs((r_h0**2 - r_t0**2) / ((r1 - r_t0) * (1 + rho1 / rho0)))
    dH_cl = 0.6 * (clearance / bw1) * np.fabs(Vt1) * (num1 * num2)**0.5
    return dH_cl


def RecirculationLoss(alpha1, Df, U1):  # Oh
    dH_rc = abs(8e-5 * math.sinh(3.5 * np.radians(alpha1)**3) * ((Df * U1)**2))
    return dH_rc


def LeakageLoss(bw0, bw1, rmean1, rmean0, Vt0, Vt1, Z, clearance, rho1, U1, L_theta):  # Aungier
    ravg = 0.5 * (rmean1 + rmean0)
    bavg = 0.5 * (bw1 + bw0)
    dH_lk = rho1 * clearance * U1 * 1.332 * \
        (rmean1 * Vt1 - rmean0 * Vt0) / (2 * ravg * bavg)
    return dH_lk
