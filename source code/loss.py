import numpy as np
from inputs import *
from functions import *
import math

# Suffix nomenclature: 0-inlet/LE; 1-exit/TE


def SlipFactor(beta_b1, gamma1, Z, s, phi1, dBetadm1, rho1, blade_width1, thk):
    F = 1 + 2*np.sin(np.pi / Z) * np.sin((np.pi / Z) +
                                       np.radians(beta_b1)) * np.cos(np.radians(beta_b1)) * np.sin(np.radians(gamma1)) - thk/(s*np.cos(np.radians(beta_b1))) #Shape FActor
    dSlip_rad = F * np.pi * \
        np.cos(np.radians(beta_b1)) * np.sin(np.radians(gamma1)) / Z
    dSlip_turn = F * s * phi1 * dBetadm1 *0.1 / (4 * np.cos(np.radians(beta_b1)))
    #dSlip_pass = F * phi1 * s * np.sin(np.radians(beta_b1)) / (4 * rho1 * blade_width1)
    dSlip_pass = 0
    slip = 1 - dSlip_rad - dSlip_turn - dSlip_pass

    return slip, F, dSlip_rad, dSlip_turn, dSlip_pass

def IncLoss(Vm0, Wt0, betab0):
    f_inc = 0.5
    dH_inc = f_inc * (Vm0 * np.tan(np.deg2rad(betab0)) - Wt0)**2 / 2
    return dH_inc


def BladeLoadLoss(Df, U):
    dH_bld = 0.05 * (Df * U)**2
    return dH_bld


def SkinFricLoss(W1, W0, rh, rt, betaf, Z, Re, Lb):
    W_avg = (W1**2 + W0**2) ** 0.5
    cf = 0.027 / Re**(1 / 7)
    nu2 = rh / rt
    d_hyd = 0.2 / 4
    p = Z / (np.pi * np.cos(np.radians(betaf)))
    dH_sf = 2 * cf * Lb * (W_avg**2) / d_hyd
    return dH_sf


def DiskFricLoss(U1, rmean1, rho1, rho0, Re):
    dH_df = 0.01356 * rho1 * U1**3 * (2 * rmean1)**2 / (mdot * Re)
    return dH_df


def ClearanceLoss(rhub, rtip, rmean1, clearance, blade_width1, U1):
    dH_cl = 2 * clearance * ((rhub + rtip) / (2 * rmean1) - 0.275) * U1**2
    return dH_cl


def RecirculationLoss(alpha1, Df, U1):
    dH_rc = abs(8e-5 * math.sinh(3.5 * np.radians(alpha1)**3) * ((Df * U1)**2))
    return dH_rc


def LeakageLoss(blade_width0, blade_width1, rmean1, rmean0, Vt0, Vt1, Z, clearance, rho1, U1, L_theta):
    ravg = 0.5 * (rmean1 + rmean0)
    bavg = 0.5 * (blade_width1 + blade_width0)
    dH_lk = rho1 * clearance * U1 * 1.332 * \
        (rmean1 * Vt1 - rmean0 * Vt0) / (2 * ravg * bavg)
    return dH_lk
