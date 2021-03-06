import numpy as np
from scipy.optimize import fsolve
from inputs import *

npts = 50  # Number of pts on streamline
np_inlt = 5
np_inlt_eq = 2
stp = 0.1
np_outlt = 25
bw = np.zeros(nstations)
gamma = np.zeros(nstations)
#------------------------------------------------------------------------------
#==============================================================================
# Generating streamlines at hub nad tip using the krain single-stage compressor
#==============================================================================


def streamlines():
# Hub
    r_hub = np.linspace(R_hub_le, R_hub_te, npts)
    x_hub = X01 + (R**2 - (r_hub - R01)**2)**0.5
# Casing
    r_tip = np.linspace(R_tip_le, R_tip_te, npts)
    x_tip = X04 + (ae / be) * (be**2 - (r_tip - R04)**2)**0.5

# Non-dimesionalizing with leading edge hub radius as scaling factor
    bsf = r_tip[0]
    x_hub_nd = x_hub / bsf
    r_hub_nd = r_hub / bsf
    x_tip_nd = x_tip / bsf
    r_tip_nd = r_tip / bsf

#------------------------------------------------------------------------------
#==============================================================================
# Splitting the blade into 3 parts for rotor and stator
#==============================================================================

    # x = np.array((npts,nstations))
    # r = np.array((npts,nstations))

    r_s[0, :] = np.array([r_hub[0], r_tip[0]])
    x_s[0, :] = np.array([x_hub[0], x_tip[0]])
    for i in range(1, nstations):
        if i == nstations - 1:
            r_s[-1, :] = np.array([r_hub[-1], r_tip[-1]])
            x_s[-1, :] = np.array([x_hub[-1], x_tip[-1]])
        else:
            r_s[i, :] = r_s[i - 1, :] + gap[i - 1]
            x_s[i, :] = np.array([X01 + (R**2 - (r_s[i, 0] - R01)**2)**0.5, X04 + (ae / be) * (be**2 -
                                                                                               (r_s[i, 1] - R04)**2)**0.5])

#-----------------------------Stator-------------------------------------------
#            r_s[2, :] = r_s[1, :]+gap
#            x_s[2, :] = np.array([X01 + (R**2-(r_s[2, 0]-R01)**2)**0.5,X04 + (ae/be)*(be**2-
#                      (r_s[2, 1]-R04)**2)**0.5])
#
#            r_s[3, :] = r_s[2, :]+s1_len
#            x_s[3, :] = np.array([X01 + (R**2-(r_s[3, 0]-R01)**2)**0.5,X04 + (ae/be)*(be**2-
#                      (r_s[3, 1]-R04)**2)**0.5])
#
#--------------------------------Rotor 2---------------------------------------
#            r_s[4, :] = r_s[3, :]+gap
#            x_s[4, :] = np.array([X01 + (R**2-(r_s[4, 0]-R01)**2)**0.5,X04 + (ae/be)*(be**2-
#                      (r_s[4, 1]-R04)**2)**0.5])


#----------------Area calculation at inlet and oultlet-------------------------
    for i in range(nstations):
        # Station 1(Inlet)
        rm[i] = np.mean(r_s[i, :])
        xm[i] = np.mean(x_s[i, :])
        area[i] = np.pi * (sum(r_s[i, :])) * ((r_s[i, 0] -
                                               r_s[i, 1])**2 + (x_s[i, 0] - x_s[i, 1])**2)**0.5
    bw = ((x_s[:, 0] - x_s[:, 1])**2 + (r_s[:, 0] - r_s[:, 1])**2)**0.5
    phi = np.fabs(np.degrees(
        np.arctan((r_s[:, 0] - r_s[:, 1]) / (x_s[:, 0] - x_s[:, 1]))))
    gamma = 90 - phi

#-----------------------------Streamlines--------------------------------------

    # Inlet Streamlines
    x_inlt = np.zeros(np_inlt) + (-stp * np_inlt)
    r_hub_inlt = np.zeros(np_inlt)
    r_tip_inlt = np.ones(np_inlt) * r_tip_nd[0]
    for i in range(0, np_inlt - 1):
        x_inlt[i + 1] = x_inlt[i] + stp

    for j in range(np_inlt - 1, -1, -1):
        if j >= np_inlt - np_inlt_eq:
            r_hub_inlt[j] = R01 / bsf - \
                ((R / bsf)**2 - (x_inlt[j] - (X01 / bsf))**2)**0.5
        else:
            r_hub_inlt[j] = r_hub_inlt[j + 1]

    # Outlet Streamlines
    R_DiffuserExit_nd = R_DiffuserExit / bsf
    r_outlt = np.linspace(r_hub_nd[-1], R_DiffuserExit_nd, np_outlt + 1)
    r_outlt = r_outlt[1:]
    x_hub_outlt = np.zeros(np_outlt) + x_hub_nd[-1]
    x_tip_outlt = np.zeros(np_outlt)

    # Diffuser
    Aoutlt = 2 * np.pi * (x_hub_nd[-1] - x_tip_nd[-1]) * r_hub_nd[-1]
    Aexit = 0.99 * (Aoutlt)

    area_diff = np.linspace(Aoutlt, Aexit, np_outlt)
    x_tip_outlt = x_hub_outlt - area_diff / (2 * np.pi * r_outlt)

    x_hub_nd = np.hstack((x_inlt, x_hub_nd, x_hub_outlt))
    r_hub_nd = np.hstack((r_hub_inlt, r_hub_nd, r_outlt))
    x_tip_nd = np.hstack((x_inlt, x_tip_nd, x_tip_outlt))
    r_tip_nd = np.hstack((r_tip_inlt, r_tip_nd, r_outlt))

    xsl = np.zeros((nsect, len(x_hub_nd)))
    rsl = np.zeros((nsect, len(x_hub_nd)))
    xsl_mean = np.zeros((1,len(x_hub_nd)))
    rsl_mean = np.zeros((1,len(x_hub_nd)))
    for j in range(len(x_hub_nd)):
        xsl[:, j] = np.linspace(x_hub_nd[j], x_tip_nd[j], nsect)
        rsl[:, j] = np.linspace(r_hub_nd[j], r_tip_nd[j], nsect)

#Calculating for LE and TE of blade rows using rm values
    xsl_mean = xsl[mean_num, :]*bsf
    rsl_mean =  rsl[mean_num, :]*bsf
    for i in range(0, nstations):
        #Searching for lower and upper r-values
        upb = rsl_mean>rm[1]
        lpb = rsl_mean<rm[1]
        upper_bound_r = rsl_mean[upb].min()
        lower_bound_r = rsl_mean[lpb].max()
        #Finding index of arry for lower and upper r-values
        mini = np.where(rsl_mean == lower_bound_r)
        maxi = np.where(rsl_mean == upper_bound_r)
        #Finding lower and upper x-values corresponding to r-values through index search
        upper_bound_x = xsl_mean[maxi[0][0]]
        lower_bound_x = xsl_mean[mini[0][0]]
        points = [(lower_bound_x, lower_bound_r), (upper_bound_x, upper_bound_r)]
        x_coords, r_coords = zip(*points)
        A = np.vstack([x_coords, np.ones(len(x_coords))]).T
        m, c = np.linalg.lstsq(A, r_coords)[0]
        def eq_hub(p):
            x, r = p
            return (r-m*x-c, x-X01 + (R**2 - (r - R01)**2)**0.5)
        x, r = fsolve(eq_hub, (lower_bound_x, lower_bound_r))
    return xm, rm, area, x_s, r_s, bsf, xsl, rsl, bw, gamma, m,c

