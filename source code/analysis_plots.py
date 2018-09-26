import pylab as py
from matplotlib.pyplot import *
import numpy as np
import subprocess
import sys, os
from inputs import *
from functions import *
from adjustText import adjust_text

workdir = os.getcwd()
if plot:
    if os.name == 'nt':
        plotdir = workdir + '\plots'
    elif os.name == 'posix':
        plotdir = workdir + '/plots'
    make_sure_path_exists(plotdir)
    ## Check if the plot directory exist and if not then make one
#------------------------------------------------------------------------------
#==============================================================================
#Plotting meridional view, velocity triangle and blade angles
#==============================================================================

#-----------------------------Meridional view=---------------------------------
def plots(xsl, rsl, Vm, Vt, W, Wt, alpha, beta, span, nstns, bsf, r_id):
    print()
    print()
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print()
    print("Generating Plots...")
    print()
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    fignum = 1
    sz=15
    py.figure(fignum, figsize=(16,9))
    py.plot(xsl[0],rsl[0],color='g',marker='.',lw=1)
    py.plot(xsl[-1],rsl[-1],color='g',marker='.',lw=1)
    for i in range(0,nstations):
        lines=py.plot(x_s[i,:]/bsf,r_s[i,:]/bsf)
        py.setp(lines,color = 'k',lw=0.75)
    for i in range(1,nstations-1):
        lines=py.plot(xsl[i],rsl[i])
        py.setp(lines,color = 'blue',linestyle='--',lw=0.75)
    py.xlabel('x-coordinate')
    py.ylabel('r-coordinate')
    py.text(0.25,0.41,'R1',size=sz)
    py.text(0.75,0.81,'S1',size=sz)
    py.text(1,1.3,'R2',size=sz)
    py.text(-0.05,1.01,'1',size=sz)
    py.text(0.3,1.075,'2',size=sz)
    py.text(0.37,1.1,'3',size=sz)
    py.text(0.61,1.25,'4',size=sz)
    py.text(0.67,1.3,'5',size=sz)
    py.text(0.87,1.778,'6',size=sz)
    py.text(0.92,2.65,'7',size=sz)
    py.axis('equal')
    py.grid(True)
    py.tight_layout()
    py.savefig(os.path.join(plotdir, "streamlines.png"))

    fignum+=1

#---------------------------Velocity triangles---------------------------------

    vel_tri = {0:"hub", mean_num:"mean", -1:"tip"}
    for j in vel_tri.keys():
        fig = py.figure(fignum, figsize=(15, 10))
        texts= []
        for i in range(nstns):
            ax = fig.add_subplot(nrows//2+1,2,i+1)
            soa = np.array([[0,0,Vm[j][i],0],[0,0,Vm[j][i],Vt[j][i]],[Vm[j][i],0,0,
                             Vt[j][i]],[0,0,Vm[j][i],Wt[j][i]],[Vm[j][i],0,0,Wt[j][i]]])
            X, Y, U, V = zip(*soa)
            ax = py.gca()
            ax.set_xlim([-10, max(Vm[:,i])+10])
            ax.set_ylim([max(Vt[:,i])+10, min(Wt[:,i])-10])
            if Wt[j][i]<0:
                ax.text(25,-30,r'$\beta$ = %.2f$\degree$'%(beta[j][i]), size=16)
                ax.text(25,10,r'$\alpha$ = %.2f$\degree$'%(alpha[j][i]), size=16)
            else:
                ax.text(25,10,r'$\beta$ = %.2f$\degree$'%(beta[j][i]), size=16)
                ax.text(25,75,r'$\alpha$ = %.2f$\degree$'%(alpha[j][i]), size=16)
            ax.invert_yaxis()
            ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
            ax.set_title("Station "+str(i+1))
        py.tight_layout()
        py.savefig(os.path.join(plotdir, "veltri_" + vel_tri[j] + ".png"))
        fignum+=1

#----------------------------Blade Angles--------------------------------------
    count = 0
    for i in range(0, nrows):
        py.figure(fignum, figsize=(16,9))
        if i%2==0:
            py.plot(beta[:,count],span,'k',marker='.',label=r'$\beta_{in}$')
            py.plot(beta[:,count+1],span,'r',marker='.',label=r'$\beta_{out}$')
        else:
            py.plot(alpha[:,count],span,'k',marker='.',label=r'$\alpha_{in}$')
            py.plot(alpha[:,count+1],span,'r',marker='.',label=r'$\alpha_{out}$')
        count+=2
        py.legend()
        py.ylabel('Normalized span',size='24')
        py.xlabel(r'$\beta_z$',size='24')
        leg = py.gca().get_legend()
        ltext = leg.get_texts()
        py.setp(ltext,fontsize='16')
        py.tight_layout()
        py.savefig(os.path.join(plotdir, r_id[i] + ".png"))
        fignum+=1

#py.show()
