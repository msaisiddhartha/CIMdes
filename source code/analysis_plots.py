import matplotlib.pyplot as py
import numpy as np
import subprocess
from inputs import *
#------------------------------------------------------------------------------
#==============================================================================
#Plotting meridional view, velocity triangle and blade angles
#==============================================================================

#-----------------------------Meridional view=---------------------------------
def plots(xsl, rsl, Vm, Vt, W, Wt, alpha, beta, span, nstns, bsf):
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
    py.savefig("streamlines.png")

    fignum+=1

#---------------------------Velocity triangles---------------------------------
    fig = py.figure(figsize=(15, 10))
    for i in range(nstns):
        ax = fig.add_subplot(2,2,i+1)
        soa = np.array([[0,0,Vm[0][i],0],[0,0,Vm[0][i],Vt[0][i]],[Vm[0][i],0,0,
                         Vt[0][i]],[0,0,Vm[0][i],Wt[0][i]],[Vm[0][i],0,0,Wt[0][i]]])
        X, Y, U, V = zip(*soa)
        ax = py.gca()
        ax.set_xlim([-10, Vm[0][i]+10])
        if Wt[0][i]<0:
            ax.set_ylim([Vt[0][i]+10, Wt[0][i]-10])
            ax.text(25,-30,r'$\beta$ = %.2f$\degree$'%(beta[0][i]),size='16')
            ax.text(25,10,r'$\alpha$ = %.2f$\degree$'%(alpha[0][i]),size='16')
        else:
            ax.set_ylim([Vt[0][i]+10, -10])
            ax.text(25,10,r'$\beta$ = %.2f$\degree$'%(beta[0][i]),size='16')
            ax.text(25,75,r'$\alpha$ = %.2f$\degree$'%(alpha[0][i]),size='16')
        ax.invert_yaxis()
        ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
        ax.set_title("Station "+str(i+1))
    py.tight_layout()
    py.savefig("veltri.png")

    fignum+=1

#----------------------------Blade Angles--------------------------------------
    py.figure(fignum, figsize=(16,9))
    py.plot(beta[:,0],span,'k',marker='.',label=r'$\beta_{in}$')
    py.plot(beta[:,1],span,'r',marker='.',label=r'$\beta_{out}$')
    py.legend()
    py.ylabel('span',size='24')
    py.xlabel(r'$\beta_z$',size='24')
    leg = py.gca().get_legend()
    ltext = leg.get_texts()
    py.setp(ltext,fontsize='16')
    py.tight_layout()
    py.savefig("rotor1.png")

    fignum+=1

    py.figure(fignum, figsize=(16,9))
    py.plot(alpha[:,2],span,'k',marker='.',label=r'$\alpha_{in}$')
    py.plot(alpha[:,3],span,'r',marker='.',label=r'$\alpha_{out}$')
    py.ylabel('span',size='24')
    py.xlabel(r'$\alpha_z$',size='24')
    py.legend()
    leg = py.gca().get_legend()
    ltext = leg.get_texts()
    py.setp(ltext,fontsize='16')
    py.tight_layout()
    py.savefig("stator1.png")

    fignum+=1
    py.figure(fignum, figsize=(16,9))
    py.plot(beta[:,4],span,'k',marker='.',label=r'$\beta_{in}$')
    py.plot(beta[:,5],span,'r',marker='.',label=r'$\beta_{out}$')
    py.ylabel('span',size='24')
    py.xlabel(r'$\beta_r$',size='24')
    py.legend()
    leg = py.gca().get_legend()
    ltext = leg.get_texts()
    py.setp(ltext,fontsize='16')
    py.tight_layout()
    py.savefig("rotor2.png")
#py.show()
