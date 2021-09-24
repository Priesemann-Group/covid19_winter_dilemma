import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl


# reset defaults
mpl.rcParams.update(mpl.rcParamsDefault)

mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"

#mpl.rcParams['axes.linewidth'] = 0.3
mpl.rcParams["axes.labelcolor"] = "black"
mpl.rcParams["axes.edgecolor"] = "black"
mpl.rcParams["xtick.color"] = "black"
mpl.rcParams["ytick.color"] = "black"

mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False

mpl.rcParams["xtick.labelsize"] = 8
mpl.rcParams["ytick.labelsize"] = 8
mpl.rcParams["axes.labelsize"] = 10
mpl.rcParams["axes.titlesize"]= 10
mpl.rcParams["legend.fontsize"] = 8
mpl.rcParams["legend.title_fontsize"] = 10

colors_paul = ["#233954", "#ea5e48", "#1e7d72", "#f49546", "#e8bf58", # dark
               "#5886be", "#f3a093", "#53d8c9", "#f2da9c", "#f9c192"] # light


def overview(model, path=None, silent=False):
    t = model.times
    data = model.chopped_data()

    fig = plt.figure(figsize=(7, 4), constrained_layout=True)
    grid = fig.add_gridspec(ncols=3, nrows=2, hspace=0.2, wspace=0.15)
    
    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1],sharex=ax1)
    ax3 = fig.add_subplot(grid[2],sharex=ax1)
    ax4 = fig.add_subplot(grid[3],sharex=ax1)
    ax5 = fig.add_subplot(grid[4],sharex=ax1)
    ax6 = fig.add_subplot(grid[5],sharex=ax1)
    
    
    ax1.plot(t,data[:,0]/1e6,label="S")
    ax1.plot(t,data[:,1]/1e6,label="V")
    ax1.plot(t,data[:,2]/1e6,label="W")
    ax1.plot(t,data[:,8]/1e6,label="R")
    ax1.set_ylim(0,None)
    ax1.set_ylabel("fraction of population")
    ax1.legend(loc='upper left', ncol=2)
    
    ax2.plot(t,data[:,3],label="E")
    ax2.plot(t,data[:,4],label="EB")
    ax2.plot(t,data[:,5],label="I")
    ax2.plot(t,data[:,6],label="IB")
    ax2.set_ylim(0,None)
    ax2.legend(loc='upper left', ncol=2)
    
    ax3.plot(t,data[:,7],label="ICU")
    ax3.plot(t,list(map(model.H_Rt, t)), label="H_Rt")
    ax3.plot(t,list(map(model.H_vac1, t)), label="H_vac1")
    ax3.plot(t,list(map(model.H_vac2, t)), label="H_vac2")
    ax3.set_ylim(0,None)
    ax3.legend(loc='upper left', ncol=2)
    
    ax4.plot(t, list(map(model.Rt, t)), label='Rt')
    ax4.plot(t, list(map(model.Gamma, t)), label='season.')
    ax4.plot(t, list(map(model.R_0, t)), label='R_0')
    ax4.set_ylabel("reproduction number")
    ax4.legend(loc='upper left', ncol=2)
    
    ax5.plot(t, list(map(model.u_w, t)), label='u_w')
    ax5.plot(t, 1/model.M*data[:,9], label='u_c')
    ax5.plot(t, list(map(model.w_w, t)), label='w_w')
    ax5.plot(t, 1/model.M*data[:,10], label='w_c')
    ax5.set_ylim(0,None)
    ax5.set_xlabel("days")
    ax5.set_ylabel("fraction of population")
    ax5.legend(loc='upper left', ncol=2)
    
    ax6.plot(t, np.array(list(map(model.Phi, t, data[:,9])))*model.M, label='1.dose')    # avoid the list*int thingy
    ax6.plot(t, np.array(list(map(model.phi, t, data[:,10])))*model.M, label='2.dose')
    ax6.set_ylim(0,None)
    ax6.set_ylabel("daily vaccinations")
    ax6.legend(loc='upper left')
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
        l,u = ax.get_ylim()
        ax.set_ylim(l,u+0.4*(u-l))
    
    fig.align_ylabels()
    
    if not silent: plt.show()
    if path!=None: fig.savefig(path)
