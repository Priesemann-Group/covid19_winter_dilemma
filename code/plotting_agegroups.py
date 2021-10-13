import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl


def set_rcParams(arial=False):
    # reset defaults
    mpl.rcParams.update(mpl.rcParamsDefault)

    if arial:
        mpl.rcParams['font.family'] = "sans-serif"
        mpl.rcParams['font.sans-serif'] = "Arial"

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


def overview(model, path=None, silent=False, arial=False):
    set_rcParams(arial=arial)
    mpl.rcParams["legend.fontsize"] = 7
    t = model.times
    data = model.chopped_data().sum(axis=2)
    AGdata = model.chopped_data()

    fig = plt.figure(figsize=(6, 3.5), constrained_layout=True)
    grid = fig.add_gridspec(ncols=3, nrows=2, hspace=0.2, wspace=0.15)
    
    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1],sharex=ax1)
    ax3 = fig.add_subplot(grid[2],sharex=ax1)
    ax4 = fig.add_subplot(grid[3],sharex=ax1)
    ax5 = fig.add_subplot(grid[4],sharex=ax1)
    ax6 = fig.add_subplot(grid[5],sharex=ax1)
    
    
    ax1.plot(t,data[:,0]/1e6,label="S")
    ax1.plot(t,data[:,1]/1e6,label="V")
    ax1.plot(t,data[:,2]/1e6,label="Wn")
    ax1.plot(t,data[:,3]/1e6,label="Wv")
    ax1.plot(t,data[:,12]/1e6,label="R")
    ax1.plot(t,data[:,13]/1e6,label="Rv")
    ax1.set_ylim(0,None)
    ax1.set_ylabel("fraction of population")
    ax1.legend(loc='upper left', ncol=2, handlelength=1.)
    
    ax2.plot(t,data[:,4],label="E")
    ax2.plot(t,data[:,5],label="EBn")
    ax2.plot(t,data[:,6],label="EBv")
    ax2.set_ylim(0,None)
    ax2.legend(loc='upper left', ncol=2, handlelength=1.)
    
    ax3.plot(t,data[:,10],label="ICU")
    ax3.plot(t,data[:,11],label="ICUv")
#    ax3.plot(t,list(map(model.H_Rt, t)), label="H_Rt")
    ax3.plot(t,list(map(model.H_vac1, t)), label="H_vac1")
#    ax3.plot(t,list(map(model.H_vac2, t)), label="H_vac2")
    ax3.set_ylim(0,None)
    ax3.legend(loc='upper left', ncol=2, handlelength=1.)
    
    ax4.plot(t, list(map(model.Rt, t)), label='Rt')
    ax4.plot(t, list(map(model.Gamma, t)), label='season.')
    ax4.plot(t, list(map(model.R_0, t)), label='R_0')
    ax4.set_ylabel("reproduction number")
    ax4.legend(loc='upper left', ncol=2, handlelength=1.)
    
    ax5.plot(t, list(map(model.u_w, t)), label='u_w')
    ax5.plot(t, data[:,14]/model.M.sum(), label='u_c')
    ax5.plot(t, list(map(model.w_w, t)), label='w_w')
    ax5.plot(t, data[:,15]/data[:,14], label='w_c')
    ax5.set_ylim(0,None)
    ax5.set_xlabel("days")
    ax5.set_ylabel("fraction of population")
    ax5.legend(loc='upper left', ncol=2, handlelength=1.)
    
#    d1a = (np.array(list(map(model.Phi, t+model.tau_vac1, AGdata[:,14,:])))).sum(axis=1)
#    d1b = (np.array(list(map(model.Phi, t+model.tau_vac1/2., AGdata[:,14,:])))).sum(axis=1)
#    d2 = (np.array(list(map(model.phi, t+model.tau_vac2, AGdata[:,14,:], AGdata[:,15,:])))).sum(axis=1)
#    ax6.plot(t, d1a, label='1.dose A')
#    ax6.plot(t, d1b, label='1.dose B')
    d1a = model.X1[model.time2index(0):-100,:].sum(axis=1)
    d1b = model.X1[model.time2index(0-model.tau_vac1/2.):model.time2index(model.t_max-model.tau_vac1/2.),:].sum(axis=1)
    d2 = model.X2[model.time2index(0):-100,:].sum(axis=1)
    ax6.plot(t, d1a, label='1.dose A')
    ax6.plot(t, d1b, label='1.dose B')
    ax6.plot(t, d2, label='2.dose')
    ax6.set_ylim(0,None)
    ax6.set_ylabel("daily vaccinations")
    ax6.legend(loc='upper left', ncol=2, handlelength=1.)
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
        l,u = ax.get_ylim()
        ax.set_ylim(l,u+0.4*(u-l))
    
    fig.align_ylabels()
    
    if not silent: plt.show()
    if path!=None: fig.savefig(path)




def sixpanels(models, path=None, silent=False, arial=False, ICUcap=None, full_wave=None):
    set_rcParams(arial=arial)
    mpl.rcParams["legend.fontsize"] = 7

    m1, m2, m3 = models
    t = m1.times
    data = [m1.chopped_data().sum(axis=2), m2.chopped_data().sum(axis=2), m3.chopped_data().sum(axis=2)]
    AGdata = [m1.chopped_data(), m2.chopped_data(), m3.chopped_data()]

    fig = plt.figure(figsize=(6., 3.5), constrained_layout=True)
    grid = fig.add_gridspec(ncols=3, nrows=2, wspace=0.1)

    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1], sharex=ax1)
    ax3 = fig.add_subplot(grid[2], sharex=ax1)
    ax4 = fig.add_subplot(grid[3], sharex=ax1)
    ax5 = fig.add_subplot(grid[4])
    ax6 = fig.add_subplot(grid[5])

    colors = {
        'low':'#0099B4FF', 'mid':'#00468BFF', 'high':'#1B1919FF', 'free':'#1B1919FF',
        'lowL':'#FFFFFFFF', 'midL':'#FFFFFFFF', 'highL':'#FFFFFFFF',
#        'lowL':'#0099B499', 'midL':'#00468B99', 'highL':'#1B191999',
        'line':'#ADB6B6FF', 'ICUcap':'#FFAAAA',
        'now':'#ADB6B6FF', 'nowL':'#FFFFFFFF',
#        'now':'#93dfedFF', 'nowL':'#93dfed99',
    }
    main_colors = [colors['low'],colors['mid'],colors['high']]
    main_colors_L = [colors['lowL'],colors['midL'],colors['highL']]



#    ax1.plot(t[1800:], np.ones(1800), color=colors['free'])
    
    for i,m in enumerate([m1,m2,m3]):
        ax1.plot(t, np.array(list(map(m.R_0, t)))/m.Rt_free, color=main_colors[i])
        ax2.plot(t, m.rho*(data[i][:,4]+data[i][:,5]+data[i][:,6]), color=main_colors[i])
        ax3.plot(t, data[i][:,10]+data[i][:,11], color=main_colors[i])
        d1a = m.X1[m.time2index(0):-100,:].sum(axis=1)
        d1b = m.X1[m.time2index(0-m.tau_vac1/2.):m.time2index(m.t_max-m.tau_vac1/2.),:].sum(axis=1)
        d2 = m.X2[m.time2index(0):-100,:].sum(axis=1)
        ax4.plot(t, d1a+d1b+d2, color=main_colors[i])

    ax5.bar(1, data[0][0,1]/1e6, 0.5,
        align='center', color=colors['now'], edgecolor='black', zorder=-1)
    ax5.bar(1, (data[0][0,12]+data[0][0,13])/1e6, 0.5,
        align='center', color=colors['nowL'], edgecolor='black', zorder=-1, bottom=data[0][0,1]/1e6)

    offset = 0.5
    for i in [2,4]:
        for ab,m,j in zip([0.5,0,-0.5],[m1,m2,m3],[0,1,2]):
            ax5.bar(offset+i+ab, data[j][900*i-1,1]/1e6, 0.5,  
                align='center', color=main_colors[j], edgecolor='black', zorder=-1)
            ax5.bar(offset+i+ab, (data[j][900*i-1,12]+data[j][900*i-1,13])/1e6, 0.5,
                align='center', color=main_colors_L[j], edgecolor='black', zorder=-1,
                    bottom=data[j][900*i-1,1]/1e6)
            ax6.bar(offset+i+ab, data[j][900*i-1,16], 0.5, 
                align='center', color=main_colors[j], edgecolor='black', zorder=-3)
 


    for ax in [ax1,ax2,ax3,ax4]:
        ax.axvline(180, ls=':', color=colors['line'], zorder=0)
        ax.set_xlabel('2021            2022')


    ax1.set_ylim(0,1.1)
    ax2.set_ylim(0,None)
    ax3.set_ylim(0,None)
    ax4.set_ylim(0,None)
    ax5.set_ylim(0,1)
    ax6.set_ylim(0,None)

    ax5.set_yticks([0.,0.25,0.5,0.75,1.])

    ax1.set_ylabel("Contact levels\ninfluenced by NPIs")
    ax2.set_ylabel("Daily new cases\nper million")
    ax3.set_ylabel("ICU occupancy\nper million")
#    ax4.set_ylabel("Reproduction number")
    ax4.set_ylabel("Daily vaccinations\nper million")
    ax5.set_ylabel("Immune fraction\nof the population")
    ax6.set_ylabel("Total deaths\nper million")

    #Panel 1
    ax1.text(0,m1.Rt_base/5+0.05,'Scenario 3', size=7, color=colors['low'])
    ax1.text(0,m2.Rt_base/5+0.05,'Scenario 2', size=7, color=colors['mid'])
    ax1.text(0,m3.Rt_base/5+0.05,'Scenario 1', size=7, color=colors['high'])
    ax1.text(200,m3.Rt_base/5+0.05,'No restrictions', size=7, color=colors['line'])
    
    #Lifting of restrictions
    ax1.text(0.54,0.05,'Lifting of\nrestrictions', size=7, color=colors['line'], transform=ax1.transAxes)
    for ax in [ax2,ax3,ax4]:
        l,u = ax.get_ylim()
        ax.set_ylim(l,u+0.15*(u-l))
        ax.text(0.54,0.9,'Lifting of\nrestrictions', size=7, color=colors['line'], transform=ax.transAxes)

    for ax, label in zip([ax1,ax2,ax3,ax4,ax5,ax6], ['A','B','C','D','E','F']):
        ax.text(-.12,1.1,label, size=12, weight='bold', color='black', transform=ax.transAxes)

    if ICUcap != None:
        ax3.axhspan(ICUcap-2,ICUcap+2, xmax=0.92, facecolor=colors['ICUcap'], edgecolor=None, zorder=-1)
        ax3.text(1.0,0.3,'ICU capacity', size=7, color='red', rotation=-90, transform=ax3.transAxes)
        ax3.scatter(380,ICUcap, marker="<", color='grey')

    if full_wave != None:
        m = full_wave
        d = m.chopped_data().sum(axis=2)
        ax2.plot(t, m.rho*(d[:,4]+d[:,5]+d[:,6]), color='red', ls=':')
        ax3.plot(t, d[:,10]+d[:,11], color='red', ls=':')

    ax1.set_xticks([45, 135, 45+2*90, 45+3*90])
    ax1.set_xticklabels(['Oct.','Jan.','Apr.','July'])


    ax5_x = [1,offset+2,offset+4]
    ax5labels=['Initial','After\nwinter', 'After\none year']
    ax5.set_xticks(ax5_x)
    ax5.set_xticklabels(ax5labels)

    ax6_x = [offset+2,offset+4]
    ax6labels=['After\nwinter', 'After\none year']
    ax6.set_xticks(ax6_x)
    ax6.set_xticklabels(ax6labels)

    handles = [mpl.patches.Patch(facecolor=colors['highL'], edgecolor='black', label='By infection'),
               mpl.patches.Patch(facecolor=colors['high'], edgecolor='black', label='By vaccination')]
    ax5.legend(handles=handles, bbox_to_anchor=(0.1,0.9), ncol=1, frameon=False)

    fig.align_ylabels()

    if not silent: plt.show()
    if path!=None: fig.savefig(path)
