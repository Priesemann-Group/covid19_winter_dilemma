import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import cosmodata 


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
    phis = np.array(list(map(model.get_phis, t, AGdata))).sum(axis=(2,3))
    shift = round(model.tau_vac1/model.step_size)
    d1a = np.roll(phis[:,0], -shift)
    d1a[-shift:] = 0
    shift = round(model.tau_vac1/2./model.step_size)
    d1b = np.roll(phis[:,0], -shift)
    d1b[-shift:] = 0
    shift = round(model.tau_vac2/model.step_size)
    d2 = np.roll(phis[:,1], -shift)
    d2[-shift:] = 0
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



def overview_agegroups(model, path=None, silent=False, arial=False):
    set_rcParams(arial=arial)
    t = model.times
    M = model.M
    data = model.chopped_data().sum(axis=2)
    AGdata = model.chopped_data()
    ags = AGdata.shape[2]

    colors = mpl.cm.viridis_r(np.linspace(0.,1.,ags))


    fig = plt.figure(figsize=(8, 9), constrained_layout=True)
    grid = fig.add_gridspec(ncols=3, nrows=5, hspace=0.1, wspace=0.1)
    

    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1],sharex=ax1)
    ax3 = fig.add_subplot(grid[2])
    ax4 = fig.add_subplot(grid[3],sharex=ax1)
    ax5 = fig.add_subplot(grid[4],sharex=ax1)
    ax6 = fig.add_subplot(grid[5],sharex=ax1)
    ax7 = fig.add_subplot(grid[6],sharex=ax1)
    ax8 = fig.add_subplot(grid[7],sharex=ax1, sharey=ax7)
    ax9 = fig.add_subplot(grid[8],sharex=ax1, sharey=ax7)
    ax10 = fig.add_subplot(grid[9],sharex=ax1)
    ax11 = fig.add_subplot(grid[10],sharex=ax1)
    ax12 = fig.add_subplot(grid[11],sharex=ax1)
    ax13 = fig.add_subplot(grid[12],sharex=ax1)
    ax14 = fig.add_subplot(grid[13],sharex=ax1)
    ax15 = fig.add_subplot(grid[14],sharex=ax1)

    axs = [ax1,ax2,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15]
#    titles = ['Incidence','Breakthrough share\nof incidence','',
#              'ICU occupancy','Breakthrough share\nnof ICU occupancy','Daily new deaths',
#              'Susceptible','Immune (vac.+natural)','Waned immunity',
#              'vac 1a','vac 2','vac total']

    # (S,V,Wn,Wv,E,EBn,EBv,I,IBn,IBv,ICU,ICUv,R,Rv,UC,WC,D,C)

    # Theta*I + (1-kappa)*Theta*(IBn+IBv) + Theta_ICU*(ICU+ICUv)
    dD = model.Theta*AGdata[:,7,:] \
         + (1-model.kappa)*model.Theta*(AGdata[:,8,:]+AGdata[:,9,:]) \
         + model.Theta_ICU*(AGdata[:,10,:]+AGdata[:,11,:])

    phis = np.array(list(map(model.get_phis, t, AGdata))).sum(axis=(2))
    shift = round(model.tau_vac1/model.step_size)
    d1a = np.roll(phis[:,0,:], -shift, axis=0)
    d1a[-shift:,:] = 0
    shift = round(model.tau_vac1/2./model.step_size)
    d1b = np.roll(phis[:,0,:], -shift, axis=0)
    d1b[-shift:,:] = 0
    shift = round(model.tau_vac2/model.step_size)
    d2 = np.roll(phis[:,1,:], -shift, axis=0)
    d2[-shift:,:] = 0
    def Rt(t):
        #Stapel aus 6x6 Matrizen, 4 lang. Elementweise gewichtete Addition und anschließende Summation über Zeilen
        return np.matmul( (np.moveaxis(model.Cs,0,2) * model.new_sr(t)).sum(axis=2), np.ones(6) )


    for ag in range(ags):
        ax1.plot(t,model.rho*(AGdata[:,4,ag]+AGdata[:,5,ag]+AGdata[:,6,ag])/(M[ag]/1e6), color=colors[ag])
        ax2.plot(t,(AGdata[:,6,ag])/(AGdata[:,4,ag]+AGdata[:,5,ag]+AGdata[:,6,ag]), color=colors[ag])
        # ax3 legend
        ax4.plot(t,(AGdata[:,10,ag]+AGdata[:,11,ag])/(M[ag]/1e6), color=colors[ag])
        ax5.plot(t,AGdata[:,11,ag]/(AGdata[:,10,ag]+AGdata[:,11,ag]), color=colors[ag])
#        ax6.plot(t,dD[:,ag]/(M[ag]/1e6), color=colors[ag])
        ax6.plot(t,AGdata[:,16,ag], color=colors[ag])
        ax7.plot(t,AGdata[:,0,ag]/M[ag], color=colors[ag])
        ax8.plot(t,(AGdata[:,1,ag]+AGdata[:,12,ag]+AGdata[:,13,ag])/M[ag], color=colors[ag])
        ax9.plot(t,(AGdata[:,2,ag]+AGdata[:,3,ag])/M[ag], color=colors[ag])
        ax10.plot(t,d1a[:,ag]/(M[ag]/1e6), color=colors[ag])
        ax11.plot(t,d2[:,ag]/(M[ag]/1e6), color=colors[ag])
        ax12.plot(t,(d1a[:,ag]+d1b[:,ag]+d2[:,ag]), color=colors[ag])
        ax13.plot(t,np.array(list(map(Rt,t)))[:,ag], color=colors[ag])
        ax14.plot(t,model.Gamma(t), color='black')
        ax15.plot(t,model.R0 * model.Gamma(t) * np.array(list(map(Rt,t)))[:,ag], color=colors[ag])

    for i,ax in enumerate(axs):
#        ax.set_title(titles[i])
        ax.set_ylim(0,None)

    ax1.set_ylabel("Daily new cases\nper million in age group")
    ax2.set_ylabel("Breakthrough share\nof daily new cases")
    #ax3 legend
    ax4.set_ylabel("ICU occupancy\nper million in age group")
    ax5.set_ylabel("Breakthrough share\nof ICU occupancy")
    ax6.set_ylabel("Cumulative deaths\nper million of population")
    ax7.set_ylabel("Susceptible fraction\nof the population")
    ax8.set_ylabel("Immune fraction\nof the population")
    ax9.set_ylabel("Waned immune fraction\nof the population")
    ax10.set_ylabel("Daily first-time vac.\nper million in age group")
    ax11.set_ylabel("Daily booster vac.\nper million in age group")
    ax12.set_ylabel("Daily total vac.\nper million of population")
    ax13.set_ylabel("Contacts")
    ax14.set_ylabel("Seasonality")
    ax15.set_ylabel("Total gross Rt")

    ax1.set_xticks([45, 135, 45+2*90, 45+3*90])
    ax1.set_xticklabels(['Oct.','Jan.','Apr.','July'])

    for ax in [ax2,ax5,ax7,ax8,ax9]:
        ax.set_ylim(0,1)

    for ax in axs:
        ax.axvline(180, ls=':', color='#ADB6B6FF', zorder=0)
        ax.set_xlabel('2021            2022')

#    for ax in [ax2,ax3,ax5,ax6,ax8,ax9,ax11,ax12]:
#        plt.setp(ax.get_yticklabels(), visible=False)

    # Build Legend in Panel 3
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    handles = [mpl.lines.Line2D([], [], color=colors[0], label='0-19'),
               mpl.lines.Line2D([], [], color=colors[1], label='20-39'),
               mpl.lines.Line2D([], [], color=colors[2], label='40-59'),
               mpl.lines.Line2D([], [], color=colors[3], label='60-69'),
               mpl.lines.Line2D([], [], color=colors[4], label='70-79'),
               mpl.lines.Line2D([], [], color=colors[5], label='80+')]
    ax3.legend(handles=handles, title='Age groups', loc='upper center', ncol=1, frameon=True)
    
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
        'FW':'#ED0000FF',
    }
    main_colors = [colors['high'],colors['mid'],colors['low']]
    main_colors_L = [colors['highL'],colors['midL'],colors['lowL']]



#    ax1.plot(t[1800:], np.ones(1800), color=colors['free'])
    def Rt(m,t):
        return max(np.linalg.eigvals((np.moveaxis(m.Cs,0,2) * m.fit_LR(t)).sum(axis=2)))
    
    for i,m in enumerate([m1,m2,m3]):
#        ax1.plot(t, np.array(list(map(m.Rt, t)))/m.R0, color=main_colors[i])
#        ax1.plot(t, np.array(list(map(m.Rt, t))), color=main_colors[i])
        ax1.plot(t, list(map(Rt, [m]*len(t), t)), color=main_colors[i])
        ax2.plot(t, m.rho*(data[i][:,4]+data[i][:,5]+data[i][:,6]), color=main_colors[i])
        ax3.plot(t, data[i][:,10]+data[i][:,11], color=main_colors[i])

        phis = np.array(list(map(m.get_phis, t, AGdata[i]))).sum(axis=(2,3))
        shift = round(m.tau_vac1/m.step_size)
        d1a = np.roll(phis[:,0], -shift)
        d1a[-shift:] = 0
        shift = round(m.tau_vac1/2./m.step_size)
        d1b = np.roll(phis[:,0], -shift)
        d1b[-shift:] = 0
        shift = round(m.tau_vac2/m.step_size)
        d2 = np.roll(phis[:,1], -shift)
        d2[-shift:] = 0
        ax4.plot(t, d1a+d1b+d2, color=main_colors[i])

    ax5.bar(1, data[0][0,1]/1e6, 0.5,
        align='center', color=colors['now'], edgecolor='black', zorder=-1)
    ax5.bar(1, (data[0][0,12]+data[0][0,13])/1e6, 0.5,
        align='center', color=colors['nowL'], edgecolor='black', zorder=-1, bottom=data[0][0,1]/1e6)

    offset = 0.5
    for i in [2,4]:
        for ab,m,j in zip([-0.5,0,0.5],[m1,m2,m3],[0,1,2]):
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


    ax1.set_ylim(0,None)
    ax2.set_ylim(0,None)
    ax3.set_ylim(0,None)
    ax4.set_ylim(0,None)
    ax5.set_ylim(0,1)
    ax6.set_ylim(0,None)

#    ax5.set_yticks([0.,0.25,0.5,0.75,1.])

    ax1.set_ylabel("Contact levels\ninfluenced by NPIs")
    ax2.set_ylabel("Daily new cases\nper million")
    ax3.set_ylabel("ICU occupancy\nper million")
#    ax4.set_ylabel("Reproduction number")
    ax4.set_ylabel("Daily vaccinations\nper million")
    ax5.set_ylabel("Immune fraction\nof the population")
    ax6.set_ylabel("Total deaths\nper million")

    #Panel 1
    ax1.text(0,1+0.05,'Scenario 1', size=7, color=colors['high'])
    ax1.text(0,0.75+0.05,'Scenario 2', size=7, color=colors['mid'])
    ax1.text(0,0.5+0.05,'Scenario 3', size=7, color=colors['low'])
    ax1.text(200,0.5+0.05,'No restrictions', size=7, color=colors['line'])
    
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
        ax2.plot(t, m.rho*(d[:,4]+d[:,5]+d[:,6]), color=colors['FW'], ls=':')
        ax3.plot(t, d[:,10]+d[:,11], color=colors['FW'], ls=':')
        ax2.text(0.58,0.72,'Full wave', size=7, color=colors['FW'], transform=ax2.transAxes)
        ax3.text(0.58,0.72,'Full wave', size=7, color=colors['FW'], transform=ax3.transAxes)

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
        
        
        
        
#-------------------------------------------------------------------------------------------------------------------        
        
        
        
def sixpanels_flexible(models, path=None, silent=False, arial=False, ICUcap=None, full_wave=None):
    set_rcParams(arial=arial)
    mpl.rcParams["legend.fontsize"] = 7

    m1, m2, m3 = models
    t = m1.times
    data = [m1.chopped_data().sum(axis=2), m2.chopped_data().sum(axis=2), m3.chopped_data().sum(axis=2)]
    AGdata = [m1.chopped_data(), m2.chopped_data(), m3.chopped_data()]
    
    ags = AGdata[0].shape[2]
    colors = mpl.cm.viridis_r(np.linspace(0.,1.,ags))

    fig = plt.figure(figsize=(6., 3.5), constrained_layout=True)
    grid = fig.add_gridspec(ncols=3, nrows=2, wspace=0.1)

    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1])
    ax3 = fig.add_subplot(grid[2])
    ax4 = fig.add_subplot(grid[3])
    ax5 = fig.add_subplot(grid[4])
    ax6 = fig.add_subplot(grid[5])

    colors = {
        'low':'#0099B4FF', 'mid':'#00468BFF', 'high':'#1B1919FF', 'free':'#1B1919FF',
        'lowL':'#FFFFFFFF', 'midL':'#FFFFFFFF', 'highL':'#FFFFFFFF',
#        'lowL':'#0099B499', 'midL':'#00468B99', 'highL':'#1B191999',
        'line':'#ADB6B6FF', 'ICUcap':'#FFAAAA',
        'now':'#ADB6B6FF', 'nowL':'#FFFFFFFF',
#        'now':'#93dfedFF', 'nowL':'#93dfed99',
        'FW':'#ED0000FF',
    }
    main_colors = [colors['high'],colors['mid'],colors['low']]
    main_colors_L = [colors['highL'],colors['midL'],colors['lowL']]


    def Rt(m,t):
        return max(np.linalg.eigvals((np.moveaxis(m.Cs,0,2) * m.fit_LR(t)).sum(axis=2)))
    
    def plot_axes_winter(ax):
        ax.axvline(180, ls=':', color=colors['line'], zorder=-4)
        ax.set_xlabel('2021           2022')
        ax.set_ylim(0,None)
        ax.set_xticks([45, 135, 45+2*90, 45+3*90])
        ax.set_xticklabels(['Oct.','Jan.','Apr.','July'])
        return None
    
    def plot_cosmo(ax):
        ax.scatter(cosmodata.cosmotimeline,cosmodata.avggroup, color='red', alpha=0.5, s=8)
        ax.plot(cosmodata.t,cosmodata.NPItime, label='Stringency', color='green', zorder=5)
        ax.set_ylabel("Stringency and\ncompliance")
        ax.set_xlabel('2020            2021')
        ax.set_xticks([cosmodata.datesdict['2020-04-15'], cosmodata.datesdict['2020-10-15'], 
                       cosmodata.datesdict['2021-04-15'], cosmodata.datesdict['2021-10-15']])
        ax.set_xticklabels(['Apr.', 'Oct.', 'Apr.', 'Oct.'])
        ax.axvspan(cosmodata.datesdict['2020-10-15'], cosmodata.datesdict['2020-12-15'], ymin=0, ymax=1, color='gray', alpha=0.2) 
        return None 
    
    def plot_romania(ax):
        startdate='2021-03-15'
        startpoint = cosmodata.datesdict[startdate]
        ax.plot(cosmodata.ROU_t[startpoint:], cosmodata.ROU_ICUtime[startpoint:])
        ax.plot(cosmodata.ROU_t[startpoint:], np.array(cosmodata.ROU_vaccinetime)[startpoint:]/100000*80)
        ax.set_ylabel("ICU occupancy and\ndaily vaccinations")
        ax.set_xlabel('2021')
        ax.set_xticks([cosmodata.ROU_datesdict['2021-04-15'], cosmodata.ROU_datesdict['2021-07-15'], 
                       cosmodata.ROU_datesdict['2021-10-15']])
        ax.set_xticklabels(['Apr.', 'July', 'Oct.'])
        ax.set_yticks([])
        return None 
    
    def plot_NPI(ax):
        for i, m in enumerate([m1,m2,m3]):
            ax.plot(t, list(map(Rt, [m]*len(t), t)), color=main_colors[i], zorder=-i)
        plot_axes_winter(ax)
        ax.set_ylabel("Contact levels\ninfluenced by NPIs")
        ax.text(0.54,0.15,'Lifting of\nrestrictions', size=7, color=colors['line'], transform=ax.transAxes)
        return None
    
    def plot_incidence(ax):
        for i, m in enumerate([m1,m2,m3]):
            ax.plot(t, m.rho*(data[i][:,4]+data[i][:,5]+data[i][:,6]), color=main_colors[i])
        ax.set_ylabel("Daily new cases\nper million")
        plot_axes_winter(ax)
        return None
    
    def plot_ICU(ax):
        for i, m in enumerate([m1,m2,m3]):
            ax.plot(t, data[i][:,10]+data[i][:,11], color=main_colors[i])
        ax.set_ylabel("ICU occupancy\nper million")
        plot_axes_winter(ax)
        return None
    
    def plot_vaccinations(ax):
        for i,m in enumerate([m1,m2,m3]):
            phis = np.array(list(map(m.get_phis, t, AGdata[i]))).sum(axis=(2,3))
            shift = round(m.tau_vac1/m.step_size)
            d1a = np.roll(phis[:,0], -shift)
            d1a[-shift:] = 0
            shift = round(m.tau_vac1/2./m.step_size)
            d1b = np.roll(phis[:,0], -shift)
            d1b[-shift:] = 0
            shift = round(m.tau_vac2/m.step_size)
            d2 = np.roll(phis[:,1], -shift)
            d2[-shift:] = 0
            ax.plot(t, d1a+d1b+d2, color=main_colors[i])
        ax.set_ylabel("Daily vaccinations\nper million")
        plot_axes_winter(ax)
        return None 
    
    def plot_bar_immunity(ax):
        ax.bar(1, data[0][0,1]/1e6, 0.5, align='center', color=colors['now'], edgecolor='black', zorder=-1)
        ax.bar(1, (data[0][0,12]+data[0][0,13])/1e6, 0.5, align='center', 
               color=colors['nowL'], edgecolor='black', zorder=-1, bottom=data[0][0,1]/1e6)
        offset = 0.5
        for i in [2,4]:
            for ab,m,j in zip([-0.5,0,0.5],[m1,m2,m3],[0,1,2]):
                ax.bar(offset+i+ab, data[j][900*i-1,1]/1e6, 0.5,  
                    align='center', color=main_colors[j], edgecolor='black', zorder=-1)
                ax.bar(offset+i+ab, (data[j][900*i-1,12]+data[j][900*i-1,13])/1e6, 0.5,
                    align='center', color=main_colors_L[j], edgecolor='black', zorder=-1,
                        bottom=data[j][900*i-1,1]/1e6)
        ax.set_ylim(0,1)
        ax.set_ylabel("Immune fraction\nof the population")
        ax_x = [1,offset+2,offset+4]
        axlabels=['Initial','After\nwinter', 'After\none year']
        ax.set_xticks(ax_x)
        ax.set_xticklabels(axlabels)
        handles = [mpl.patches.Patch(facecolor=colors['highL'], edgecolor='black', label='By infection'),
               mpl.patches.Patch(facecolor=colors['high'], edgecolor='black', label='By vaccination')]
        
        #ax.legend(handles=handles, loc='upper right', ncol=1, frameon=False)
        return None 
    
    def plot_bar_deaths(ax):
        offset = 0.5
        for i in [2,4]:
            for ab,m,j in zip([-0.5,0,0.5],[m1,m2,m3],[0,1,2]):
                ax.bar(offset+i+ab, data[j][900*i-1,16], 0.5, 
                    align='center', color=main_colors[j], edgecolor='black', zorder=-3)    
        ax.set_ylim(0,None)
        ax.set_ylabel("Total deaths\nper million")
        ax_x = [offset+2,offset+4]
        axlabels=['After\nwinter', 'After\none year']
        ax.set_xticks(ax_x)
        ax.set_xticklabels(axlabels)
        return None
    
        
    
    
    # -----------------------------------------Plotting ----------------------------------------------
    
    plotting_dict={
    ax1: plot_cosmo,
    ax2: plot_NPI,
    ax3: plot_NPI,
    ax4: plot_romania,
    ax5: plot_incidence,
    ax6: plot_bar_immunity
    }
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:   
        plotting_dict[ax](ax)


    
    # --------------------------------------- General Axis Things -----------------------------------
    

    for ax, label in zip([ax1,ax2,ax3,ax4,ax5,ax6], ['A','B','C','D','E','F']):
        ax.text(-.12,1.1,label, size=12, weight='bold', color='black', transform=ax.transAxes)
        
    # Not touched yet:

    if ICUcap != None:
        ax3.axhspan(ICUcap-2,ICUcap+2, xmax=0.92, facecolor=colors['ICUcap'], edgecolor=None, zorder=-1)
        ax3.text(1.0,0.3,'ICU capacity', size=7, color='red', rotation=-90, transform=ax3.transAxes)
        ax3.scatter(380,ICUcap, marker="<", color='grey')

    if full_wave != None:
        m = full_wave
        d = m.chopped_data().sum(axis=2)
        ax2.plot(t, m.rho*(d[:,4]+d[:,5]+d[:,6]), color=colors['FW'], ls=':')
        ax3.plot(t, d[:,10]+d[:,11], color=colors['FW'], ls=':')
        ax2.text(0.58,0.72,'Full wave', size=7, color=colors['FW'], transform=ax2.transAxes)
        ax3.text(0.58,0.72,'Full wave', size=7, color=colors['FW'], transform=ax3.transAxes)


    fig.align_ylabels()

    if not silent: plt.show()
    if path!=None: fig.savefig(path)
