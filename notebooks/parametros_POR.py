import parametros


params = {
    'delta': 0.0019083,
    'Theta':0.000523277,
    'Theta_ICU':0.097566211,
}




# Values from OWD
V_raw = 755500
R_raw =  102200.0
cases =    191.4
ICU =       12.88
W_V =    18550
W_R =     18400.0
darkfigure = 2
y0 = {
    'V': V_raw - W_V,
    'R': darkfigure * (R_raw-W_R),
    'W': W_V + darkfigure*W_R,
    'ICU': ICU,
}


E_stay = 1./parametros.params['rho']
I_stay = 1./(parametros.params['gamma']+parametros.params['delta']+parametros.params['Theta'])

immune_infected = (y0['R']+y0['V'])*(1-parametros.params['eta'])
S_approx = 1e6 - sum(list(y0.values())) - darkfigure*(E_stay+I_stay)*cases
breakthrough = immune_infected / (immune_infected + S_approx)

y0.update({
    'E': darkfigure*cases*E_stay*(1-breakthrough),
    'EB': darkfigure*cases*E_stay*breakthrough,
    'I': darkfigure*cases*I_stay*(1-breakthrough),
    'IB': darkfigure*cases*I_stay*breakthrough,
})

y0['S'] = 1e6 - sum(list(y0.values()))

y0.update({
    'UC': V_raw,
    'WC': 0.,
    'D': 0.,
    'C': 0.,
})


y0_array = [y0['S'],y0['V'],y0['W'],y0['E'],y0['EB'],y0['I'],y0['IB'],y0['ICU'],y0['R'],y0['UC'],y0['WC'],y0['D'],y0['C']]
params['y0'] = y0_array
