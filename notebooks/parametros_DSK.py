import parametros


params = {
    'delta': 0.0017,
    'Theta':0.000429676,
    'Theta_ICU':0.09786686,
}

#December 2020
#serop=0.067
#30. December
#reported=0.0278
serop=0.173


# Values from OWD
V_raw = 725400.0
R_raw =  59800.0
cases =    147.03
ICU =       4.3
W_V =    22450.0
W_R =     11610.0
darkfigure = serop/(R_raw/1e6)

y0 = {
    'V': V_raw - W_V,
    'R': darkfigure * (R_raw-W_R),
    'W': W_V + darkfigure*W_R,
    'ICU': ICU,
}
y0['V']= V_raw - W_V - y0['R']*(1-(V_raw-W_V)/1e6)

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


