import numpy as np
import pandas as pd

DEUparams = pd.read_csv('../parameters/DEU_params.csv', sep=';',header=0)


def get_params(country='Germany', modeltype='base', scenario='scenario2', nr_ages=6):
    #mapping_raw = {'Germany':raw_deu, 'Czech':raw_cze, 'Denmark':raw_dsk, 'Portugal':raw_por}
    #mapping_params = {'Germany':{}, 'Czech':params_cze, 'Denmark':params_dsk, 'Portugal':params_por}
    params = params_base.copy()
    #params.update(mapping_params[country])

    mapping_modelparams = {'base':{}, 'ramp':params_ramp, 'logistic':params_logistic}
    params.update(mapping_modelparams[modeltype])

    params['Rt_base'] = Rtbase[scenario]
    
    allinitials = initials(params, nr_ages)
    
    params.update({'y0':allinitials})
    
    for i in ['delta', 'Theta', 'Theta_ICU', 'gamma', 'gamma_ICU', 'alpha_R', 'alpha_u', 'alpha_w', 'chi_0', 'chi_1', 'u_base']:
        params.update({i: np.array(DEUparams[i])})
        
    params.update({'CM': CM_uniform})
    
    return params 


def initials(params, nr_ages):
    y0allages = []
    
    for ag in range(nr_ages):
        newparams = updateparams(params,ag)
        inits = calc_y0(newparams, ag)
        y0allages.append(inits)

    
    
    return np.array(y0allages).flatten(order='F')


def updateparams(params, agegroup):
    agespecificparams = {
        'delta': DEUparams.values[agegroup, 3],
        'Theta': DEUparams.values[agegroup, 4],
        'Theta_ICU': DEUparams.values[agegroup, 5],
        'gamma': DEUparams.values[agegroup, 6],
        'gamma_ICU': DEUparams.values[agegroup, 7],
        'alpha_R': DEUparams.values[agegroup, 8],
        'alpha_u': DEUparams.values[agegroup, 9],
        'alpha_w': DEUparams.values[agegroup, 10],
        'chi_0': DEUparams.values[agegroup, 11],
        'chi_1': DEUparams.values[agegroup, 12],
        'u_base': DEUparams.values[agegroup, 13],
    }
    
    params.update(agespecificparams)
    
    return params
    

def calc_y0(params, agegroup):
    Mi = DEUparams.values[agegroup, 1]
    
    V_raw = DEUparams.values[agegroup, 14]
    R_raw = DEUparams.values[agegroup, 15]
    cases = DEUparams.values[agegroup, 16]
    ICU_raw = DEUparams.values[agegroup, 17]
    W_V = DEUparams.values[agegroup, 18]
    W_R = DEUparams.values[agegroup, 19]
    serop = DEUparams.values[agegroup, 20]
    
    darkfigure = serop / (R_raw/1e6)


    RRv = darkfigure*(R_raw-W_R)

    V = V_raw - W_V - RRv*((V_raw-W_V)/(Mi*1e6))


    Wn = darkfigure*W_R
    Wv = W_V

    Etot = 1./params['rho']*cases*darkfigure
    Itot = 1./(params['gamma']+params['delta']+params['Theta'])* cases*darkfigure

    S = 1e6*Mi - Etot - Itot - ICU_raw - V - RRv - Wn - Wv

    print('Sapprox:', S)

    Eimmune =  (1-params['eta'])*(V)

    En = (Wn)/(S+Wn+Wv+Eimmune)*Etot
    Ev = (Wv+Eimmune)/(S+Wn+Wv+Eimmune)*Etot
    E = Etot - En - Ev

    In = (Wn)/(S+Wn+Wv+Eimmune)*Itot
    Iv = (Wv+Eimmune)/(S+Wn+Wv+Eimmune)*Itot
    I = Itot - In -Iv


    ICUimmune = (1-params['kappa'])*(In+Iv)

    ICUv = Iv*(1-params['kappa'])/(I+ ICUimmune)*ICU_raw
    ICU = ICU_raw - ICUv

    #Nur ne Approximation da Leute wegsterben in D
    Rimmune = S + Wn + Wv + (1-params['eta'])*V

    R = RRv * (S + Wn)/Rimmune
    Rv = RRv - R


    y0= {
            'V': V,
            'Wn': Wn,
            'Wv': Wv,
            'E': E,
            'EBn': En,
            'EBv': Ev,
            'I': I,
            'IBn': In,
            'IBv': Iv,
            'ICU': ICU,
            'ICUv': ICUv,
            'R': R,
            'Rv': Rv,
            'UC': V_raw,
            'WC': 0.,
            'D': 0.,
            'C': 0.,
        }


    S = Mi*1e6 - (sum(y0.values())-y0['UC'])


    y0.update({'S':S})
    
    y0_array = [y0['S'],y0['V'],y0['Wn'],y0['Wv'],y0['E'],y0['EBn'],y0['EBv'],y0['I'],y0['IBn'],y0['IBv'],
                y0['ICU'],y0['ICUv'],y0['R'],y0['Rv'],y0['UC'],y0['WC'],y0['D'],y0['C']]
    
    return y0_array


# Values from OWD
raw_deu = {
    'V_raw': 603900.0,
    'R_raw':  47400.0,
    'cases':    116.25,
    'ICU':       13.46,
    'W_V':    18000.0,
    'W_R':     9030.0,
    'serop':      0.112,
}

raw_cze = {
    'V_raw': 536500.0,
    'R_raw': 156600.0,
    'cases':     18.56,
    'ICU':        0.84,
    'W_V':    14650.0,
    'W_R':    29100.0,
    'serop':      0.51, #Feb and March 2021; 15. March reported = 0.13
}

raw_dsk = {
    'V_raw': 725400.0,
    'R_raw':  59800.0,
    'cases':    147.03,
    'ICU':        4.3,
    'W_V':    22450.0,
    'W_R':    11610.0,
    'serop':      0.173, #Dec 2020 serop=0.067, 30. Dec reported=0.0278
}

raw_por = {
    'V_raw': 755500.0,
    'R_raw': 102200.0,
    'cases':    191.4,
    'ICU':       12.88,
    'W_V':    18550.0,
    'W_R':    18400.0,
    'serop':      0.173,
}


params_base = {
    'Rt_base': 3.25,
    'Rt_free': 5,
    'eta': 0.75,                # vac effect against transmission
    'kappa': 0.8,               # vac effect against severe disease
    'sigma': 0.5,
    'gamma': 0.1,
    'gamma_ICU': 0.13,
    'delta': 0.0019,
    'rho': 0.25,        
    'omega_v_b': 1./(8*30),
    'omega_n_b': 1./(12*30),
    'chi_0': 0.1,       #
    'chi_1': 0.2,       #
    'alpha_w': 0.03,    #
    'alpha_u': 0.02,    #
    'alpha_R': 0.01,    #
    'e_R': 0.,
    'e_u': 0.,
    'e_w': 0.,
    'vac_max': 0.005,   # maximal vaccination rate?
    'u_base': 0.5,      # add w_base?
    'mu': 0.267,        # seasonality strength
    'd_0': 8*30.,       # seasonality start point
    'd_mu': 0.,         # day of seasonality maximum
    'a_Rt': 4.,
    'b_Rt': 0.7,
    'a_vac': 6.,
    'b_vac': 0.4,
    'gamma_cutoff': 45.,
    'tau_vac1': 6*7.,
    'tau_vac2': 2*7.,
    'Theta':0.0005366,      ###
    'Theta_ICU':0.09755,    ###
    't_max': 360.,
    'step_size': 0.1,
    'influx': 1,
    'epsilon_free':14.,
    'w_ICU': 0.,
}

params_cze = {
    'delta': 0.001717,
    'Theta':0.00040588,
    'Theta_ICU':0.0978747,
}

params_dsk = {
    'delta': 0.0017,
    'Theta':0.000429676,
    'Theta_ICU':0.09786686,
}

params_por = {
    'delta': 0.0019083,
    'Theta':0.000523277,
    'Theta_ICU':0.097566211,
}

params_ramp = {
    'time_u':7.,
    'time_w':7.,
}

params_logistic = {
    'r_base':0.015,
}


Rtbase = {
    'scenario1':5,
    'scenario2':3.75,
    'scenario3':2.5,
}

tmp = np.array([[0.11379518, 0.12928566, 0.09578086, 0.07168928, 0.04753613, 0.03652391],
[0.20671978, 0.13455128, 0.17471492, 0.09173655, 0.05707505, 0.0366261],
[0.09233393, 0.13444288, 0.17963023, 0.18572315, 0.11123063, 0.0440081],
[0.06793128, 0.05804828, 0.14403942, 0.38214947, 0.28185652, 0.12505345],
[0.05425213, 0.04132636, 0.10075469, 0.34296702, 0.35412762, 0.15467227],
[0.05417458, 0.03423004, 0.05199785, 0.19549228, 0.19778876, 0.69455006]])
CM_DEU = np.flip(np.flip(tmp, axis=0), axis=1)

CM_uniform = 1/6*np.ones([6,6])
