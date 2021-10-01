y0 = {
    'S': 285497.,
    'V': 525350.,
    'W':  30850.,
    'R':  156800.,
    'ICU': 3.,
}
y0['UC'] = y0['V']+y0['W']
y0['WC'] = y0['V']/100


params = {
    'delta': 0.001717,
    'Theta':0.00040588,
    'Theta_ICU':0.0978747,
}
