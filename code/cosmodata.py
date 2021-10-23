import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import gamma as gammafunc
from scipy.optimize import curve_fit

import parameters_agegroups


# Timeseries of the ICU 
germany_owid = pd.read_csv('../parameters/germany_ourworldindata.csv', sep=',', header=None, usecols=[3,18])
ICUtime = germany_owid[18]

#Cosmo data 
gespraeche = pd.read_csv('../parameters/Gespraeche.csv', sep=';', header=None, usecols=[0,3,7,11,15])

dates = germany_owid[3]


datesCosmo = gespraeche[0]


age1 = gespraeche[3]
age2 = gespraeche[7]
age3 = gespraeche[11]
age4 = gespraeche[15]


#GAMMA DISTRIBUTION----------------
alpha = parameters_agegroups.params_base['a_Rt']
beta = parameters_agegroups.params_base['b_Rt']
gamma_cutoff = round(parameters_agegroups.params_base['gamma_cutoff'])
x = -np.arange(-gamma_cutoff, 0, 1)


probdensity = beta**alpha*x**(alpha-1)*np.e**(-beta*x)/gammafunc(alpha)

#Day that the cosmo data starts (13.10.2020) compared to the start date of ourworldindata 
cosmostart = 264


# For loop calculates a timeseries of the perceived ICU 
Htime=[]
for times in range(len(dates[cosmostart:])):
    H = 0
    for i in range(gamma_cutoff):
        H = H + probdensity[i]*float(ICUtime[cosmostart+times-gamma_cutoff+i])
    Htime.append(H)
    

    

dict={}
for i,j in zip(dates[cosmostart:], Htime):
    dict[i] = j

ICU=[]
for i in datesCosmo[1:]:
    #print(i, dict[f'{i}'])
    ICU.append(dict[f'{i}'])
    
#Max ICU occupancy in the timeseries:
xmax=70


# COSMO DATA (because its a string so far...)

group1 = []
group2 = []
group3 = []
group4 = []

# Create arrays 
for i in range(len(datesCosmo)-1):
    group1.append(float(age1[i+1]))
    group2.append(float(age2[i+1]))
    group3.append(float(age3[i+1]))
    group4.append(float(age4[i+1]))

#------------------- AGE RELATED THINGS------------------------------


# Age structure of Germany (to calculate the average age of each age group)
germany = np.loadtxt('../parameters/germany.csv', delimiter=',')

# AVERAGE AGE 4 AGEGROUPS

avg_age1 = (np.array(germany[18:30,0])*np.array(germany[18:30,1])).sum()/np.sum(germany[18:30,1])
avg_age2 = (np.array(germany[30:50,0])*np.array(germany[30:50,1])).sum()/np.sum(germany[30:50,1])
avg_age3 = (np.array(germany[50:65,0])*np.array(germany[50:65,1])).sum()/np.sum(germany[50:65,1])
avg_age4 = (np.array(germany[65:75,0])*np.array(germany[65:75,1])).sum()/np.sum(germany[65:75,1])


# AVERAGE AGE 6 AGEGROUPS

avg_age1_6 = (np.array(germany[0:20,0])*np.array(germany[0:20,1])).sum()/np.sum(germany[0:20,1])
avg_age2_6 = (np.array(germany[20:40,0])*np.array(germany[20:40,1])).sum()/np.sum(germany[20:40,1])
avg_age3_6 = (np.array(germany[40:60,0])*np.array(germany[40:60,1])).sum()/np.sum(germany[40:60,1])
avg_age4_6 = (np.array(germany[60:70,0])*np.array(germany[60:70,1])).sum()/np.sum(germany[60:70,1])
avg_age5_6 = (np.array(germany[70:80,0])*np.array(germany[70:80,1])).sum()/np.sum(germany[70:80,1])
avg_age6_6 = (np.array(germany[80:,0])*np.array(germany[80:,1])).sum()/np.sum(germany[80:,1])

#Array with the average ages
ages = [avg_age1,avg_age2,avg_age3,avg_age4]

ages6 = [avg_age1_6,avg_age2_6,avg_age3_6,avg_age4_6,avg_age5_6,avg_age6_6]




grouparray = [group1,group2,group3,group4]
guesses = [[4.1,0.017], [4.3,0.015], [4.5,0.013], [4.7,0.011]] #initial guesses for the scipy fit 
agegroups = ['18-29', '30-49', '50-64', '65-74']

agegroups6=['0-20','20-40','40-60','60-70','70-80','80+']




ICUcap = parameters_agegroups.params_base['fit_ICUcap']  #the ICU beyond which no further reduction happens
epsilon = parameters_agegroups.params_base['fit_epsilon'] #smoothness of the function

#Shape of the function (first linear, then constant with a smooth transition)
def f(ICU, saturation, slope):
    return saturation - slope*epsilon*np.log(np.e**(1/epsilon*(ICUcap-ICU))+1)


plateaus = []
slopes = []


for g,p0,ag in zip(grouparray,guesses, agegroups): #Fits the function 
    popt, pcov = curve_fit(f, ICU, g, p0)
    
    
    plateaus.append(popt[0])
    slopes.append(popt[1])
    
    
    
# CONVERTING TO 6 AGEGROUPS 
    
def linf(ages,m,c): #linear function (see below) for the conversion to 6 agegroups 
    return ages*m + c

#Linear fits in order to transit from 4 to 6 agegroups (looks quite linear, see the notebook) 
popt_plateau, pcov = curve_fit(linf, ages, plateaus)

#calculating the plateaus for the new agegroups 
plateaus_6agegroups = linf(np.array(ages6),popt_plateau[0],popt_plateau[1])


popt_slopes, pcov = curve_fit(linf, ages, slopes)

#calculating the slopes for the new agegroups 
slopes_6agegroups = linf(np.array(ages6),popt_slopes[0],popt_slopes[1])
