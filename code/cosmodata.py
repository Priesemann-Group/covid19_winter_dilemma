import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import gamma as gammafunc


germany_owid = pd.read_csv('../parameters/germany_ourworldindata.csv', sep=',', header=None, usecols=[3,18,45])
owid_ICU = germany_owid[18]
owid_NPI = germany_owid[45]
owid_dates = germany_owid[3]


ICUtime=[]
NPItime =[]
datesdict = {}
for times, dates in zip(range(len(owid_dates)-1), owid_dates[1:]):
    ICUtime.append(float(owid_ICU[1+times]))
    NPItime.append(float(owid_NPI[1+times])/100)
    datesdict[dates] = times
    
t=np.linspace(0,len(NPItime),len(NPItime))


feiernnoages = pd.read_csv('../parameters/PrivateFeiern_no_ages.csv', sep=',', header=None, usecols=[0,3])

datesCosmo = feiernnoages[0]

cosmotimeline=[]
for i in datesCosmo[1:]:
    cosmotimeline.append(datesdict[i])

avggroup = []

for i in range(len(datesCosmo)-1):
    avggroup.append((float(feiernnoages[3][i+1])-1)/4)
    
average_stringency = np.nanmean(NPItime)
average_cosmo = np.nanmean(avggroup)

NPItime_aligned = np.array(NPItime)/average_stringency*average_cosmo
    
#Plotten:
#ax.scatter(cosmotimeline,avggroup, color='red', alpha=0.5)
#ax.plot(t,NPItime, label='Stringency', color='green', zorder=5)



# ----------------------------------------- ROMANIA -----------------------------------------------

romania_owid = pd.read_csv('../parameters/romania_ourworldindata.csv', sep=',', header=None, usecols=[3,18,39])
                           
ROU_owid_ICU = romania_owid[18]
ROU_owid_vaccines = romania_owid[39]
ROU_owid_dates = romania_owid[3]

ROU_t=np.linspace(0,len(NPItime),len(NPItime))

ROU_ICUtime=[]
ROU_vaccinetime=[]
ROU_datesdict = {}
for times, dates in zip(range(len(ROU_owid_dates)-1), owid_dates[1:]):
    ROU_ICUtime.append(float(ROU_owid_ICU[1+times]))
    ROU_vaccinetime.append(float(ROU_owid_vaccines[1+times]))
    ROU_datesdict[dates] = times
    
#ax.plot(ROU_t, ROU_ICUtime)
#ax.plot(ROU_t, ROU_vaccinetime)
