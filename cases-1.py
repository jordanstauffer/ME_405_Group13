import pyromat as pm
import numpy as np
import matplotlib.pyplot as plt
import numpy_financial as npfin

Tmax = 625 #celcius
Tsat = 295.01
Pmax = 15*1000 #kPa
Pcond = 8.209798435 #kPa
Preheat = 3*1000 #kPa

power = 35 #megaWatts
plantprice = 50000000 #dollars
#over 25 years

fuelcost = [3.66,3.54,3.34,2.96,2.86,2.72,2.55,2.91,3.16,2.35,3.93,4.07] ###natural gas cost $/1000ft^3
sunhours = [5,6,7,9,14,16,15,13,10,8,7,6] # at 800W/m^2
saleprice = [47,47.3,48.6,47.7,46.8,48.3,53.5,52.8,58.2,57.2,56.8,53.8] #$/MWh

days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
hours = [ds*24 for ds in days]
hourssun = [[ds*sunhours[i] for i in range(len(sunhours))] for ds in days] #hours of operation per month

inflationsale = 0.03
inflationcost = 0.05

nt = [0.8,0.85,0.9,0.95]
highpress = [1500000,2000000,2500000,3000000]
midpress = highpress
lowpress = midpress
#combine all costs for no reheat, combine high and mid for 1st reheat and low for 2nd reheat

Tcond = [45,35,30,25]
condensorcost = [12000000,25000000,30000000,40000000]

npump = [0.75,0.8,0.85,0.9]
pumpcost = [25000,30000,35000,40000]

#####84 turbine options#######
#4 options with no reheat: use combined cost of high, mid, and low
#16 options with 1 reheat: use 0.1 to 4 MPa for Preheat. use high and mid cost for 1st stage (boiler, turbine), and low cost for 2nd stage (boiler, turbine)
                           # 
#64 options with 2 reheat: 1st stage is 1-5 MPa, 2nd stage is 0.1-0.5MPa. use high and mid cost for 1st stage (boiler,trurbine), and low cost for 2nd stage (boiler,turbine)
                           #
#open feedwater heaters, if included, must be same pressure as reheat. can include condensor temp and pump efficiency

steam = pm.get('mp.H2O')
units = {'temperature':'C','pressure':'kPa'} #change units
for k,v in units.items():
    pm.config['unit_'+k] = v

#boiler
def isobaric_heat_addition(T,P):
    Temp = T
    Pressure = P
    h = steam.h(T=Temp,p=Pressure)
    s = steam.s(T=Temp,p=Pressure)
    v = steam.v(T=Temp,p=Pressure)
    x=None
    return Temp,Pressure,h,s,v,x

#turbine
def adiabatic_expansion(nt,h,P,s):
    Pressure = P
    sf,sg = steam.ss(p=Pressure)
    hf,hg = steam.hs(p=Pressure)
    snew=s 
    if np.any(snew>sg):
        hs = steam.h(s=snew,p=P)
        hnew = h-nt*(h-hs)
        Temp = steam.T(h=hnew,p=Pressure)
    else:
        x = (snew-sf)/(sg-sf)
        hs = hf+x*hg
        hnew = h-nt*(h-hs)
        Temp = steam.Ts(p=Pressure)
              
    if np.any(hnew>hg):
        v = steam.v(T=Temp,p=Pressure)
        snew = steam.s(T=Temp,p=Pressure)
        x = None

    else:
        x = (hnew-hf)/(hg-hf)
        v = steam.v(p=Pressure,x=x)
        snew = steam.s(p=Pressure,x=x)
    return Temp, Pressure, hnew,snew,v,x

#condensor
def get_Saturated_state(T,P):
    Temp=T
    Pressure=P
    hf,hg = steam.hs(Temp)
    sf,sg = steam.ss(Temp)
    vf,vg = steam.vs(Temp)
    s=sf
    h=hf
    x=0
    v=vf
    return Temp, Pressure,h,s,v,x

#pump
def adiabatic_compression(Pf,Pmax,npump,v,h):
    Pressure=Pmax
    hs = h+v*(Pressure-Pf) 
    hnew = h+(hs-h)/npump 
    Temp = steam.T(h=hnew,p=Pressure) #may be wrong
    vnew = v
    s=steam.s(T=Temp,p=Pressure)
    x=None
    return Temp,Pressure,hnew,s,vnew,x

def mixing(T,P):   #for feedwater
    Pressure = P
    Temp = T
    hf,hg= steam.hs(T=Temp)
    sf,sg= steam.ss(T=Temp)
    vf,vg = steam.vs(T=Temp)
    h = hf
    s=sf
    v=vf
    x=0
    return Temp,Pressure,h,s,v,x

#1 reheat stage - nt will change    
T1,P1,h1,s1,v1,x1 = get_Saturated_state(Tcond,Pcond) #condensor
T2,P2,h2,s2,v2,x2 = adiabatic_compression(P1,Preheat,npump,v1,h1) #pump 1
T3,P3,h3,s3,v3,x3 = isobaric_heat_addition(Tmax,Pmax) #boiler 1            STAGE 1
T4,P4,h4,s4,v4,x4 = adiabatic_expansion(nt,h3,Preheat,s3) #turbine 1
T5,P5,h5,s5,v5,x5 = isobaric_heat_addition(Tmax,Preheat)#boiler 2
T6,P6,h6,s6,v6,x6 = adiabatic_expansion(nt,h5,Pcond,s5)#turbine 2

def special_calculations(power,h1,h2,h3,h4,h5,h6):
    specific1 = h2-h1 #work, pump
    specific2 = h3-h2 #heat, boiler
    specific3 = h4-h3 #work, turbine
    specific4 = h5-h4 #heat, boiler

    flow = (power*100)/(-1*specific3)
    
    total1 = flow*specific1
    total2 = flow*specific2
    total3 = flow*specific3
    total4 = flow*specific4

    totalheat = total4+total2

    n_total = power/totalheat*1000
    return n_total,totalheat
n_total,totalheat= special_calculations(power,h1,h2,h3,h4,h5,h6)

#ask prof how to plot in month form. struggling for some reason
#ADD TURBINE, CONDENSOR, and PUMP COSTS ONCE CHOSEN
cos_fuel= [[((fuelcost[i]*hours[i]*totalheat*1000*3600)/(10**9*1.055))*(1+inflationcost)**(j) for i in range(len(hours))] for j in range(25)]

cos_fuel80 = [item[0] for item in cos_fuel]
cos_fuel85 = [item[1] for item in cos_fuel]
cos_fuel90 = [item[2] for item in cos_fuel]
cos_fuel95 = [item[3] for item in cos_fuel]

totalrevenue = [[power*saleprice[i]*8760*(1+inflationsale)**(j) for i in range(len(saleprice))] for j in range(25)]


#####purchase cost is the dollar amount for the turbines. will need to be changed to include condensor and pump when applicable#####
# netcost80 = [-1*purchasecost[0]]
# netcost80.extend([float(-1*cos_fuel80[i]+totalrevenue[i]) for i in range(len(totalrevenue))])
# netcost85 = [-1*purchasecost[1]]
# netcost85.extend([float(-1*cos_fuel85[i]+totalrevenue[i]) for i in range(len(totalrevenue))])
# netcost90 = [-1*purchasecost[2]]
# netcost90.extend([float(-1*cos_fuel90[i]+totalrevenue[i]) for i in range(len(totalrevenue))])
# netcost95 = [-1*purchasecost[3]]
# netcost95.extend([float(-1*cos_fuel95[i]+totalrevenue[i]) for i in range(len(totalrevenue))])

# internalrr80 = npfin.irr(netcost80)
# internalrr85 = npfin.irr(netcost85)
# internalrr90 = npfin.irr(netcost90)
# internalrr95 = npfin.irr(netcost95)
# print(internalrr80,internalrr85,internalrr90,internalrr95)


#######order of operations#########
# 1. best condenser temperature
# 2. best pump efficiency
# 3. best number of reheats/efficiencies

#####or #####
# 1. best number of reheats/efficiencies
# 2. best condensor temp, 
# 3. best pump efficiencies
# 3. redo number of reheats/efficiencies