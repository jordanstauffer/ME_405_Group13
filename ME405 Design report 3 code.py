import numpy as np
import pyromat as pm
import math 

n = 10  # number of discrete sections
pp = 53 # number of parallel pipes
k_pipe= 160 #this is for copper, feel free to change for other materials 
Tw_in = 12 # C
Tw_out = 30 # C
T_steam_in = 120 # C
T_steam_sat = 40 # C
rho_steam_in = 0.1 #kg/m3
cp_cw= 4.186 #kJ/kgK
m_steam = 33 #33 #kg/s
Wall_thickness=0.01
D_inner_pipe=0.09 #meters total guess just need something feel free to change
Dout_inner_pipe=D_inner_pipe+Wall_thickness
D_outer_pipe=0.335
Dout_outer_pipe=D_outer_pipe+Wall_thickness
h_pipe=0.9 #found online must change for different materials 
R_w_unit = .05 # K/W cooling water surface
R_c_unit = 0.0003792 # K/W radial solid conduction
R_s_unit = .02 # K/W steam surface convection
R_ax_unit = (math.log(Dout_outer_pipe/D_inner_pipe))/(2*3.14*k_pipe*1) # K/W  axial conduction resistance 

#load the pyromat object for mixed-phase water
steam = pm.get('mp.H2O')
units = {'temperature':'C','pressure':'kPa'}
#Switch to Celsius & kPa
for k,v in units.items():
    pm.config['unit_'+k] = v


def condensing_hx(L_total,n,R_w_unit,R_c_unit,R_s_unit,R_ax_unit,m_water,cp_cw,T_steam_sat,Tw_in):
    L= L_total/n
    # surface resistance
    # The length scales are constant and the area gets smaller with more nodes
    # Thus the resistance value should be larger with increasing n
    Rw = [(R_w_unit+.5*R_c_unit)/L for i in range(n)]
    Rs = [(R_s_unit+.5*R_c_unit)/L for i in range(n)]
    # axial resistance
    # note that this resistance scales differently with n
    # cross sectional area stays the same and length gets 
    # shorter, thus resistance is smaller for smaller nodes 
    Rsolid = [R_ax_unit*L for i in range(n-1)] 

    A = np.zeros([2*n,2*n])
    b = np.zeros(2*n)

    C1 = m_water*cp_cw #fluid energy capacitance
    C2 = [1/Rw[i] for i in range(n)] # water surface heat transfer
    if abs(R_ax_unit)>0:
        C3 = [1/Rsolid[i] for i in range(n-1)] #solid conduction to next node
        C4 = [1/Rsolid[i] for i in range(n-1)] #solid conduction to prvious node
    else:
        C3 = [0 for i in range(n-1)]
        C4 = [0 for i in range(n-1)]
    C5 = [1/Rs[i] for i in range(n)] # air surface heat transfer


    b[0] = -C1*Tw_in
    for i in range(n):
        b[2*i+1] = -C5[i]*T_steam_sat
        A[2*i,2*i] += -C1
        A[2*i,2*i] += -C2[i]
        A[2*i,2*i+1] += C2[i]
        A[2*i+1,2*i] += C2[i]
        A[2*i+1,2*i+1] += -C2[i]
        A[2*i+1,2*i+1] += -C5[i]
        if i<(n-1): #all except last node
            A[2*i+1,2*i+1] += -C3[i]
            A[2*i+1,2*i+3] += C3[i]
        if i>0: #all except first node
            A[2*i,2*i-2] += C1
            A[2*i+1,2*i-1] += C4[i-1]
            A[2*i+1,2*i+1] += -C4[i-1]
    T = np.linalg.solve(A,b)
    T_cooling_water = [T[2*i] for i in range(n)]
    T_solid_pipe =[T[2*i+1] for i in range(n)]
    Q = m_water*cp_cw*(T_cooling_water[-1]-Tw_in)
    return Q, T_cooling_water, T_solid_pipe

def concentric_hx(L_total,m,R_w_unit,R_c_unit,R_s_unit,R_ax_unit,m_water,cp_cw,m_steam,cp_steam,T_steam_in,Tw_in):
    n = 1
    L= L_total/m
    # surface resistance
    # The length scales are constant and the area gets smaller with more nodes
    # Thus the resistance value should be larger with increasing n
    Rw = [(R_w_unit+.5*R_c_unit)/L for i in range(n)]
    Rs = [(R_s_unit+.5*R_c_unit)/L for i in range(n)]
    # axial resistance
    # note that this resistance scales differently with n
    # cross sectional area stays the same and length gets 
    # shorter, thus resistance is smaller for smaller nodes 
    Rc = R_ax_unit*L

    # Order of states is Tw_1, Ts_1, Tw2, Ts2....Ta1, Ta2, Ta3
    A = np.zeros([2*n*m+m,2*n*m+m])
    b = np.zeros(2*n*m+m)

    C1 = m_water*cp_cw #fluid energy capacitance
    C2 = [1/Rw[i] for i in range(n)] # water surface heat transfer
    if abs(Rc)<1e-6:
        C3 = [0 for i in range(m-1)]
        C4 = [0 for i in range(m-1)]
    else:
        C3 = [1/Rc for i in range(m-1)] #solid conduction to next node
        C4 = [1/Rc for i in range(m-1)] #solid conduction to previous node
    C5 = [1/Rs[i] for i in range(n)] # air surface heat transfer
    C6 = m_steam*cp_steam #ratio of energy capacity


    b[0] = -C1*Tw_in
    for j in range(m):
        #if counter-flow air
        k = 2*n*m+m-j-1
        #Co flow air
        #k = 2*n*m+j
        for i in range(n):
            row = 2*n*j+2*i
            col = 2*n*j+2*i
            #energy balance for the water at node i and pass j
            A[row,col] += -C1 #fluid capacitance
            if i>0 or j>0: #all except first node of first pass (handled in b[0])
                A[row,col-2] += C1 #fluid capacitance
            A[row,col] += -C2[i]#Heat transfer to/from solid
            A[row,col+1] += C2[i]#Heat transfer to/from solid
            
            #Energy balance for the solid at node i and pass m
            A[row+1,col+1] += -C2[i]#Heat transfer to/from water
            A[row+1,col] += C2[i]#Heat transfer to/from water
            A[row+1,col+1] += -C5[i]#Heat transfer to/from air
            A[row+1,k] += .5*C5[i]
            if j == (m-1): #last pass uses the inlet air temperature (no previous air temperature to use)
                b[row+1] += -.5*C5[i]*T_steam_in
            else:
                A[row+1,k-1] += .5*C5[i]
                
            #Axial conduction terms
            if i<(n-1) and j == (m-1): #all except last node of the last pass
                A[row+1,col+1] += -C3[i]#solid conduction to next node
                A[row+1,col+3] += C3[i]#solid conduction to next node
            if i>0 or j>0: #all except first node of first pass
                A[row+1,col+1] += -C4[j-1]#solid conduction to previous node
                A[row+1,col-1] += C4[j-1] #solid conduction to previous node

        #Energy balance for the steam temperature: Ts_i = Ts_i-1 + ((Tw_start_row - Tw_end_row)M_water*Cp)/M_steam*Cp
        A[k,k] += -C6 #Ta_i
        if j == (m-1): #last pass uses the inlet air temperature
            b[k] += -C6*T_steam_in
        else:
            A[k,k-1] += C6# Ts_i-1
        if j == 0:
            b[k] += -C1*Tw_in
        else:
            A[k,2*n*j-2] += C1#last node of previous pass
        A[k,2*n*(j+1)-2] += -C1 #last node of this pass

    T = np.linalg.solve(A,b)
    Twater = [T[2*i] for i in range(n*m)]
    Tsolid = [T[2*i+1] for i in range(n*m)]
    Tsteam = [T[2*n*m+i] for i in range(m)]
    Q = m_steam*cp_steam*(T_steam_in - Tsteam[-1])
    return Q, Twater, Tsolid,Tsteam

cp_inlet = float(steam.cp(T = T_steam_in, d= rho_steam_in)) # kJ/kg*K
lat_heat = float(steam.h(T = T_steam_sat,x=1)) - float(steam.h(T = T_steam_sat,x=0)) # kJ/kg

Q_superheat = m_steam*cp_inlet*(T_steam_in - T_steam_sat)#kg/s * kJ/kgK = kW
Q_latent = m_steam*lat_heat
Q_total = Q_superheat + Q_latent  
m_water = Q_total/(cp_cw*(Tw_out - Tw_in)) #kg/s



# Solve for length of superheated steam section (steam changes temperature)
T_cw_into_superheat_section = Tw_out-Q_superheat/(m_water*cp_cw)
dT1 = (T_steam_in - T_cw_into_superheat_section)
dT2 = (T_steam_sat - Tw_out)
LMTD_superheat = (dT1 - dT2)/np.log(dT1/dT2)
R_total_unit = R_w_unit+R_c_unit+R_s_unit
L_est = Q_superheat/(pp*LMTD_superheat/R_total_unit/1000) #m, could get closer initial guess using LMTD for the delta T instead of max delta T
L_guess = []
Q = []
for j in range(7,14):
    L = 0.1*j*L_est
    Qi, _, _,_ = concentric_hx(L,n,R_w_unit,R_c_unit,R_s_unit,R_ax_unit,m_water/pp*1000,cp_cw,m_steam/pp*1000,cp_inlet,T_steam_in,T_cw_into_superheat_section)
    L_guess.append(L)
    Q.append(Qi/1000)
#interpolate to find length
L_superheat = np.interp(Q_superheat/pp,Q,L_guess)
print (L_superheat)

# Solve for length of condensing section (steam remains constant)
dT1 = (T_steam_sat - Tw_in)
dT2 = (T_steam_sat - T_cw_into_superheat_section)
LMTD_condense = (dT1 - dT2)/np.log(dT1/dT2)
R_total_unit = R_w_unit+R_c_unit+R_s_unit
L_est = Q_latent/(pp*LMTD_condense/R_total_unit/1000)#m
L_guess = []
Q = []
for j in range(7,18):
    L = 0.1*j*L_est
    Qi, T_cooling_water, T_solid_pipe = condensing_hx(L,n,R_w_unit,R_c_unit,R_s_unit,R_ax_unit,m_water/pp*1000,cp_cw,T_steam_sat,Tw_in)
    L_guess.append(L)
    Q.append(Qi/1000)
#interpolate to find length
L_condense = np.interp(Q_latent/pp,Q,L_guess)
print (L_condense)

# put these two sections together
Q1, T_cooling_water, T_solid_pipe = condensing_hx(L_condense,n,R_w_unit,R_c_unit,R_s_unit,R_ax_unit,m_water/pp*1000,cp_cw,T_steam_sat,Tw_in)
Q2, Twater, Tsolid,Tsteam = concentric_hx(L_superheat,n,R_w_unit,R_c_unit,R_s_unit,R_ax_unit,m_water/pp*1000,cp_cw,m_steam/pp*1000,cp_inlet,T_steam_in,T_cooling_water[-1])

error = (((Q1+Q2)*pp/1000 - Q_total)/Q_total)*100
print('Length error is ' + str(error) + '%')

#Total heat transfer calculations 


print (LMTD_condense)

U_o=1/((((Dout_inner_pipe/2)/(h_pipe*D_inner_pipe))+(Dout_inner_pipe*math.log(Dout_inner_pipe/D_inner_pipe))/k_pipe)+(1/h_pipe))
Total_Q=U_o*3.14*(D_inner_pipe+D_outer_pipe)*L_condense*LMTD_condense
print (Total_Q)

Mass_flow_water=Total_Q/(cp_cw*LMTD_condense)
print (Mass_flow_water)



