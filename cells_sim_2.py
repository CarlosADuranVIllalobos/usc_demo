"""
cells_sim_2.py
Script contains functions to run cell simulations
Created on Tue Jun 20 16:42:31 2023
@author:Carlos Alberto Duran Villalobos-The University of Manchester
"""
import numpy as np
from scipy.special import expit
from scipy.integrate import odeint
import math

def tcellmodel(y, t, us, F_air, rgvar, bivar): #ODE Model 

    #Variables
    X = y[0]  # Viable cell count
    X_dead = y[1]  # Dead cell count
    G = y[2]  # Glucose conc. 
    I = y[3]  # Inhibitory factor
    V = y[4]  # Volume
    DO = (y[5]*6.5/100)/1e3  # %->mg/mL 
    pH = y[6]
    F=us    # Feed rate for inputs [glucose + ]
    
    #Constant rates
    Co = 40      #Glucose concentration in the the feeding vessel [mg/mL]
    Vs = 0     #Sampling volume [mL]
    if I<0:
        rg = 0.02/5 + rgvar/5  #growth rate lag phase h-1
        rd = 0
    else:
        rg = 0.02 + rgvar  #growth rate exponential phase h-1
        rd = 0.007 #Death rate by exhaustion h-1
        
    
    rdg = 0.083 #death rate caused by lack of glucose  h-1
    rc = 0.05     #constant specific consumption rate of glucose mg *1 x 10^6 cell/hr  
    ri = 1     #inhibition units created 1x10^6 cell*h-1
    #growth inhibition
    bi = 95.260 + bivar
    ai =  0.111
    #Death promotion
    ad = 150
    #bp = 80 #Threshold of death promotion
    #ap = 0.193 #Sensitivity of death promotion
    #Constants glucose
    kp = 2.710    #Monod promotion of T-cell growth rate by glucose half-velocity mg/mL
    bgluc = 0.05     #Decay Glc Thresh mg/mL
    agluc = 30       #Decay Glc Taut 
    kgu = 0.2     #Monod Glucose use half-velocity mg/mL 
    #Constants DO
    a_kla = 0.0001#0.0002    #Vessel geometry constant
    m_o2 = 0.0012    #Maintenance rate
    a = 0.38     #exponent of kla for superficial gas velocity V_s
    b = 0.34     #exponent of kla for required power over volume P_t1 / V_m
    c = -0.38    #exponent of kla for viscosity
    r = 0.01405                 #Tank radius [m]
    n_imp = 1                   #number of impellers
    Po = 2.15                      #unaerated impeller power number
    N = 450/60                     #impeller [rpm] to [rps]
    D = 0.0114                  #Diameter of the impeller [m]
    R = 8.314                   #Universal gas constant [J/Kmol ]
    T = 310.15                        #Temperature [K] assume constant 37[C]
    V_m = V/1e6                         # Volume ml to [m^3]
    V_s = (F_air)/(3.1416*(r**2));      #superficial gas velocity [m/s]
    eps = 0.1                           #gas hold-up coefficient, assumed to be constant
    pressure = 1                        # asume constant [bar]
    viscosity = 0.95                    #asume constant [cp = mPa/s]
    vis_scaled =(viscosity / 100)      # Scaled viscosity [cP]
    Henry_c = 0.0251                    #Henry equilibrium constant [bar L/mg]
    O_2_in = 0.21                       #oxygen concentration in the air entering the tanh air= 0.21
    #O_2_in = 0.21*6.5                  #.% to g/mL
    #gamma = 1.85e-6-1.5e-6                     #propotionality constant pH [ml/cells]
    gamma = -5.899e-09                 #propotionality constant pH [ml/cells]
    k_o2 =  1.145e7                   #
    #Cab = 3e-3                                #concentration acid and mase mol/l
    bf=1e-14     
    
    #Equations glucose
    mx = G/(kp+G)                       #Promotion in growth by glucose   
    migluc = 1 - expit(agluc*(G-bgluc)) #Inhibition in decay by glucose, is defined as expit(x) = 1/(1+exp(-x))
    mg =G/(kgu+G)                       #promotion of glucose use 
    dG = -rc*mg*X/1e6 - G*F/V + Co*us/V  #Equations inhibition
    m_o2=DO/(k_o2+DO)                    #Promotion of oxygen use
 
   
    #mp = expit(ap*(I-bp))       #Promotion in cell death by unknown Inhibitor
    mi = 1 - expit(ai*(I-bi))  
    mp = I/(ad+I)
    #Equations DO
    
    ##Calculating liquid height in vessel
    Z = V_m/(3.1416*(r**2))   #[m]
    Z = Z*(1-eps)        # ungassed height
    ##broth density
    pho_b = 1.005e3 + X*20e-12*1e6*1e-3  + X_dead*20e-12*1e6*1e-3 # [cells/ml *20picograms  ml/m3 kg/g] 
                                                            #Broth density RPMI 10% is 1.005 g/cm3 = 1.006e6 g/m3
    ##necessary power 
    unaerated_power = n_imp * Po * pho_b * (N**3) * (D**5)
    P_g = 0.706 * (((unaerated_power**2) * N * D**3) / (F_air**0.56))**0.45 # gassed power consumption (Pag)
    P_n = P_g/unaerated_power; 
    variable_power = (n_imp * Po * pho_b * (N**3) * (D**5) * P_n)/1000;
    
    ##Calculating log mean pressure of vessel [bar]= 100 000 [N/m²]=[kg/mseg2]   density[g/m3]*1e-3=[kg/m3]
    pressure_bottom  =  1 + pressure + pho_b * Z * 9.81*1e-5 #* 1e-5 #density*[m2/seg2]    # [bar]
    pressure_top = 1 + pressure     # [bar]
    log_mean_pressure = (pressure_bottom - pressure_top)/(math.log(pressure_bottom/pressure_top))
    total_pressure = log_mean_pressure
    do_star =  (((total_pressure) * O_2_in) / Henry_c)/1e3     #in mg/l -> mg/mL(Henry's constant has the units of (bar L mg-1)
    #do_star = do_star*100/6.5       ##in mg/l -> % , saturation in RPMI at 37C approx. do_star= 6.71 mg/L=.00671 mg/mL
    ##Calculate required power
    P_air = (V_s * R * T * V_m  / 22.4 * Z) * math.log(1 + pho_b * 9.81 * Z /(pressure_top*10^5)) # aeration power rate [kW]
    P_t1 = variable_power + P_air #required power [kW]
    ##calculate kLa: Oxygen transfer coefficientin AMBR15 manual 2-6
    kla = a_kla * V_s**a * (P_t1 / V_m)**b * vis_scaled**c #[1/h]
    #Equation H+ effect in cell growth, other change is not included due to lack of experimental data
    #f_H = 1/(1+(k1/H)+(H/k2))
    
    #Equation hydrogen ion concentration [H+] in liters
 
    #B = ((bf/H-H)*1e-3*V - (Cab*(1e-3*Fa -1e-3*Fb)*dt))/(1e-3*V+(1e-3*Fa+1e-3*Fb)*dt) #(mL/h)(1L/1000ml)
    #dH = 10**(-gamma*(X-F*X/V))+ ((-B+math.sqrt(B**2+4*bf))/2-H)/dt
    dH = 10**(-gamma*(X-F*X/V))#
 
    F=0 #Consider no changes in volume, perfusion
    dV = F-Vs 
    #critical values of oxygen (vardar and lilly)
    #mo = 0.5*(1−np.tanh(Ainhib(DO2critP−DO2)) #no data but it can be included as a modulator
    dDO =  - m_o2*X  + kla*(do_star-DO) - DO*dV/V  #[mg/mL] 
    dDO = dDO * 100 * 1e3 /6.5   #[mg/mL] to %
 
    #dX = rg*X*mx - rd*X*mp - rdg*X*migluc -X*dV/V
            
    dX = rg*X*mx*mi  - rd*X*mp -rdg*X*migluc-dV*X/V
    dX_dead = rd*X*mp + rdg*X*migluc -X_dead*dV/V
    dI = ri*(X)/1e6 
    #dH = 10**(-X*gamma)
    dpH = -np.log10(dH)#+Cab*Fb*24/(.03*40+24)*350-Cab*Fa*.03*40/(.03*40+24)*50
    
    return [dX, dX_dead, dG, dI, dV, dDO, dpH]

# Define the control and simulation functions
def sat(qc, qc_min, qc_max):                          # function to return feasible value of qc
    return max(qc_min,min(qc_max,qc))

def PID(Kp, Ki, Kd, u_min, u_max, MV_bar=0):
    
    # initialize stored data
    e_prev = 0
    t_prev = -100
    I = 0
    
    # initial control
    MV = MV_bar
    
    while True:
        # yield MV, wait for new t, PV, SP
        t, PV, SP = yield MV
        
        # PID calculations
        e = SP - PV
        
        P = Kp*e
        I = I + Ki*e*(t - t_prev)
        D = Kd*(e - e_prev)/(t - t_prev)
        
        MV = MV_bar + P + I + D
        MV = sat(MV, u_min, u_max)
        
        # update stored data for next iteration
        e_prev = e
        t_prev = t

def tcell_loop(model,ctrl_param):   
    #Function that runs the ODE loop 
    tf = ctrl_param['hours']        
    z0 = ctrl_param['Y0'].copy() 
    dilution = ctrl_param['dilution'].copy()  
    ## DO PID control
    DO_sp = 100 # setpoint
    DO = z0[5] #Initial DO
    # do simulation at fixed time steps dt
    dt = 1
    ti = 0.0   
    # control saturation
    u_min = 1e-12                            # minimum possible flowrate
    u_max = 0.91 #* 60                         # 0.91 ml/min *60 min/1h maximum possible flowrate
    u = 0                                    #Initial flowrate
    # control parameters
    kp = 0.015 #0.00015#0.015
    ki = 0.05#0.005#.05
    kd = 0
    do_pid = PID(kp, ki, kd, u_min, u_max)        # create pid control
    do_pid.send(None)              # initialize

    ## Solve ODE
    Y = np.array([np.append(z0[:],u)]) #Create an extra dimension in the array [[]]
    for t in np.linspace(ti,tf,int((tf-ti)/dt)): #get simulation time   
        # PID control calculations
        u = do_pid.send([t, DO, DO_sp]) # compute manipulated variable

        if dilution[1]==t: #Add dilution
            z0[0]=dilution[0]*z0[0] + (1 - dilution[0])*0 #Y0[0] if cells in the new media
            z0[1]=dilution[0]*z0[1] + (1 - dilution[0])*0 #Y0[0] if cells in the new media
            z0[2]=dilution[0]*z0[2] + (1 - dilution[0])*z0[2] #glucose
            #z0[3]=dilution[0]*z0[1] + (1 - dilution[1])*0  #Inhibitor
            z0[4]=z0[4]+(1-dilution[0])*z0[4]
            z0[5]=dilution[0]*z0[5] + (1 - dilution[0])*z0[5] #DO
            #z0[6]=dilution[0]*z0[6] + (1 - dilution[0])*Y0[6] #pH
         
        z = odeint(model, z0,[t,t+dt], args= (ctrl_param['Recipe_Fs_sp'][int(t)], u,ctrl_param['rg_var'],ctrl_param['bi_var']))
        #physical constraints
        if z[-1][0]<1000:      #viable cells
            z[-1][0]=0
        if z[-1][1]<0:      #dead cells
            z[-1][1]=0
        if z[-1][2]<0:      #glucose conc.
            z[-1][2]=0
        if z[-1][5]<0:      #DO
            z[-1][5]=0
        if z[-1][6]<0:      #pH
            z[-1][6]=0
        
        Y = np.concatenate((Y[:],np.array([np.append(z[-1],u)])),axis=0) #Create an extra dimension in the array [[]]
        z0 = z[-1] 
        DO = z0[5]
    return Y

def control_tcell(model,hours,us,variability):
    #variability: [viability, inhibitor, rg, bi]                
    dilution=[0,160] #[Proportion,time] , no dilution
    Y0 = np.array ([1.1e6, 0, 0, -24, 50, 100, 7.4 ]) #[viable cells, dead cells, glucose c.,Inhibitor, Volume, DO, pH]
    viab = variability [0]
    Y0[1] = Y0[0]*(1-viab)
    Y0[0] = Y0[0]*viab
    Y0[3] = Y0[3] +  variability[1]*48
    rg_var = 0.02189*variability[2]
    bi_var = 95.26*variability[3]+48
    #Function for control purpose/simulates sensors   
    controlstep = 1 
    num = int(hours/controlstep)+1
    ctrl_param = {
            'Y0':Y0,
            'Recipe_Fs_sp':us ,
            'hours':hours,
            'dilution':dilution,
            'rg_var':rg_var,
            'bi_var':bi_var,         
            } 
    sol = tcell_loop(model,ctrl_param)       
    sy=sol[:,0]
    sx=sol[:,1:]
    su=us  # Extra variable for future control
    return sx,sy,su

