#!/usr/bin/env python
# coding: utf-8

# The spiking neuron network presented here has been adapted from (Torao-Angosto et al., 2021).
# The network is made of interconnected leaky integrate-and-fire excitatory and inhibitory neurons.
# Excitatory neurons incorporate the activity-dependent fatigue mechanism modeling spike-frequency adaptation (SFA), 
# which is a key ingredient to generate Up-Down slow oscillations typical of deep sleep and anesthesia.
# What follows is the python code, and associated explanations, to generate dynamics for several adaptation and excitation levels.
# 

import nest
import time
from numpy import exp
import numpy
import math
import random
import sys
from scipy.io import savemat
import matplotlib.pyplot as plt
import numpy as np 


# We first set a few parameters
RelaxTime = 5000;                               # [ms] duration of the initial spontaneous dynamics. This duration should be
                                                #      sufficiently long to ensure both excitatory and inhibitory neuron firing
                                                #      rates reaches the equilibrium.
StimTime  = 2;                                  # [ms] duration of an external stimulation
RepNum    = 3;                                  #      number of stimulations, i.e., number of trials
RepTime   = 3000;                               # [ms] duration of each trial
simtime   = RelaxTime+ RepTime*RepNum;          # [ms] total duration of the simulation
dt        = 0.1;                                # [ms] integration time step
Tau_W     = 150;                                # [ms] relaxation time of adaptation
T_bins    = 5;                                  # [ms] temporal extension of the binning for the calculation of the firing rate
Delta     = 0.25;                               #      variance of the distribution of couplings \sigma^2= float(J[i][j])*float(Delta))) 
tot       = 100;                                #      this parameter must remain fixed, it represents the number of intervals
                                                #      in which we divide the nest.simulate.


# Adaptation level (g) and excitation level (C_ext) are two key parameters we focus on to generate qualitatively different spontaneous and stimulus-evoked dynamics.
# These two parameters must be chosen in accordance with the diagram of figure 1A.
from IPython.display import Image
Image(filename='Diagramma.jpg')


g      = 6*5+2.5;      # adaptation
C_Ext  = 12*5+3160;    # excitation
Stim   = 4             # stimulus intensity


# general partameter and nest kernel parameters setting,
B =g/1000.;
tempo_i=time.time();
start=0.0  
origin=0.0                          

# set the kernel parameters
nest.SetKernelStatus({"local_num_threads":6})
nest.SetKernelStatus({"resolution": dt, "print_time": False,"overwrite_files": True})

# define variable
#NU  =np.zeros((2,1, int(round(float(simtime)/T_bins-1))));


# Building network
#############################------------------------------------------------------------------------
print("Building network")
#############################------------------------------------------------------------------------

startbuild   = time.time()  # initialize the calculation of the time used to simulate

# define variables
NeuronPop    =[]            # Neuronal population
NoisePop     =[]            # poissonian generator 
DetectorPop  =[]            # Detectors device
MultimeterPop=[]            # Multimeter device


# build neuronal excitatory population
app2 = nest.Create("aeif_psc_delta", int(round(6300)),params={"C_m":     1.0,
                                                                           "g_L":     1.0/float(20),
                                                                           "t_ref":   float(2),
                                                                           "E_L":     0.0,
                                                                           "V_reset": float(15),
                                                                           "V_m":     0.,
                                                                           "V_th":    float(20),
                                                                           "Delta_T": 0.,
                                                                           "a":       0.0,
                                                                           "b":       float(B),
                                                                           "tau_w":   float(Tau_W),
                                                                           "V_peak":  float(20)})
NeuronPop.append(app2)

# build neuronal inhibitory population
app2 = nest.Create("aeif_psc_delta", int(round(2580)),params={"C_m":     1.0,
                                                                           "g_L":     1.0/float(10),
                                                                           "t_ref":   float(1),
                                                                           "E_L":     0.0,
                                                                           "V_reset": float(15),
                                                                           "V_m":     0.,
                                                                           "V_th":    float(20),
                                                                           "Delta_T": 0.,
                                                                           "a":       0.0,
                                                                           "b":       float(0),
                                                                           "tau_w":   float(Tau_W),
                                                                           "V_peak":  float(20)})
NeuronPop.append(app2)


# define and initialize the poisson generators
app3 = nest.Create("poisson_generator",params = {"rate"  : float(0.5*C_Ext*0.75),#3#2.54608
                                                 'origin':0.,
                                                 'start' :0.})
NoisePop.append(app3)

# define and initialize the spike recorder 
app4 = nest.Create("spike_recorder",params    = {"start":0.})
DetectorPop.append(app4)

# define and initialize the poisson generators
app3= nest.Create("poisson_generator",params  = {"rate": float(0.5*733*0.75),
                                                 'origin':0.,
                                                 'start':0.})
NoisePop.append(app3)

# define and initialize the spike recorder 
app4 = nest.Create("spike_recorder",params    = {"start":0.})
DetectorPop.append(app4)

# define and initialize the multimeter 
app5 = nest.Create("multimeter", params       = {'interval': dt,'record_from': ['V_m','w']})                                                          
MultimeterPop.append(app5)


endbuild = time.time()


# Connecting devices
#############################------------------------------------------------------------------------
print("Connecting devices")
#############################------------------------------------------------------------------------

startconnect = time.time()

J_Ext=[0.476279,2.20634]; #the efficacy of the excitatory synapses of neurons external to the network with internal neurons

for i in range(0,2):
    
    #connect the poisson processes external with excitatory synapses with the neurons of the network
    nest.Connect(NoisePop[i], NeuronPop[i], syn_spec={'synapse_model': 'static_synapse', 
                                              'delay': dt,        
                                              'weight': float(J_Ext[i]*(1+Delta))}) 
    nest.Connect(NoisePop[i], NeuronPop[i], syn_spec={'synapse_model': 'static_synapse', 
                                              'delay': dt,        
                                              'weight': float(J_Ext[i]*(1-Delta)) }) 
    
    #we connect the neuronal populations of the network with the measuring instruments
    nest.Connect(NeuronPop[i], DetectorPop[i], syn_spec={"weight": 1.0, "delay": dt})


# We  define the intrapopulation (E-E, I-I) and interpopulations connections (E-I, I-E) efficacy, connectivity, trasmission delay and neuronal population quantity
J     =  [[1.95 ,2.20634 ],[-1.09887, -1.09887]];             # the efficacy of the excitatory and inhibitory synapses 
c     =  [[0.00634921 ,0.00206349],[0.0498809, 0.0174867]];   # neuronal connectivity
Delay =  [[80 ,80 ],[20, 20 ]];                               # neuronal connectivity transmission delay
N     =  [6300,2580];                                         # neuronal pool extencion

#connect the neuronal pool
for i in range(0,2):
    for j in range(0,2):
        
        print("Connecting devices ",i,j,int(round(c[i][j]*N[int(i)])),float(1./float(2.99573227355/(float(Delay[i][j])-float(0.1)))),float(Delay[i][j]),float(J[i][j]),math.fabs(float(J[i][j])*float(Delta)))
        
        nest.Connect(NeuronPop[int(i)], NeuronPop[int(j)], {'rule': 'fixed_indegree', 'indegree':int(round(c[i][j]*N[int(i)]))},                                                
                                                  syn_spec={'synapse_model': 'static_synapse_hpc', 
                                                                    'delay':  nest.math.redraw(nest.random.exponential(beta=float(1./float(2.99573227355/(float(Delay[i][j])-float(0.1))))),
                                                                              min= float(0.1)-dt/2,
                                                                              max= float(Delay[i][j])),
                                                                    'weight': nest.random.normal(mean=float(J[i][j]),std=math.fabs(float(J[i][j])*float(Delta))) })
        print("Connecting devices ",i,j)

endconnect = time.time()


# Simulating.
# Takes time (minutes)

#############################------------------------------------------------------------------------
print("Simulating")
#############################------------------------------------------------------------------------
sys.stdout.flush()


# define temporary variables
Nu_temp_E=[];
Nu_temp_I=[];

# split and execute the RelaxTime simulation 
for ii in range(0,tot):
    nest.Simulate(np.round(RelaxTime/tot,1))
print("Relaxing")
sys.stdout.flush()


# download the data from the detectors and free the memory
for k in range(0,len(DetectorPop),1):
    y=nest.GetStatus(DetectorPop[k], "events")[0]['times'];
    if k==0:
        Nu_temp_E.append(y)
    if k==1:
        Nu_temp_I.append(y)
    nest.SetStatus(DetectorPop[k],{'n_events':0})


# split and execute the stimulation trial simulation 
for x in range(0,RepNum):
    nest.SetStatus(NoisePop[0],{'rate':float(0.5*C_Ext*Stim)})
    nest.Simulate(StimTime)
    nest.SetStatus(NoisePop[0],{'rate':float(0.5*C_Ext*0.75)})
    appo=0;
    for ii in range(0,tot):
        nest.Simulate(np.round((RepTime-StimTime)/tot,1))
        appo=appo+np.round((RepTime-StimTime)/tot,1);
    nest.Simulate(np.round((appo-StimTime)/tot,1))
    sys.stdout.flush()
    
    # download the data from the detectors and free the memory
    for k in range(0,len(DetectorPop),1):
        y=nest.GetStatus(DetectorPop[k], "events")[0]['times'];
        if k==0:
            Nu_temp_E.append(y)
        if k==1:
            Nu_temp_I.append(y)
        nest.SetStatus(DetectorPop[k],{'n_events':0})            


# data acquisition and analysis
T1=0;
T2=int(simtime)
n_bins=int(round((T2-T1)/T_bins));
Nu=np.zeros((len(DetectorPop),n_bins-1))


for k in range(0,len(DetectorPop),1):
    if k==0:
        Nu_temp_e = [item for sublist in Nu_temp_E for item in sublist]
        x=Nu_temp_e
    if k==1:
        Nu_temp_i = [item for sublist in Nu_temp_I for item in sublist]
        x=Nu_temp_i  

    x=np.sort(x)
    x=x[x>=T1 ];
    x=x[x<= T2 ];
    his,bins=np.histogram(x,bins=range(T1,T2,int(np.round(float(T2-T1)/n_bins)) ))
    if k==0:
        Nu[k]=his*1000./(N[int(k)]*float(T2-T1)/n_bins)
        V=nest.GetStatus(MultimeterPop[k])[0]['events']['V_m']
        W=nest.GetStatus(MultimeterPop[k])[0]['events']['w']
        Time=nest.GetStatus(MultimeterPop[k])[0]['events']['times']  
        MultiTarget=nest.GetStatus(MultimeterPop[k])[0]['events']['senders']  
        SpikeTime=nest.GetStatus(DetectorPop[k], "events")[0]['times']
        SpikeSender=nest.GetStatus(DetectorPop[k], "events")[0]['senders']
    elif k==1:
        Nu[k]=his*1000./(N[int(k)]*float(T2-T1)/n_bins)
    elif k==2:
        Nu[k]=his*1000./(N[int(k)]*float(T2-T1)/n_bins)

tempo_f=time.time();
print("completed in:",tempo_f-tempo_i)


# data plot
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'widget')
plt.figure()
plt.plot(bins[0:-1],Nu[0],label='ExcFg g='+str(B));
plt.plot(bins[0:-1],Nu[1],label='Inh g='+str(0));

plt.xlabel('t [ms]')
plt.ylabel('Firing Rate \nu [Hz]')
plt.title('Tau_w='+str(Tau_W)+' [mS]C_ext='+str(C_Ext))
plt.legend( loc='upper left')
plt.show()
