import numpy as np
from brian2 import *

defaultclock.dt = 0.05*ms

Cm = 1*uF # /cm**2
Iapp = 2.0*uA
gL = 0.1*mS #msiemens
EL = -65*mV
ENa = 55*mV
EK = -90*mV
gNa = 35*msiemens
gK = 9*msiemens

sigma = 0.1*mV / sqrt(ms)


eqs = '''
dv/dt = (-gNa*m**3*h*(v-ENa)-gK*n**4*(v-EK)-gL*(v-EL)+Iapp + Isyn)/Cm + sigma*xi: volt
m = alpha_m/(alpha_m+beta_m) : 1
alpha_m = 0.1/mV*10*mV/exprel(-(v+35*mV)/(10*mV))/ms : Hz
beta_m = 4*exp(-(v+60*mV)/(18*mV))/ms : Hz
dh/dt = 5*(alpha_h*(1-h)-beta_h*h) : 1
alpha_h = 0.07*exp(-(v+58*mV)/(20*mV))/ms : Hz
beta_h = 1./(exp(-0.1/mV*(v+28*mV))+1)/ms : Hz
dn/dt = 5*(alpha_n*(1-n)-beta_n*n) : 1
alpha_n = 0.01/mV*10*mV/exprel(-(v+34*mV)/(10*mV))/ms : Hz
beta_n = 0.125*exp(-(v+44*mV)/(80*mV))/ms : Hz
Isyn : ampere
'''

neuron = NeuronGroup(100, eqs, method='heun', threshold='v > -10*mV')
neuron.v = '10*mV*randn() - 70*mV'        # -70*mV
neuron.h = 'alpha_h / (alpha_h + beta_h)' # 1
neuron.n = 'alpha_n / (alpha_n + beta_n)' # 1
M = StateMonitor(neuron, 'v', record=np.arange(10))

theta_syn = 0*mV
alpha_s = 12/ms
beta_s = 0.1/ms
syn_eqs = """
ds/dt = alpha_s * F * (1 - s) - beta_s * s : 1
F = 1/(1 + exp(-(v_pre - theta_syn)/(2*mV))) : 1
Isyn_post = gsyn*s*(Erev - v_post) : ampere (summed)
"""
# ds/dt = -s/(5*ms) : 1



gsyn = 0.1*mS
Erev = -75*mV
# on_pre = """
# s += 1
# """
synapses = Synapses(neuron, neuron,  model=syn_eqs,  method='rk4') #on_pre=on_pre
synapses.connect(p=0.6) # i=0, j=1

SpM = SpikeMonitor(neuron)

run(500*ms, report='text')

scatter(SpM.t/ms, SpM.i)


# for idx in range(2):
#     plot(M.t/ms, M[idx].v/mV, label=str(idx))
# legend(loc='upper right')
show()
