import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
exp = np.exp

def alpha_m(V):
    alpha = 5.0 * -0.1*(V + 35) / (exp(-0.1*(V + 35)) - 1)
    return alpha

def beta_m(V):
    beta = 5.0 * 4 * exp(-(V + 60)/18)
    return beta


def alpha_h(V):
    alpha =  5.0 * 0.07 * exp(-(V + 58)/20)
    return alpha

def beta_h(V):
    beta = 5.0 * 1/(exp(-0.1*(V + 28)) + 1)
    return beta

def alpha_n(V):
    alpha = 5.0 * -0.01*(V + 34)/(exp(-0.1*(V + 34)) - 1)
    return alpha

def beta_n(V):
    beta = 5.0 * 0.125 * exp(-(V + 44)/80)
    return beta


Cm = 1

EL = -65.0
EK = -90.0
ENa = 55.0

gL = 0.1
gK = 9.0
gNa = 35.0

Iapp = 1.9


V = -60
m = alpha_m(V) / (alpha_m(V) + beta_m(V))
h = alpha_h(V) / (alpha_h(V) + beta_h(V))
n = alpha_n(V) / (alpha_n(V) + beta_n(V))


delta_t = 0.1

duration = 200
N = int(duration / delta_t)
Vhist = [V, ]
mhist = [m, ]
hhist = [h, ]
nhist = [n, ]


slope = -0.7729495609898452
bias = 0.6136305014913137

for i in range(N):
    # z_n = z_n-1 + delta t * f(z, t)
    INa = gNa * m**3 * h * (ENa - V)
    IK = gK * n**4 * (EK - V)
    V = V + delta_t * ( gL * (EL - V) + INa + IK + Iapp) / Cm


    m_tau = 1.0 / (alpha_m(V) + beta_m(V))
    h_tau = 1.0 / (alpha_h(V) + beta_h(V))
    #n_tau = 1.0 / (alpha_n(V) + beta_n(V))

    m_inf = alpha_m(V) * m_tau
    h_inf = alpha_h(V) * h_tau
    #n_inf = alpha_n(V) * n_tau

    m = m_inf # m + delta_t * ( (m_inf - m) / m_tau )
    #m = m_inf - (m_inf - m) * exp( -delta_t / m_tau)

    # h = h + delta_t * ( (h_inf - h) / h_tau )
    # n = n + delta_t * ( (n_inf - n) / n_tau )
    h = h_inf - (h_inf - h) * exp( -delta_t / h_tau)
    #n = n_inf - (n_inf - n) * exp( -delta_t / n_tau)

    n = h*slope + bias


    Vhist.append(V)
    mhist.append(m**3)
    hhist.append(h)
    nhist.append(n)

t = np.linspace(0, duration, N+1)

# hhist = np.asarray(hhist)
# nhist = np.asarray(nhist)
# res = linregress(hhist, nhist)
# print(res.slope, res.intercept)
#
# x = np.array([0.1, 0.7])
# y = x * res.slope + res.intercept
#
# fig, axes = plt.subplots(nrows=1, sharex=True)
# axes.scatter(nhist, hhist, label="n vs h")
# axes.plot(x, y, color='red')
# axes.set_xlabel("n")
# axes.set_ylabel("h")


fig, axes = plt.subplots(nrows=2, sharex=True)
axes[0].plot(t, Vhist)
axes[1].plot(t, mhist, label="m")
axes[1].plot(t, hhist, label="h")
axes[1].plot(t, nhist, label="n")
axes[1].legend()
plt.show()