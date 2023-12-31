import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def func(t, Y, sigma, rho, beta):
    x, y, z = Y
    dxdt = sigma*(y - x)
    dydt = x*(rho - z) - y
    dzdt = x*y - beta*z

    dYdt = [dxdt, dydt, dzdt]
    return dYdt
#############################

class Diff_Eqs:
    def __init__(self, sigma, rho, beta):
        self.sigma = sigma
        self.rho = rho
        self.beta = beta

    def func(self, t, Y):
        x, y, z = Y
        dxdt = self.sigma * (y - x)
        dydt = x * (self.rho - z) - y
        dzdt = x * y - self.beta * z

        dYdt = [dxdt, dydt, dzdt]
        return dYdt
##############################

sigma = 10
rho = 28
beta = 8/3
y0 = [2.0, 1.0, 1.0]
t_span = [0, 100]
t_eval = np.linspace(0, 100, 1000)
#sol = solve_ivp(func, t_span, y0, method='LSODA', args=(sigma, rho, beta), t_eval=t_eval, rtol=10e-6)

system = Diff_Eqs(sigma, rho, beta)
sol = solve_ivp(system.func, t_span, y0, method='LSODA', t_eval=t_eval, rtol=10e-6)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.plot(sol.y[0], sol.y[1], sol.y[2])
plt.show()
