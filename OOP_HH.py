import numpy as np
import matplotlib.pyplot as plt


# Инкапсуляция
# Полиформизм
# Наследование

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


class Channel:
    def __init__(self, args):
        self.E = args["E"]
        self.gmax = args["gmax"]
        self.degrees = args["degrees"]

        V = np.asarray([-60.0, ])
        self.x, _ = self.get_x_inftau(V)

    def get_x_inftau(self):
        return 0.0, 0.0

    def update(self, dt, V):
        x_inf, x_tau = self.get_x_inftau(V)
        self.x = x_inf - (x_inf - self.x) * exp(-dt / x_tau)
        I = self.gmax * np.prod(self.x**self.degrees) * (self.E - V)
        return I

class Kdr_Channel(Channel):
    def __init__(self, args):
        super(Kdr_Channel, self).__init__(args)

    def get_x_inftau(self, V):
        alpha = alpha_n(V)
        beta = beta_n(V)

        n_tau = 1.0 / (alpha + beta)
        n_inf = alpha * n_tau
        return n_inf, n_tau

class Na_Channel(Channel):
    def __init__(self, args):
        super(Na_Channel, self).__init__(args)

    def get_x_inftau(self, V):
        alpham = alpha_m(V)
        betam = beta_m(V)
        m_tau = 1.0 / (alpham + betam)
        m_inf = alpham * m_tau

        alphah = alpha_h(V)
        betah = beta_h(V)
        h_tau = 1.0 / (alphah + betah)
        h_inf = alphah * h_tau

        x_inf = np.append(m_inf, h_inf)
        x_tau = np.append(m_tau, h_tau)
        return x_inf, x_tau


class Neuron:
    def __init__(self, args):
        self.Cm = args['Cm']
        self.EL = args['EL']

        self.gL = args['gL']

        self.Iapp = args['Iapp']
        self.V = np.array( [-60.0, ])

        na_ch = Na_Channel({'gmax': args['gNa'], 'E':args['ENa'], 'degrees':np.asarray([3, 1])})
        kdr_ch = Kdr_Channel({'gmax': args['gK'], 'E':args['EK'], 'degrees':np.asarray([4, ])})
        self.channels = [na_ch, kdr_ch]





    def update(self, delta_t, duration):
        N = int(duration / delta_t)
        Vhist = [self.V, ]

        for i in range(N):
            I = 0
            for ch in self.channels:
                I += ch.update(delta_t, self.V)

            self.V = self.V + delta_t * (self.gL * (self.EL - self.V) + I + self.Iapp) / self.Cm
            Vhist.append(self.V)

        return Vhist
################################################################
params = {
    'Cm' : 1,
    'EL' : -65.0,
    'EK' : -90.0,
    'ENa' : 55.0,
    'gL' : 0.1,
    'gK' : 9.0,
    'gNa' : 35.0,
    'Iapp' : 0.5,
}
neuron = Neuron(params)

Vhist = neuron.update(0.1, 200)

fig, axes = plt.subplots()
axes.plot(Vhist)
axes.set_title('V')
plt.show()

