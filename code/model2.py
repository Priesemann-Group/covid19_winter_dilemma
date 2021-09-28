from scipy.integrate import solve_ivp
from scipy.stats import gamma as gamma_dist
import numpy as np


class Model:
    def __init__(
        self,
        y0,
        Rt_base, Rt_free,
        eta, kappa, sigma,
        gamma, gamma_ICU,
        delta,
        rho,
        omega_v_b, omega_n_b,
        chi_0, chi_1,
        alpha_w, alpha_u, alpha_R,
        e_R, e_u, e_w,
        Phi_0,phi_0,
        u_base,
        mu, d_0, d_mu,
        a_Rt, a_vac, gamma_cutoff,
        tau_vac1, tau_vac2,
        t_max, step_size,
        Theta, Theta_ICU, influx
    ):
        self.y0 = y0
        self.Rt_base = Rt_base
        self.Rt_free = Rt_free
        self.eta = eta
        self.kappa = kappa
        self.sigma = sigma
        self.gamma = gamma
        self.gamma_ICU = gamma_ICU
        self.delta = delta
        self.rho = rho
        self.omega_v_b = omega_v_b
        self.omega_n_b = omega_n_b
        self.chi_0 = chi_0
        self.chi_1 = chi_1
        self.alpha_w = alpha_w
        self.alpha_u = alpha_u
        self.alpha_R = alpha_R
        self.e_R = e_R
        self.e_u = e_u
        self.e_w = e_w
        self.Phi_0 = Phi_0
        self.phi_0 = phi_0
        self.u_base = u_base
        self.mu = mu
        self.d_0 = d_0
        self.d_mu = d_mu
        self.a_Rt = a_Rt
        self.a_vac = a_vac
        self.gamma_cutoff = gamma_cutoff
        self.tau_vac1 = tau_vac1
        self.tau_vac2 = tau_vac2
        self.Theta = Theta
        self.Theta_ICU = Theta_ICU
        self.t_max = t_max
        self.step_size = step_size
        self.influx = influx

        self.M = sum(self.y0[:-4])
        self.u_max = 1-chi_0
        self.w_max = 1-chi_1
        self.c_v = 2. / ( self.M*self.omega_v_b*(1/(self.gamma+self.delta)) )
        self.c_n = 2. / ( self.M*self.omega_n_b*(1/(self.gamma+self.delta)) )
        self.t_min = self.gamma_cutoff + max(self.tau_vac1,self.tau_vac2)

        self.n = 1


    def time2index(self, t):
        return round((t+self.t_min)/self.step_size)

    def H_Rt(self, t):
        x = -np.arange(-self.gamma_cutoff,0,self.step_size)
        d = self.data[self.time2index(t-self.gamma_cutoff):self.time2index(t),7]
        g = gamma_dist.pdf(x,self.a_Rt)
        return (d*g).sum()*self.step_size

    def H_vac1(self, t):
        x = -np.arange(-self.gamma_cutoff,0,self.step_size)
        d = self.data[self.time2index(t-self.gamma_cutoff-self.tau_vac1):self.time2index(t-self.tau_vac1),7]
        g = gamma_dist.pdf(x,self.a_vac)
        return (d*g).sum()*self.step_size

    def H_vac2(self, t):
        x = -np.arange(-self.gamma_cutoff,0,self.step_size)
        d = self.data[self.time2index(t-self.gamma_cutoff-self.tau_vac2):self.time2index(t-self.tau_vac2),7]
        g = gamma_dist.pdf(x,self.a_vac)
        return (d*g).sum()*self.step_size

    def I_eff(self, I, IB):
        return I + self.sigma*IB + self.influx

    def Gamma(self, t):
        return 1 + self.mu*np.cos(2*np.pi*(t+self.d_0-self.d_mu)/360.)

    def R_0(self, t):
        return self.Rt_base if t<180 else self.Rt_free

    def Rt(self, t):
        return self.R_0(t)*np.exp(-self.alpha_R*self.H_Rt(t)-self.e_R) * self.Gamma(t) /self.Gamma(360-self.d_0)

    def Phi(self, t):
        return 4*self.Phi_0/(1-self.chi_0)*(1-np.exp(-self.alpha_u*self.H_vac1(t)))

    def phi(self, t):
        return 2*4*self.phi_0/(1-self.chi_1)*(1-np.exp(-self.alpha_w*self.H_vac2(t)))

    def omega_v(self, t, I, IB):
        return 2*self.omega_v_b*(1-1/(1+np.exp(-self.c_v*self.Rt(t)*self.I_eff(I,IB))))

    def omega_n(self, t, I, IB):
        return 2*self.omega_n_b*(1-1/(1+np.exp(-self.c_v*self.Rt(t)*self.I_eff(I,IB))))

    def fun(self, t, y):
        (S,V,W,E,EB,I,IB,ICU,R,UC,WC,D,C) = y
        
        dUC = self.Phi(t)*UC*(1-UC/(self.M*(1-self.chi_0)))
        dWC = self.phi(t)*WC*(1-WC/(V*(1-self.chi_1)))
        
        dS = -self.gamma*self.Rt(t)*S/self.M*self.I_eff(I,IB) - dUC
        dV = -(1-self.eta)*self.gamma*self.Rt(t)*V/self.M*self.I_eff(I,IB) + dUC + dWC - self.omega_v(t,I,IB)*V
        dW = self.omega_v(t,I,IB)*V - self.gamma*self.Rt(t)*W/self.M*self.I_eff(I,IB) - dWC + self.omega_n(t,I,IB)*R
        dE = self.gamma*self.Rt(t)*(S)/self.M*self.I_eff(I,IB) - self.rho*E
        dEB = (1-self.eta)*self.gamma*self.Rt(t)*V/self.M*self.I_eff(I,IB) + self.gamma*self.Rt(t)*W/self.M*self.I_eff(I,IB)- self.rho*EB
        dI = self.rho*E - (self.gamma+self.delta+self.Theta)*I
        dIB = self.rho*EB - (self.gamma + (self.Theta+self.delta)*(1-self.kappa))*IB
        dICU = self.delta*(I+(1-self.kappa)*IB) - (self.Theta_ICU+self.gamma_ICU)*ICU
        dR = self.gamma*(I+IB) - self.omega_n(t,I,IB)*R + self.gamma_ICU*ICU
        dD = self.Theta*I+(1-self.kappa)*self.Theta*IB+self.Theta_ICU*ICU
        dC = (S+W+(1-self.eta)*V)*self.gamma*self.Rt(t)*self.I_eff(I,IB)/self.M
    
        return [dS,dV,dW,dE,dEB,dI,dIB,dICU,dR,dUC,dWC,dD,dC]


    def run(self):
        times = np.arange(0,self.t_max,self.step_size)
        self.data = np.zeros((self.time2index(self.t_max)+100,len(self.y0)))    # solve_ivp tends to look in the future, H_* needs values
        self.data[:self.time2index(0)+1,:] = [self.y0]*(self.time2index(0)+1)

        for i in range(len(times)-1):
            res = solve_ivp(self.fun, (times[i],times[i+1]), self.data[self.time2index(i*self.step_size)])
            self.data[self.time2index(i*self.step_size)+1,:] = res["y"][:,-1:].reshape(len(self.y0))

        self.times = times

        return self.times, self.chopped_data()

    def chopped_data(self):
        return self.data[self.time2index(0):-100,:]
    