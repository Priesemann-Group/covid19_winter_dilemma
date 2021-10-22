from scipy.integrate import solve_ivp
from scipy.special import gamma as gamma_func
import numpy as np
import pickle


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
        vac_max,
        u_base,
        mu, d_0, d_mu,
        a_Rt, b_Rt, a_vac, b_vac, gamma_cutoff,
        tau_vac1, tau_vac2,
        t_max, step_size,
        Theta, Theta_ICU,
        influx,
        time_u, time_w,
        epsilon_free, CM,
        w_ICU, fractions, plateaus, slopes
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
        self.vac_max = vac_max
        self.u_base = u_base
        self.mu = mu
        self.d_0 = d_0
        self.d_mu = d_mu
        self.a_Rt = a_Rt
        self.b_Rt = b_Rt
        self.a_vac = a_vac
        self.b_vac = b_vac
        self.gamma_cutoff = gamma_cutoff
        self.tau_vac1 = tau_vac1
        self.tau_vac2 = tau_vac2
        self.Theta = Theta
        self.Theta_ICU = Theta_ICU
        self.t_max = t_max
        self.step_size = step_size
        self.influx = influx
        self.time_u = time_u
        self.time_w = time_w
        self.epsilon_free = epsilon_free
        self.CM = CM
        self.w_ICU = w_ICU
        self.fractions = fractions
        self.plateaus = plateaus
        self.slopes = slopes


        self.eqs = 18           # number equations
        self.ags = len(y0)//self.eqs  # number agegroups

        self.M = self.y0.reshape([self.eqs,self.ags])[:-4,:].sum(axis=0)
        self.u_max = 1-chi_0
        self.w_max = 1-chi_1
        self.c_v = 2. / ( self.M*self.omega_v_b*(1/(self.gamma+self.delta+self.Theta)) )
        self.c_n = 2. / ( self.M*self.omega_n_b*(1/(self.gamma+self.delta+self.Theta)) )
        self.t_min = self.gamma_cutoff + max(self.tau_vac1,self.tau_vac2)

        self.n = 1


    def time2index(self, t):
        return round((t+self.t_min)/self.step_size)

    def H(self,t,tau,a,b):
        x = -np.arange(-self.gamma_cutoff,0,self.step_size)
        d = self.data[self.time2index(t-self.gamma_cutoff-tau):self.time2index(t-tau),10:12].sum(axis=1)
        d = (1-self.w_ICU) * d.sum(axis=1) / (self.M.sum()/1e6) + (self.w_ICU * d / (self.M/1e6)).transpose()
        g = b**a * x**(a-1) * np.e**(-b*x)/gamma_func(a)
        return (d*g).sum(axis=1)*self.step_size

    def H_Rt(self, t):
        return self.H(t, 0, self.a_Rt, self.b_Rt)

    def H_vac1(self, t):
        return self.H(t, self.tau_vac1, self.a_vac, self.b_vac)

    def H_vac2(self, t):
        return self.H(t, self.tau_vac2, self.a_vac, self.b_vac)

    def I_eff(self, I, IBn, IBv):
        return (I + self.sigma*(IBn+IBv)) + self.influx*self.M/self.M.sum()

    def Gamma(self, t):
        return 1 + self.mu*np.cos(2*np.pi*(t+self.d_0-self.d_mu)/360.)

    def R_0(self, t):
        if t<180-self.epsilon_free:
            return self.Rt_base
        if t>180+self.epsilon_free:
            return self.Rt_free
        else:
            return self.Rt_base + (self.Rt_free-self.Rt_base)*(t-180+self.epsilon_free)/2/self.epsilon_free

    #Reproduction Number Old    
    #def Rt(self, t):
     #   return self.R_0(t)*np.exp(-self.alpha_R*self.H_Rt(t)-self.e_R) * self.Gamma(t) /self.Gamma(360-self.d_0)
    
    #Reproduction number new 
    
    #Cosmo Data fit 
    def selfregulation(self, t):
        return self.plateaus - self.slopes*1*np.log(np.e**(1/1*(35-self.H_Rt(t)))+1)
        
    # We subtract 2.5 and divide by 2.5 such that we project (2.5,5) -> (0,1) 
    def Rt(self,t):
        return self.R_0(t) * (1-self.fractions*(self.selfregulation(t)-2.5)/2.5) * self.Gamma(t) /self.Gamma(360-self.d_0)
            
    def u_w(self, t):
        return self.u_base + (self.u_max-self.u_base)*(1-np.exp(-self.alpha_u*self.H_vac1(t)-self.e_u))

    def Phi(self, t, UC, frac):
        ppl = (self.u_w(t))*self.M - UC
        return ( ppl > 0 ) * ppl * frac / self.time_u

    def w_w(self, t):
        return self.w_max*(1-np.exp(-self.alpha_w*self.H_vac2(t)-self.e_w))

    def phi(self, t, UC, WC, frac):
        ppl = ((self.w_w(t))*UC - WC)
        return ( ppl > 0 ) * ppl * frac / self.time_w

    def get_phis(self, t, y):   # for plotting only right now; TODO later: check speed and use if not to slow
        (S,V,Wn,Wv,E,EBn,EBv,I,IBn,IBv,ICU,ICUv,R,Rv,UC,WC,D,C) = np.split(y, self.eqs)
        Phi = self.Phi(t, UC, (S+Wn)/(self.M-UC))
        phi = self.phi(t, UC, WC, (Wv)/(UC-WC))
#        ratio = (Phi+phi).sum() / (self.vac_max*self.M.sum())
#        if ratio > 1:
#            Phi /= ratio
#            phi /= ratio
        return Phi, phi

    def omega_v(self, t, I, IBn, IBv):
        return self.omega_v_b*(1-1/(1+np.exp(-self.c_v*self.Rt(t)*self.I_eff(I,IBn,IBv))))*2

    def omega_n(self, t, I, IBn, IBv):
        return self.omega_n_b*(1-1/(1+np.exp(-self.c_v*self.Rt(t)*self.I_eff(I,IBn,IBv))))*2

    def fun(self, t, y):
        y.reshape([self.eqs,self.ags])
        (S,V,Wn,Wv,E,EBn,EBv,I,IBn,IBv,ICU,ICUv,R,Rv,UC,WC,D,C) = np.split(y, self.eqs)

        # definitions to make DEs more readable
        M = self.M
        eta = self.eta
        rho = self.rho
        kappa = self.kappa
        gamma = self.gamma
        gamma_ICU = self.gamma_ICU
        delta = self.delta
        Theta = self.Theta
        Theta_ICU = self.Theta_ICU
        Rt = self.Rt(t)
        I_eff = self.I_eff(I, IBn, IBv)
        omega_n = self.omega_n(t, I, IBn, IBv)
        omega_v = self.omega_v(t, I, IBn, IBv)
        Phi = self.Phi(t, UC, (S+Wn)/(M-UC))
        phi = self.phi(t, UC, WC, (Wv)/(UC-WC))
        CM = self.CM
        
        #infect = gamma*Rt*I_eff/M.sum()
        infect = gamma*Rt * np.matmul(CM.transpose(), I_eff/M )


        # differential equations
        dS = -S*infect - Phi*(S/(S+Wn))
        dV = -(1-eta)*V*infect + Phi + phi - omega_v*V
        dWn = -Wn*infect + omega_n*R - Phi*(Wn/(S+Wn))
        dWv = -Wv*infect + omega_v*V + omega_n*Rv - phi
        dE = S*infect - rho*E
        dEBn = Wn*infect - rho*EBn
        dEBv = ((1-eta)*V+Wv)*infect - rho*EBv
        dI = rho*E - (gamma+delta+Theta)*I
        dIBn = rho*EBn - (gamma + (Theta+delta)*(1-kappa))*IBn
        dIBv = rho*EBv - (gamma + (Theta+delta)*(1-kappa))*IBv
        dICU = delta*(I + (1-kappa)*IBn) - (Theta_ICU+gamma_ICU)*ICU
        dICUv = delta*(1-kappa)*IBv - (Theta_ICU+gamma_ICU)*ICUv
        dR = gamma*(I+IBn) - omega_n*R + gamma_ICU*ICU
        dRv = gamma*IBv - omega_n*Rv + gamma_ICU*ICUv
        dUC = Phi
        dWC = phi
        dD = Theta*I + (1-kappa)*Theta*(IBn+IBv) + Theta_ICU*(ICU+ICUv)
        dC = (S+Wn+Wv+(1-eta)*V)*infect

        return np.concatenate([dS,dV,dWn,dWv,dE,dEBn,dEBv,dI,dIBn,dIBv,dICU,dICUv,dR,dRv,dUC,dWC,dD,dC]).flatten()

    def build_data(self):
        times = np.arange(0,self.t_max,self.step_size)
        self.data = np.zeros((self.time2index(self.t_max)+100,self.eqs,self.ags))    # solve_ivp tends to look in the future, H_* needs values
        self.data[:self.time2index(0)+1,:,:] = [self.y0.reshape([self.eqs,self.ags])]*(self.time2index(0)+1)

    def run(self):
        times = np.arange(0,self.t_max,self.step_size)
        self.data = np.zeros((self.time2index(self.t_max)+100,self.eqs,self.ags))    # solve_ivp tends to look in the future, H_* needs values
        self.data[:self.time2index(0)+1,:,:] = [self.y0.reshape([self.eqs,self.ags])]*(self.time2index(0)+1)

        for i in range(len(times)-1):
            res = solve_ivp(self.fun, (times[i],times[i+1]), self.data[self.time2index(i*self.step_size)].flatten())
            self.data[self.time2index(i*self.step_size)+1,:,:] = res["y"][:,-1:].reshape(self.eqs,self.ags)

        self.times = times

        return self.times, self.chopped_data()

    def chopped_data(self):
        return self.data[self.time2index(0):-100,:,:]
    
    def save(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path):
        with open(path, "rb") as f:
            return pickle.load(f)

