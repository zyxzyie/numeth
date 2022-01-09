from numpy import exp
import matplotlib.pyplot as plt

# basic parameter
class para:
    
    def __init__(self, a, b, N):
        
        self.a = a
        self.b = b
        self.N = N
        self.h = (b-a) / N
    
    def Plot(self):
        
        plt.plot(self.lst_t, self.lst_y, "r", linewidth=2.5, label="Exact Solution")
        plt.plot(self.lst_t, self.lst_w[:6], "bo-", linewidth=2, label="Numerical Solution")
        plt.legend(loc="best", shadow=True)
        plt.grid(1)

# Define Adam-Fourth order predictor corrector method
class AF_method(para):
    
    # parameter
    def __init__(self, a, b, N, w):
        
        para.__init__(self, a, b, N)
        self.w = w
        self.t = a
        self.y = w
        self.w0 = w
        self.w1 = 0
        self.w2 = 0
        self.w3 = 0
        self.w_pred = 0
        self.w_core = 0
        self.lst_t = []
        self.lst_w = []
        self.lst_y = []
    
    # ODE
    def f(self, t, w):
        
        ODE = t*exp(3*t) - 2*w
        return (ODE)
    
    # real solution
    def y_exact(self):
        
        func = 1/5*self.t*exp(3*self.t) - 1/25*exp(3*self.t) + 1/25*exp(-2*self.t)
        return (func)
    
    # recursive formula of predictor for AF
    def predictor_AF(self):
        
        recu = (55*self.f(self.t+3*self.h, self.w3)-
                59*self.f(self.t+2*self.h, self.w2)+
                37*self.f(self.t+self.h, self.w1)-
                9*self.f(self.t, self.w0))
        return (recu)
    
    # recursive formula of corrector for AF
    def corrector_AF(self):
        
        recu = (9*self.f(self.t+4*self.h, self.w_pred)+
                19*self.f(self.t+3*self.h, self.w3)-
                5*self.f(self.t+2*self.h, self.w2)+
                self.f(self.t+self.h, self.w1))
        return (recu)
    
    def RK4(self, t, w, h):
        k1 = h * self.f(t, w)
        k2 = h * self.f(t+h/2, w+k1/2)
        k3 = h * self.f(t+h/2, w+k2/2)
        k4 = h * self.f(t+h, w+k3)
        recu = (k1 + 2*k2 + 2*k3 + k4)/6
        return (recu)
    
    # AF method
    def Method(self, show=True, figPlot=False, figSave=False):
        
        if show == True:
            print("t_i %3s w_i %7s y_i %7s Error" %("", "", ""))
        
        self.w1 = self.w0 + self.RK4(self.t, self.w0, self.h)
        self.w2 = self.w1 + self.RK4(self.t+self.h, self.w1, self.h)
        self.w3 = self.w2 + self.RK4(self.t+2*self.h, self.w2, self.h)
        self.lst_w = [self.w0, self.w1, self.w2, self.w3]
        for i in range(1, self.N+2):
            self.lst_t.append(self.t)
            self.lst_y.append(self.y_exact())
            x1 = self.w1
            x2 = self.w2
            x3 = self.w3
            self.w_pred = self.w3 + self.h/24*self.predictor_AF()
            self.w_core = self.w3 + self.h/24*self.corrector_AF()
            self.lst_w.append(self.w_core)
            self.t = self.a + i * self.h
            self.w3 = self.w_core
            self.w2 = x3
            self.w1 = x2
            self.w0 = x1
    
        for i in range(0, 6):
            print("%.1f   %.7f   %.7f   %.7f" 
                   %(self.lst_t[i], self.lst_w[i],
                     self.lst_y[i], abs(self.lst_w[i]-self.lst_y[i]))) 
        
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("AF Method")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)
