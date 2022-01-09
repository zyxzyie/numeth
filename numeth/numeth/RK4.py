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
        plt.plot(self.lst_t, self.lst_w, "bo-", linewidth=2, label="Numerical Solution")
        plt.legend(loc="best", shadow=True)
        plt.grid(1)

# Define Runge-Kutta method of order 4
class RK4(para):
    
    # parameter
    def __init__(self, a, b, N, w):
        
        para.__init__(self, a, b, N)
        self.w = w
        self.t = a
        self.y = w
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
    
    # recursive formula for RK4
    def formula_RK4(self):
        
        k1 = (self.h*self.f(self.t, self.w))
        k2 = (self.h*self.f(self.t+self.h/2, self.w+k1/2))
        k3 = (self.h*self.f(self.t+self.h/2, self.w+k2/2))
        k4 = (self.h*self.f(self.t+self.h,self.w+k3))
        
        recu = (k1 + 2*k2 + 2*k3 + k4) / 6
        return (recu)
    
    # modified Euler's method
    def Method(self, show=True, figPlot=False, figSave=False):
        
        if show == True:
            print("t_i %3s w_i %7s y_i %7s Error" %("", "", ""))
        
        for i in range(1, self.N+2):
            if show == True:
                print("%.1f   %.7f   %.7f   %.7f"
                      %(self.t, self.w, self.y_exact(), abs(self.y_exact()-self.w)))
            self.lst_t.append(self.t)
            self.lst_w.append(self.w)
            self.lst_y.append(self.y_exact())
            self.w = self.w + self.formula_RK4()
            self.t = self.a + i * self.h
            
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("RK4")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)