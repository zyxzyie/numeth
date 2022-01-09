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

# Define Adam-Bashforth method
class AB_method(para):
    
    # parameter
    def __init__(self, a, b, N, w):
        
        para.__init__(self, a, b, N)
        self.w = w
        self.t = a
        self.y = w
        self.w0 = w
        self.w1 = 0
        self.w2 = 0
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
    
    # recursive formula for AB2
    def formula_AB2(self):
        
        recu = (3*self.f(self.t+self.h, self.w1) - self.f(self.t, self.w0))
        return (recu)
    
    # AB2 method
    def AB2(self, show, figPlot, figSave):
        
        self.w1 = self.w0 + self.h * self.f(self.t, self.w0)
        self.lst_w = [self.w1]
        for i in range(1, self.N+2):
            self.lst_t.append(self.t)
            self.lst_w.append(self.w)
            self.lst_y.append(self.y_exact())
            x = self.w1
            self.w = self.w1 + self.h/2 * self.formula_AB2()
            self.t = self.a + i * self.h
            self.w1 = self.w
            self.w0 = x
    
        for i in range(0, 6):
            if show == True:
                print("%.1f   %.7f   %.7f   %.7f" 
                       %(self.lst_t[i], self.lst_w[i],
                         self.lst_y[i], abs(self.lst_w[i]-self.lst_y[i])))
        
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("AB2")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)
    
    # recursive formula for AB3
    def formula_AB3(self):
        
        recu = (23*self.f((self.t+2*self.h), self.w2)-
                16*self.f(self.t+self.h, self.w1)+5*self.f(self.t, self.w0))
        return (recu)
    
    # AB3 method
    def AB3(self, show, figPlot, figSave):
        
        self.w1 = self.w0 + self.h * self.f((self.t+self.h/2),
                                            (self.w0+self.h/2*self.f(self.t, self.w0)))
        self.w2 = self.w1 + self.h * self.f((0.2+self.h/2),
                                            (self.w1+self.h/2*self.f(0.2, self.w1)))
        self.lst_w = [self.w1, self.w2]
        for i in range(1, self.N+2):
            self.lst_t.append(self.t)
            self.lst_w.append(self.w)
            self.lst_y.append(self.y_exact())
            x1 = self.w1
            x2 = self.w2
            self.w = self.w2 + self.h/12 * self.formula_AB3()
            self.t = self.a + i * self.h
            self.w2 = self.w
            self.w1 = x2
            self.w0 = x1
    
        for i in range(0, 6):
            if show == True:
                print("%.1f   %.7f   %.7f   %.7f" 
                       %(self.lst_t[i], self.lst_w[i],
                         self.lst_y[i], abs(self.lst_w[i]-self.lst_y[i])))
        
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("AB3")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)
            
    # AB method
    def Method(self, deg, show=True, figPlot=False, figSave=False):
        
        if show == True:
            print("t_i %3s w_i %7s y_i %7s Error" %("", "", ""))
        
        if deg == 2:
            self.AB2(show, figPlot, figSave)
        
        elif deg == 3:
            self.AB3(show, figPlot, figSave)
test = AB_method(a=0, b=1, N=5, w=0).Method(deg=2, figPlot=1)