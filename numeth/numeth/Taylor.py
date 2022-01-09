from numpy import exp
import matplotlib.pyplot as plt

# basic parameter
class para:
    
    # parameter
    def __init__(self, a, b, N):
        
        self.a = a
        self.b = b
        self.N = N
        self.h = (b-a) / N

# Define Taylor's method
class taylor(para):
    
    # parameter
    def __init__(self, a, b, N, w):
        
        para.__init__(self, a, b, N)
        self.t = a
        self.w = w
        self.y = w
        self.lst_t = []
        self.lst_w = []
        self.lst_y = []
    
    # plot figure
    def Plot(self):
        
        plt.plot(self.lst_t, self.lst_y, "r", linewidth=2.5, label="Exact Solution")
        plt.plot(self.lst_t, self.lst_w, "bo-", linewidth=2, label="Numerical Solution")
        plt.legend(loc="best", shadow=True)
        plt.grid(1)
    
    # ODE
    def f(self):
        
        ODE = (2/self.t) * self.w + self.t**2 * exp(self.t)
        return (ODE)
    
    # real solution
    def y_exact(self):
        
        func = self.t**2*(exp(self.t) - exp(1))
        return (func)
    
    # recursive formula for Taylor's of order 2
    def formula_T2(self):

        recu = 2 / self.t*self.w + self.t**2*exp(self.t) + self.h / 2 * \
                 (self.t**2*exp(self.t) + 4*self.t*exp(self.t) + 2/self.t**2*self.w)
        return (recu)
    
    # Taylor's method of order 2
    def Tay_deg2(self, show, figPlot, figSave):

        for i in range(0, self.N+1):
            self.lst_t.append(self.t)
            self.lst_w.append(self.w)
            self.lst_y.append(self.y_exact())
            self.w = self.w + self.h * self.formula_T2()
            self.t = self.a + i * self.h
            if show == True:
                print("%.1f   %.7f   %.7f   %.7f"
                      %(self.t, self.w, self.y_exact(), abs(self.y_exact()-self.w)))
            
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("Taylor's Method of Order 2")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)
    
    # recursive formula for Taylor's of order 4
    def formula_T4(self):

        recu = 2/self.t*self.w + self.t**2*exp(self.t) +\
               self.h / 2*(2*self.t*exp(self.t) + self.t**2*exp(self.t) -\
               2/self.t**2*self.w + 2/self.t*(2/self.t*self.w+self.t**2*exp(self.t)))+ \
               self.h**2 / 6*(self.t**2*exp(self.t) + 6*self.t*exp(self.t) + 6*exp(self.t))+ \
               self.h**3 / 24*(self.t**2*exp(self.t) + 8*self.t*exp(self.t) + 12*exp(self.t))
        return (recu)
    
    # Taylor's method of order 4
    def Tay_deg4(self, show, figPlot, figSave):
            
        for i in range(1, self.N+2):
            if show == True:
                print("%.1f   %.7f   %.7f   %.7f"
                      %(self.t, self.w, self.y_exact(), abs(self.y_exact()-self.w)))
            self.lst_t.append(self.t)
            self.lst_w.append(self.w)
            self.lst_y.append(self.y_exact())
            self.w = self.w + self.h * self.formula_T4()
            self.t = self.a + i * self.h
            
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("Taylor's Method of Order 4")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)
            
    # Taylor's method
    def Method(self, deg, show=True, figPlot=True, figSave=False):
        
        if show == True:
            print("t_i %3s w_i %7s y_i %7s Error" %("", "", ""))
        
        if deg == 2:
            self.Tay_deg2(show, figPlot, figSave)
        
        elif deg == 4:
            self.Tay_deg4(show, figPlot, figSave)