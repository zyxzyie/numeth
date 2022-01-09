from numpy import log
import matplotlib.pyplot as plt

# basic parameter
class para:
    
    # parameter
    def __init__(self, a, b, N):
        
        self.a = a
        self.b = b
        self.N = N
        self.h = (b-a) / N
        
# Define Euler's method
class euler(para):
    
    # parameter
    def __init__(self, a, b, N, w):
        
        para.__init__(self, a, b, N)
        self.w = w
        
        self.t = a
        self.y = w
        
        self.lst_t = []
        self.lst_w = []
        self.lst_y = []

    # Plot figure
    def Plot(self):
        
        plt.plot(self.lst_t, self.lst_y, "r", linewidth=2.5, label="Exact Solution")
        plt.plot(self.lst_t, self.lst_w, "bo-", linewidth=2, label="Numerical Solution")
        plt.legend(loc="best", shadow=True)
        plt.grid(1)
        
    # ODE
    def f(self):
        
        ODE = ((self.w/self.t) - (self.w/self.t)**2)
        return (ODE)
    
    # real solution
    def y_exact(self):
        
        func = (self.t / (1 + log(self.t)))
        return (func)
    
    # recursive formula for Euler
    def formula_Euler(self):
        
        recu = self.f()
        return (recu)
    
    # Euler's method
    def Method(self, show=True, figPlot=True, figSave=False):
        
        if show == True:
            print("t_i %3s w_i %7s y_i %7s Error" %("", "", ""))
        
        for i in range(0, self.N+1):
            self.lst_t.append(self.t)
            self.lst_w.append(self.w)
            self.lst_y.append(self.y_exact())
            self.w = self.w + self.h * self.formula_Euler()
            self.t = self.a + i * self.h
            if show == True:
                print("%.1f   %.7f   %.7f   %.7f"
                      %(self.t, self.w, self.y_exact(), abs(self.y_exact()-self.w)))
            
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("Euler Method")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)