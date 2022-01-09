import matplotlib.pyplot as plt

# basic parameter
class para:
    
    def __init__(self, a, b, N):
        
        self.a = a
        self.b = b
        self.N = N
        self.h = (b-a) / N

# Define modified Euler's method
class mod_euler(para):
    
    # parameter
    def __init__(self, a, b, N, w):
        
        para.__init__(self, a, b, N)
        self.w = w
        self.t = a
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
    def f(self, t, w):
        
        ODE = (2-2*t*w) / (t**2+1)
        return (ODE)
    
    # real solution
    def y_exact(self):
        
        func = (2*self.t+1) / (self.t**2+1)
        return (func)
    
    # recursive formula for midpoint
    def formula_modEuler(self):

        recu = self.f(self.t, self.w) +\
               self.f(self.t+self.h, self.w+self.h*self.f(self.t, self.w))
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
            self.w = self.w + self.h/2*(self.formula_modEuler())
            self.t = self.a + i * self.h
            
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("Midpoint Method")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)