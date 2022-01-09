from numpy import exp
import matplotlib.pyplot as plt

# basic parameter
class para:
    
    # parameter
    def __init__(self, a, b):
        
        self.a = a
        self.b = b

# Define Runge-Kutta-Fehlberg method
class RKF(para):
    
    # parameter
    def __init__(self, a, b, w, TOL, h_max, h_min):
        
        para.__init__(self, a, b)
        self.t = a
        self.w = w
        self.y = w
        self.TOL = TOL
        self.h_max = h_max
        self.h_min = h_min
        self.h = h_max
        self.flag = 1
        self.R = 0
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
        
        ODE = w - t**2 + 1
        return (ODE)
    
    # real solution
    def y_exact(self):
        
        func = (self.t+1)**2 - 0.5*exp(self.t)
        return (func)
    
    # RKF method
    def Method(self, show=True, figPlot=True, figSave=False):
        
        if show == True:
            print("t_i %4s y_i %5s w_i %5s h %4s R %5s error" %("", "", "", "", ""))
        
        while self.flag == 1:
            k1 = self.h * self.f(self.t, self.w)
            k2 = self.h * self.f(self.t+1/4*self.h, self.w+1/4*k1)
            k3 = self.h * self.f(self.t+3/8*self.h, self.w+3/32*k1+9/32*k2)
            k4 = self.h * self.f(self.t+12/13*self.h, 
                                 self.w+1932/2197*k1-7200/2197*k2+7296/2197*k3)
            k5 = self.h * self.f(self.t+self.h, 
                                 self.w+439/216*k1-8*k2+3680/513*k3-845/4104*k4)
            k6 = self.h * self.f(self.t+1/2*self.h, 
                                 self.w-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)
        
            R = 1/self.h * abs(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)
        
            if R <= self.TOL:
                self.lst_t.append(self.t)
                self.lst_w.append(self.w)
                self.lst_y.append(self.y_exact())
                self.t = self.t + self.h
                self.w = self.w + 25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5
                if show == True:
                    print("%.4f %.7f %.7f %.4f %.7f %.7f"
                           %(self.t, self.y_exact(), self.w, self.h,
                             R, abs(self.y_exact()-self.w)))
                
        
            delta = 0.84*(self.TOL/R)**(1/4)
        
            if delta <= 0.1: self.h = 0.1*self.h
            elif delta >= 4: self.h = 4*self.h
            else: self.h = delta*self.h
            
            if self.h > self.h_max: self.h = self.h_max
            
            if self.t>= self.b: self.flag = 0
            elif (self.t + self.h) > self.b: self.h = self.b - self.t
            elif self.h < self.h_min: self.flag = 0
            
        if figPlot == True:
            plt.figure()
            self.Plot()
            plt.title("")
            plt.show()
            
        if figSave == True:
            plt.savefig("result.png", dpi=200)

# a, b, w = 0, 2, 0.5
# TOL = 10**(-5)
# h_max, h_min = 0.25, 0.01
# test = RKF(a, b, w, TOL, h_max, h_min).Method()