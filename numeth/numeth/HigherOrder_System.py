from numpy import exp
import pandas as pd

# basic parameter
class para:
    
    # parameter
    def __init__(self, a, b, N):
        
        self.a = a
        self.b = b
        self.N = N
        self.h = (b-a) / N

# Define the method of solving high order system
class highOrder_System(para):
    
    # parameter
    def __init__(self, a, b, N, w_u1, w_u2):
        
        para.__init__(self, a, b, N)
        self.t = a
        self.w_u1 = w_u1
        self.w_u2 = w_u2
        self.k1_1, self.k2_1, self.k3_1, self.k4_1 = [], [], [], []
        self.k1_2, self.k2_2, self.k3_2, self.k4_2 = [], [], [], []
        self.lst_t = []
        self.lst_u1, self.lst_u2 = [], []
        self.lst_w1, self.lst_w2 = [], []
        self.lst_error1, self.lst_error2 = [0], [0]
    
    # ODE of u1
    def u1_prime(self, t, u1, u2):
        return (3*u1 + 2*u2 - (2*t**2+1) * exp(2*t))
    
    # ODE of u2
    def u2_prime(self, t, u1, u2):
        return (4*u1 + u2 + (t**2+2*t-4) * exp(2*t))
    
    # real solution of u1
    def u1_exact(self):
        return (1/3*exp(5*self.t) - 1/3*exp(-self.t) + exp(2*self.t))
    
    # real solution of u2
    def u2_exact(self):
        return (1/3*exp(5*self.t) + 2/3*exp(-self.t) + self.t**2*exp(2*self.t))
    
    # the method of solving higher ordr system
    def Method(self, show=True):
        
        self.lst_t = [self.t]
        self.lst_u1, self.lst_u2 = [self.u1_exact()], [self.u2_exact()]
        self.lst_w1, self.lst_w2 = [self.w_u1], [self.w_u2]
        
        for i in range(1, self.N+1):
            for j in range(1, 3):
                if j == 1:
                    k1 = self.h * self.u1_prime(self.t, self.w_u1, self.w_u2)
                    self.k1_1.append(k1)
                else:
                    k1 = self.h * self.u2_prime(self.t, self.w_u1, self.w_u2)
                    self.k1_2.append(k1)
            for j in range(1, 3):
                if j == 1:
                    k2 = self.h * self.u1_prime(self.t+self.h/2, 
                                                self.w_u1+1/2*self.k1_1[-1], 
                                                self.w_u2+1/2*self.k1_2[-1])
                    self.k2_1.append(k2)
                else:
                    k2 = self.h * self.u2_prime(self.t+self.h/2, 
                                                self.w_u1+1/2*self.k1_1[-1], 
                                                self.w_u2+1/2*self.k1_2[-1])
                    self.k2_2.append(k2)
            for j in range(1, 3):
                if j == 1:
                    k3 = self.h * self.u1_prime(self.t+self.h/2, 
                                                self.w_u1+1/2*self.k2_1[-1], 
                                                self.w_u2+1/2*self.k2_2[-1])
                    self.k3_1.append(k3)
                else:
                    k3 = self.h * self.u2_prime(self.t+self.h/2, 
                                                self.w_u1+1/2*self.k2_1[-1], 
                                                self.w_u2+1/2*self.k2_2[-1])
                    self.k3_2.append(k3)
            for j in range(1, 3):
                if j == 1:
                    k4 = self.h * self.u1_prime(self.t+self.h, 
                                                self.w_u1+self.k3_1[-1], 
                                                self.w_u2+self.k3_2[-1])
                    self.k4_1.append(k4)
                else:
                    k4 = self.h * self.u2_prime(self.t+self.h, 
                                                self.w_u1+self.k3_1[-1], 
                                                self.w_u2+self.k3_2[-1])
                    self.k4_2.append(k4)
            for j in range(1, 3):
                if j == 1:
                    self.w_u1 = self.w_u1 + \
                                (self.k1_1[-1] + 2*self.k2_1[-1] + 
                                 2*self.k3_1[-1] + self.k4_1[-1])/6
                else:
                    self.w_u2 = self.w_u2 + \
                                (self.k1_2[-1] + 2*self.k2_2[-1] + 
                                 2*self.k3_2[-1] + self.k4_2[-1])/6
            self.t = self.a + i * self.h
            self.lst_t.append(self.t)
            self.lst_w1.append(self.w_u1)
            self.lst_w2.append(self.w_u2)
            self.lst_u1.append(self.u1_exact())
            self.lst_u2.append(self.u2_exact())
            self.lst_error1.append(abs(self.w_u1-self.u1_exact()))
            self.lst_error2.append(abs(self.w_u2-self.u2_exact()))
            
        RK4_system = pd.DataFrame({"t_i":self.lst_t,
                                   "w_1":self.lst_w1,
                                   "u_1":self.lst_u1,
                                   "Error1":self.lst_error1,
                                   "w_2":self.lst_w2,
                                   "u_2":self.lst_u2,
                                   "Error2":self.lst_error2})
        if show == True:
            print(RK4_system)
       
a, b, N, w_u1 , w_u2 = 0, 1, 5, 1, 1
test = highOrder_System(a, b, N, w_u1, w_u2).Method()