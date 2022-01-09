from numpy import exp
import matplotlib.pyplot as plt

class ode_system01:
    
    def __init__(self, t, end, u1, u2, h):
        self.t = t
        self.temp_t = t
        self.end = end
        self.u1 = u1
        self.u2 = u2
        self.temp_u2 = 0
        self.h = h
        self.u1_lst = []
        self.y_lst = []
        self.t_lst = []
    
    def y(self):
        ans = exp(-self.t)
        return ans
    
    def u2_d(self):
        ans = -51*self.u2 - 50*self.u1
        return ans

    def u1_d(self):
        ans = self.temp_u2
        return ans

    def u2i(self):
        ans = self.u2 + self.h * self.u2_d()
        return ans

    def u1i(self):
        ans = self.u1 + self.h * self.u1_d()
        return ans
         
    def Method(self, x = 0):
        while self.t <= self.end:
            if x:
                print("%.7f %.7f" %(self.y(), self.u1))
            self.u1_lst.append(self.u1)
            self.y_lst.append(self.y())
            self.t_lst.append(self.t)
            self.temp_u2 = self.u2
            self.u2 = self.u2i()
            self.u1 = self.u1i()
            self.t += self.h
        self.t = self.temp_t
    
    def plot(self):
        plt.plot(self.t_lst, self.u1_lst)
        plt.plot(self.t_lst, self.y_lst)
        plt.title("When t = %.2f" % self.h)
        plt.show()
    