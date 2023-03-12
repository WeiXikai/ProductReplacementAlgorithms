from group import Permutation_Group
from random import randint
from matplotlib import pyplot as plt
from lib import plot_bar, chi2test
from collections import deque
from math import log, pi

class Product_Replacement_Process():
    def __init__(self, N, X, S):
        self.N = N
        self.X = X
        self.S = S
        e = Permutation_Group(X[0].n)
        self.A = [e]
        self.state = "bad"
        self.good_state = None
        for i in range(N):
            if i < len(X):
                self.A.append(X[i])
            else:
                self.A.append(e)
    
    def update(self):
        for s in range(self.S):
            i = randint(1,self.N)
            j = i
            while i == j:
                j = randint(1,self.N)
            k = randint(1,self.N)
            self.A[i] = self.A[i]*self.A[j]
            self.A[0] = self.A[0]*self.A[k]
        return self.A[0]
    
    def get_element(self):
        return self.A[0]
    
    def store_good_state(self):
        self.good_state = self.A
    
    def reset_good_state(self):
        self.A = self.good_state
        self.state = "good"

def chi2(T1, T2, alpha):
    C = {}
    B = {}
    E = {}
    for x in T1:
        if x not in C:
            C[x] = 1
            B[x] = 1
            E[x] = 0
        else:
            C[x] += 1 
            B[x] += 1
    for x in T2:
        if x not in C:
            C[x] = 1
            B[x] = 0
            E[x] = 1
        else:
            C[x] += 1
            E[x] += 1
    r1 = sum(B.values())
    r2 = sum(E.values())
    r = r1+r2
    X = r*(sum([(B[j]-r1*C[j]/r)**2/(r1*C[j]) for j in C])+sum([(E[j]-r2*C[j]/r)**2/(r2*C[j]) for j in C]))

    f = len(C) - 1 
    if f == 0:
        return False
    Y = ((X/f)**(1/3)-1-2/(9*f))/(2/(9*f))**2

    erf_inv = lambda z: (1/2*(log(2/(pi*(z-1)**2))-log(log(2/(pi*(z-1)**2)))))**(1/2)
    Y_alpha = (2)**(1/2)*erf_inv(1-2*alpha)

    return Y < Y_alpha
            

def Henrik(X, N, S, m, u, c, a, alpha, t):
    Ta = deque([])
    Te = deque([])
    
    Sampler = Product_Replacement_Process(N,X,S)
    Seeker = Product_Replacement_Process(N,X,S)
    
    passed_times = 0
    l = 0
    current_state = "bad"
    previous_Ta = None
    
    while True:
        l += 1

        ga = Sampler.update()
        Ta.append(t(ga))
        if len(Ta) > m:
            Ta.popleft()
            
        ge = Seeker.update()
        Te.append(t(ge))
        if len(Te) > m:
            Te.popleft()
        
        if l % u == 0 and len(Ta) == m:
            if previous_Ta is None:
                passed_times = 0
            elif chi2(list(Ta),list(previous_Ta),alpha):
                passed_times += 1
            else:
                passed_times = 0
            previous_Ta = Ta
        
        if passed_times == c:
            break
            
        if l % u == 0 and len(Te) == m:
            if chi2(list(Ta),list(Te),alpha):
                current_state = "good"
        
        if l % a == 0 and current_state == "good":
            if Seeker.good_state is None:
                Seeker.store_good_state()
            else:
                Seeker.reset_good_state()
        
    return Seeker.get_element()


if __name__ == "__main__":
    X = [Permutation_Group(4,{1:2,2:1,3:3,4:4}), Permutation_Group(4,{1:3,3:1,2:2,4:4}), Permutation_Group(4,{1:4,4:1,2:2,3:3})]
    N = 10
    m = 10
    u = 2
    c = 2
    a = 10
    alpha = 0.05
    # t = lambda g:g.cycle_number()
    # t = lambda g:g.fix_point_number()
    t = lambda g:g.largest_cycle_size()
    S = 5
    counter = {}

    for i in range(10000):
        g = str(Henrik(X, N, S, m, u, c, a, alpha, t))
        if g not in counter:
            counter[g] = 1
        else:
            counter[g] += 1
    print(len(counter))
    print("counter =", counter)
    print("chi2 = ", chi2test(counter))
    plot_bar(counter)
