from Leedham import Leedham
from lib import chi2test
from group import Permutation_Group
from matplotlib import pyplot as plt
from random import randint

X = [Permutation_Group(4,{1:2,2:1,3:3,4:4}), Permutation_Group(4,{1:3,3:1,2:2,4:4}), Permutation_Group(4,{1:4,4:1,2:2,3:3})]

N = 5
K = 20

recorder_x = []
recorder_y = []
for Number_non_zero in range(1,30):
    recorder_x.append(Number_non_zero)
    Counter = 0
    for test_number in range(20):
        T = [[[0 for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]
        for t in range(Number_non_zero):
            while True:
                i = randint(1,N)
                j = i
                while i == j:
                    j = randint(1,N)
                k = randint(1,N)
                if T[i][j][k] == 0:
                    T[i][j][k] = 1/Number_non_zero
                    break
        counter = {'(1,2,4,3)': 0, '(1,4,3,2)': 0, '(1)': 0, '(1,4,3)': 0, '(1,4)': 0, '(1,3,2,4)': 0, '(2,3)': 0, '(1,4,2,3)': 0, '(3,4)': 0, '(1,3,2)': 0, '(1,2,3)': 0, '(1,2,3,4)': 0, '(1,4)(2,3)': 0, '(1,3)(2,4)': 0, '(1,2)(3,4)': 0, '(2,4,3)': 0, '(1,3)': 0, '(1,4,2)': 0, '(1,2,4)': 0, '(2,3,4)': 0, '(1,2)': 0, '(1,3,4,2)': 0, '(1,3,4)': 0, '(2,4)': 0}
 
        for i in range(200):
            g = str(Leedham(X,N,K,T))
            counter[g] += 1
        chi2 = chi2test(counter)
        passed = (chi2<35.172)
        if passed:
            Counter += 1
        print(f"N={N}, K={K}, Number_non_zero={Number_non_zero}, test_number={test_number}, chi2={chi2}, passed={passed}")
    recorder_y.append(Counter/20)
plt.plot(recorder_x,recorder_y)
plt.xlabel("The number of element which are non-zero")
plt.ylabel("Possibility to pass the chi^2 test")
plt.show()
