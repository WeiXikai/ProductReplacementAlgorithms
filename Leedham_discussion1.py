from Leedham import Leedham_s
from lib import chi2test
from group import Permutation_Group
from matplotlib import pyplot as plt

X = [Permutation_Group(4,{1:2,2:1,3:3,4:4}), Permutation_Group(4,{1:3,3:1,2:2,4:4}), Permutation_Group(4,{1:4,4:1,2:2,3:3})]

P_x = []
P_y = []
N_x = []
N_y = []

for N in range(3,11):
    for K in range(1,51):
        counter = {}
        for i in range(1000):
            g = str(Leedham_s(X,N,K))
            if g not in counter:
                counter[g] = 1
            else:
                counter[g] += 1
        chi2 = chi2test(counter)
        passed = (chi2<35.172)
        if passed:
            P_x.append(N)
            P_y.append(K)
        else:
            N_x.append(N)
            N_y.append(K)
        print(f"N={N}, K={K}, chi2={chi2}, passed={passed}")

plt.scatter(P_y,P_x, color="green", label="passed", marker="o")
plt.scatter(N_y,N_x, color="red", label="not passed", marker="x")
plt.legend()
plt.xlabel("K")
plt.ylabel("N")
plt.show()
