from group import Permutation_Group
from random import randint
from lib import plot_bar, chi2test

def Frank(X, N, K):
    e = Permutation_Group(X[0].n)
    S = []
    for i in range(N):
        if i < len(X):
            S.append(X[i])
        else:
            S.append(e)
    for round in range(K):
        i = randint(0,N-1)
        j = i
        while j == i:
            j = randint(0,N-1)
        if randint(1,2) == 1:
            S[i] = S[i] * S[j]
        else:
            S[i] = S[j] * S[i]
    return S[randint(0,N-1)]

if __name__ == "__main__":
    # X = [Permutation_Group(5,{1:2,2:1,3:3,4:4,5:5}), Permutation_Group(5,{1:3,3:1,2:2,4:4,5:5}), Permutation_Group(5,{1:4,4:1,2:2,3:3,5:5}), Permutation_Group(5,{1:5,5:1,2:2,3:3,4:4})]
    X = [Permutation_Group(4,{1:2,2:1,3:3,4:4}), Permutation_Group(4,{1:3,3:1,2:2,4:4}), Permutation_Group(4,{1:4,4:1,2:2,3:3})]
    N = 10
    K = 100
    counter = {}
    for i in range(10000):
        g = str(Frank(X,N,K))
        if g not in counter:
            counter[g] = 1
        else:
            counter[g] += 1
    print(len(counter))
    print("counter =", counter)
    print("chi2 =",chi2test(counter))
    plot_bar(counter)

