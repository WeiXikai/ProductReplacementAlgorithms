from group import Permutation_Group
from random import randint
from matplotlib import pyplot as plt
from lib import plot_bar, chi2test

def choose_random_from_tensor(T):
    N = len(T) -1
    sum_number = 0
    random_number = 1.0 * randint(1,10000) / 10000
    for i in range(1, N+1):
        for j in range(1, N+1):
            for k in range(1, N+1):
                sum_number += T[i][j][k]
                if sum_number >= random_number:
                    return (i,j,k)
    return (N-1, N, 1)

def Leedham_s(X, N, K):
    e = Permutation_Group(X[0].n)
    S = [e]
    for i in range(N):
        if i < len(X):
            S.append(X[i])
        else:
            S.append(e)
    for round in range(K):
        i = randint(1,N)
        j = i
        while i == j:
            j = randint(1,N)
        k = randint(1,N)
        S[i] = S[i]*S[j]
        S[0] = S[0]*S[k]
    return S[0]


def Leedham(X, N, K, T):
    e = Permutation_Group(X[0].n)
    S = [e]
    for i in range(N):
        if i < len(X):
            S.append(X[i])
        else:
            S.append(e)
    for round in range(K):
        (i,j,k) = choose_random_from_tensor(T)
        S[i] = S[i]*S[j]
        S[0] = S[0]*S[k]
    return S[0]


if __name__ == "__main__":
    X = [Permutation_Group(4,{1:2,2:1,3:3,4:4}), Permutation_Group(4,{1:3,3:1,2:2,4:4}), Permutation_Group(4,{1:4,4:1,2:2,3:3})]
    N = 10
    K = 100
    counter = {}
    T = [[[0 for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]
    for i in range(1,N+1):
        for j in range(1,N+1):
            for k in range(1,N+1):
                if i != j:
                    T[i][j][k] = 1.0/N/(N-1)/N

    for i in range(10000):
        g = str(Leedham(X,N,K,T))
        if g not in counter:
            counter[g] = 1
        else:
            counter[g] += 1
    print(len(counter))
    print("counter =", counter)
    print("chi2 = ", chi2test(counter))
    plot_bar(counter)
