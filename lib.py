from matplotlib import pyplot as plt
from group import Permutation_Group

def plot_bar(counter, label=""):
    plt.title(label)
    plt.bar(range(len(counter)), list(counter.values()), align='center')
    plt.xticks(range(len(counter)), list(counter.keys()))
    plt.show()

def chi2test(counter):
    fo = list(counter.values())
    sum_data = sum(fo)
    fe = [sum_data/len(fo)] * len(fo)
    chi2 = sum([(fo[i]-fe[i])**2/fe[i] for i in range(len(fo))])
    return chi2
