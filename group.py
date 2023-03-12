class Permutation_Group():
    def __init__(self, n, data=None):
        self.n = n
        if data == None:
            self.data = {i:i for i in range(1,n+1)}
        else:
            self.data = data
        self.number_cycle = n
        self.number_fix_point = n
        self.largest_cycle = 1
    
    def cycle_number(self):
        self.__str__()
        return self.number_cycle
    
    def fix_point_number(self):
        self.__str__()
        return self.number_fix_point

    def largest_cycle_size(self):
        self.__str__()
        return self.largest_cycle

    def inverse(self):
        ans = {}
        for x in self.data:
            ans[self.data[x]] = x
        return Permutation_Group(self.n,ans)

    def __str__(self):
        visited = [False for i in range(self.n+1)]
        flag = False
        ans = ""
        self.number_cycle = 0
        self.largest_cycle = 1
        self.number_fix_point = 0
        for i in range(1,self.n+1):
            if not visited[i]:
                self.number_cycle += 1
                len_cycle = 1
                visited[i] = True
                if self.data[i] != i:
                    flag = True
                    ans += f'({i}'
                    temp = self.data[i]
                    while not visited[temp]:
                        visited[temp] = True
                        ans += f',{temp}'
                        temp = self.data[temp]
                        len_cycle += 1
                    ans += ')'
                    if len_cycle > 1:
                        self.largest_cycle = len_cycle
                else:
                    self.number_fix_point += 1
        if not flag:
            ans = "(1)"
        return ans

    def __mul__(self, other):
        if self.n != other.n:
            raise "Error: Calculate the multiplication of two elements in two different permutation group."
        return Permutation_Group(self.n,{i:self.data[other.data[i]] for i in range(1,self.n+1)})
    
    def __pow__(self, k):
        if k == 0:
            return Permutation_Group(self.n)
        x = self
        if k < 0:
            x = self.inverse()
        ans = Permutation_Group(self.n)
        for i in range(abs(k)):
            ans = ans * x
        return ans

    def __eq__(self,other):
        return (self.n == other.n) and (self.data == other.data)

    def __ne__(self,other):
        return not self.__eq__(other)

    def order(self):
        counter = 1
        x = self
        e = Permutation_Group(self.n)
        while x != e:
            x = x * self
            counter += 1
        return counter
    
if __name__ == "__main__": # for test
    gL = [Permutation_Group(5, {1:2,2:3,3:5,5:1,4:4}), Permutation_Group(5, {1:3,2:4,4:2,3:5,5:1}), Permutation_Group(5), Permutation_Group(5,{1:2,2:4,4:5,5:3,3:1})]
    for g in gL:
        print(f"g = {g}")
        print(f"g^-1 = {g.inverse()}")
        print(f"#cycles = {g.cycle_number()}")
        print(f"#fix_points = {g.fix_point_number()}")
        print(f"largest_cycle_size = {g.largest_cycle_size()}")
        print()
