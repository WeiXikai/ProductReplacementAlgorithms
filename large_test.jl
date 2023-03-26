import Base.*
import Base.string
using BenchmarkTools

randint(a::Int,b::Int)::Int = round(rand()*(b-a)+a)

function chi2test(d::Dict{String, Int})::Float64
    fo = [x for x in values(d)]
    sum_data = sum(fo)
    fe = [sum_data/length(fo) for i = 1:length(fo)]
    chi2 = sum([(fo[i]-fe[i])^2/fe[i] for i = 1:length(fo)])
    return chi2
end

erf_inv(z::Float64)::Float64 = (1/2*(log(2/(pi*(z-1)^2))-log(log(2/(pi*(z-1)^2)))))^(1/2)

struct Permutation_Group
    n :: Int
    data :: Vector{Int}
end

function group_identity(n::Int)
    Permutation_Group(n, collect(1:n))
end

function inverse(g::Permutation_Group)::Permutation_Group
    new_data = Array{Int}(undef, g.n)
    for i = 1:g.n
        new_data[g.data[i]] = i
    end
    Permutation_Group(g.n, new_data)
end

(*)(g::Permutation_Group, h::Permutation_Group)::Permutation_Group = Permutation_Group(h.n, [g.data[h.data[i]] for i = 1:h.n])

function dfs_Permutation_Group(g::Permutation_Group)::Tuple{Int, Int, Int, String}
    visited = Vector{Bool}(undef, g.n)
    for i = 1:g.n
        visited[i] = false
    end

    number_cycles = 0
    largest_cycle = 1
    number_fix_points = 0
    st = ""
    flag = false

    for i = 1:g.n
        if ! visited[i]
            number_cycles += 1
            len_cycle = 1
            visited[i] = true
            if g.data[i] != i
                flag = true
                st = string(st,"(",string(i))
                temp = g.data[i]
                while ! visited[temp]
                    visited[temp] = true
                    st = string(st,",",string(temp))
                    temp = g.data[temp]
                    len_cycle += 1
                end
                st = string(st,")")
                if len_cycle > 1
                    largest_cycle = len_cycle
                end
            else
                number_fix_points += 1
            end
        end
    end

    if ! flag
        st = "(1)"
    end

    return (number_cycles, largest_cycle, number_fix_points, st)
end

function number_cycles(g::Permutation_Group)::Int
    number_cycles, largest_cycle, number_fix_points, st = dfs_Permutation_Group(g)
    return number_cycles
end

function largest_cycle(g::Permutation_Group)::Int
    number_cycles, largest_cycle, number_fix_points, st = dfs_Permutation_Group(g)
    return largest_cycle
end

function number_fix_points(g::Permutation_Group)::Int
    number_cycles, largest_cycle, number_fix_points, st = dfs_Permutation_Group(g)
    return number_fix_points
end

function string(g::Permutation_Group)::String
    number_cycles, largest_cycle, number_fix_points, st = dfs_Permutation_Group(g)
    return st
end

function Frank(X::Vector{Permutation_Group}, N::Int, K::Int)::Permutation_Group
    e = group_identity(X[1].n)
    S = Vector{Permutation_Group}(undef, N)
    for i = 1:N
        if i <= length(X)
            S[i] = X[i]
        else
            S[i] = e
        end
    end

    for r = 1:K
        i = randint(1,N)
        j = i
        while i == j
            j = randint(1,N)
        end
        if randint(1,2) == 1
            S[i] = S[i] * S[j]
        else
            S[i] = S[j] * S[i]
        end
    end

    return S[randint(1,N)]    
end

mutable struct Product_Replacement_Process
    N :: Int
    X :: Vector{Permutation_Group}
    S :: Int
    state :: String
    good_state :: Any
    data :: Vector{Permutation_Group}
end

function new_Product_Replacement_Process(N :: Int, X :: Vector{Permutation_Group}, S :: Int)::Product_Replacement_Process
    ans = Product_Replacement_Process(N,X,S,"bad",nothing,Vector{Permutation_Group}(undef,N+1))
    for i = 1:N+1
        if i <= length(X)
            ans.data[i] = X[i]
        else
            ans.data[i] = group_identity(X[1].n)
        end
    end
    return ans
end

function get_Product_Replacement_Process(p::Product_Replacement_Process)::Permutation_Group
    return p.data[p.N+1]
end

function update_Product_Replacement_Process!(p::Product_Replacement_Process)::Permutation_Group
    for s = 1:p.S
        i = randint(1,p.N)
        j = i
        while j == i
            j = randint(1,p.N)
        end
        k = randint(1,p.N)
        p.data[i] = p.data[i] * p.data[j]
        p.data[p.N+1] = p.data[p.N+1] * p.data[k]
    end
    return get_Product_Replacement_Process(p)
end

function store_Product_Replacement_Process!(p::Product_Replacement_Process)
    p.good_state = p.data
    return nothing
end

function reset_Product_Replacement_Process!(p::Product_Replacement_Process)
    p.data = p.good_state
    p.state = "good"
    return nothing
end

mutable struct Fix_Length_Window
    state :: Bool
    i :: Int
    N :: Int
    data :: Vector{Int}
end

function new_Fix_Length_Window(N::Int)
    return Fix_Length_Window(false, 0, N, Vector{Int}(undef,N))
end

function push_Fix_Length_Window!(W::Fix_Length_Window, x::Int)
    W.i += 1
    if W.i > W.N
        W.i -= W.N
    end
    W.data[W.i] = x
    if W.i == W.N
        W.state = true
    end
    return nothing
end

function chi2(W1::Fix_Length_Window, W2::Fix_Length_Window, alpha::Float64)::Bool
    C = Dict{Int, Int}()
    B = Dict{Int, Int}()
    E = Dict{Int, Int}()
    for x in W1.data
        if !haskey(C,x)
            C[x] = 1
            B[x] = 1
            E[x] = 0
        else
            C[x] += 1
            B[x] += 1
        end
    end
    for x in W2.data
        if !haskey(C,x)
            C[x] = 1
            B[x] = 0
            E[x] = 1
        else
            C[x] += 1
            E[x] += 1
        end
    end
    r1 = sum(values(B))
    r2 = sum(values(E))
    r = r1+r2
    X = r*(sum([(B[j]-r1*C[j]/r)^2/(r1*C[j]) for j in keys(C)])+sum([(E[j]-r2*C[j]/r)^2/(r2*C[j]) for j in keys(C)]))
    f = length(C) - 1
    if f == 0
        return false
    else
        Y = ((X/f)^(1/3)-1-2/(9*f))/(2/(9*f))^2
        Y_alpha = (2)^(1/2)*erf_inv(1-2*alpha)
        return Y < Y_alpha
    end
end

function Henrik(X::Vector{Permutation_Group}, N::Int, S::Int, m::Int, u::Int, c::Int, a::Int, alpha::Float64, t)::Permutation_Group
    Ta = new_Fix_Length_Window(m)
    Te = new_Fix_Length_Window(m)

    Sampler = new_Product_Replacement_Process(N,X,S)
    Seeker = new_Product_Replacement_Process(N,X,S)

    passed_times = 0
    l = 0
    current_state = "bad"
    previous_Ta = nothing

    while true
        l += 1

        ga = update_Product_Replacement_Process!(Sampler)
        push_Fix_Length_Window!(Ta,t(ga))

        ge = update_Product_Replacement_Process!(Seeker)
        push_Fix_Length_Window!(Te,t(ge))

        if (l % u == 0) && (Ta.state)
            if previous_Ta == nothing
                passed_times = 0
            elseif chi2(Ta, previous_Ta, alpha)
                passed_times += 1
            else
                passed_times = 0
            end
            previous_Ta = Ta
        end

        if passed_times == c
            break
        end

        if (l % u == 0) && (Te.state)
            if chi2(Ta, Te, alpha)
                current_state = "good"
            end
        end

        if (l % a == 0) && (current_state == "good")
            if Seeker.good_state == nothing
                store_Product_Replacement_Process!(Seeker)
            else
                reset_Product_Replacement_Process!(Seeker)
            end
        end
    end

    return get_Product_Replacement_Process(Seeker)
end

function random_permutation(len_X::Int)::Permutation_Group
    ans = group_identity(len_X)
    for i in 1:(len_X-1)
        j = randint(i,len_X)
        if i != j
            temp = ans.data[i]
            ans.data[i] = ans.data[j]
            ans.data[j] = temp
        end
    end
    return ans
end

function compare_Frank_Henrik()
    for len_X in [20,50,100,500,1000] # Length of permutation
        for len_S in [3,5,10,20,50] # number of generators
            S = Vector{Permutation_Group}(undef, len_S)
            for i in 1:len_S
                S[i] = random_permutation(len_X)
            end
            for t in 1:10
                println("Frank ",len_X," ",len_S)
                @time Frank(S, 5*len_S, 50)
                println("Henrik-Number_cycles ",len_X," ",len_S)
                @time Henrik(S, 5*len_S, 5, 5*len_S, 2, 2, len_S, 0.05, number_cycles)
                println("Henrik-number_fix_points ",len_X," ",len_S)
                @time Henrik(S, 5*len_S, 5, 5*len_S, 2, 2, len_S, 0.05, number_fix_points)
                println("Henrik-largest_cycle ",len_X," ",len_S)
                @time Henrik(S, 5*len_S, 5, 5*len_S, 2, 2, len_S, 0.05, largest_cycle)
            end
        end
    end
end

compare_Frank_Henrik()
