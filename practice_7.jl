# 1
function next_repit_placement!(p::Vector{T}, n::T) where T<:Integer
    i = findlast(x->(x < n), p) 
    isnothing(i) && (return nothing)
    p[i] += 1
    p[i+1:end] .= 1 
    return p
end

println("Генерация следующего размещения с повторениями из n элементов {1,2,...,n} по k")
println(next_repit_placement!([1, 1, 1], 3))

println("Генерация всех размещений с повторениями из n элементов {1,2,...,n} по k")
n = 2; k = 3
prev = ones(Int,k)
println(prev)
while !isnothing(prev)
    global prev = next_repit_placement!(prev, n)
    !isnothing(prev) && println(prev)
end

# 2
function next_permute!(p::AbstractVector)
    n = length(p)
    k = 0 
    for i in reverse(1:n - 1) 
        if p[i] < p[i + 1]
            k = i
            break
        end
    end
    k == firstindex(p) - 1 && return nothing 

    i = k + 1
    while i < n && p[i + 1] > p[k]
        i += 1
    end
    p[k], p[i] = p[i], p[k]
    reverse!(@view p[k + 1:end])
    return p
end

println("Генерация следующей перестановки из n элементов {1,2,...,n} по k")
println(next_permute!([1, 2, 3, 4]))
println("Генерация всех перестановок из n элементов {1,2,...,n} по k")
p=[1, 2, 3, 4]
println(p)
while !isnothing(p)
    global p = next_permute!(p)
    !isnothing(p) && println(p)
end

# 3
println("Генерация всех всех подмножеств n-элементного множества {1,2,...,n} 1 способ")

indicator(i::Integer, n::Integer) = reverse(digits(Bool, i; base=2, pad=n))

println("1 способ")
println(indicator(12, 5))
println("1 способ (Все)")
for i in range(0, 2^5 - 1)
    println(indicator(i, 5))
end

function next_indicator!(indicator::AbstractVector{Bool})
    i = findlast(x->(x==0), indicator)
    isnothing(i) && return nothing
    indicator[i] = 1
    indicator[i+1:end] .= 0
    return indicator 
end

println("2 способ")
println(next_indicator!(indicator(12, 5)))
println("2 способ (Все)")
p = indicator(0, 5)
println(p)
while !isnothing(p)
    global p = next_indicator!(p)
    !isnothing(p) && println(p)
end


# 4
function next_indicator!(indicator::AbstractVector{Bool}, k)
    i = lastindex(indicator)
    while indicator[i] == 0
        i -= 1
    end
    m = 0 
    while i >= firstindex(indicator) && indicator[i] == 1 
        m += 1
        i -= 1
    end
    if i < firstindex(indicator)
        return nothing
    end
    indicator[i] = 1
    indicator[i + 1:end] .= 0
    indicator[lastindex(indicator) - m + 2:end] .= 1
    return indicator 
end

println("Генерация всех k-элементных подмножеств n-элементного множества {1, 2, ..., n}")
n = 6
k = 3
a = 1:6
a[findall(next_indicator!([zeros(Bool, n-k); ones(Bool, k)], k))] |> println


# 5
function next_split!(s ::AbstractVector{Int64}, k)
    k == 1 && return (nothing, 0)
    i = k-1
    while i > 1 && s[i-1]>=s[i]
        i -= 1
    end
    s[i] += 1
    r = sum(@view(s[i+1:k]))
    k = i+r-1 
    s[(i+1):(length(s)-k)] .= 1
    return s, k
end

println("Генерация всех разбиений натурального числа на положительные слагаемые")
println(next_split!(ones(Int64, 5), 5))


abstract type AbstractCombinObject
end

Base.iterate(obj::AbstractCombinObject) = (get(obj), nothing)
Base.iterate(obj::AbstractCombinObject, state) = (isnothing(next!(obj)) ? nothing : (get(obj), nothing))

# 6
struct RepitPlacement{N,K} <: AbstractCombinObject
    value::Vector{Int}
    RepitPlacement{N,K}() where {N, K} = new(ones(Int, K))
end

Base.get(p::RepitPlacement) = p.value
next!(p::RepitPlacement{N,K}) where {N, K} = next_repit_placement!(p.value, N)

println("Размещения с повторениями")
for a in RepitPlacement{2,3}() 
    println(a)
end

struct Permute{N} <: AbstractCombinObject
    value::Vector{Int}
    Permute{N}() where N = new(collect(1:N))
end

Base.get(obj::Permute) = obj.value
next!(permute::Permute) = next_permute!(permute.value)

println("Перестановки")
for p in Permute{4}()
    println(p)
end

struct Subsets{N} <: AbstractCombinObject
    indicator::Vector{Bool}
    Subsets{N}() where N = new(zeros(Bool, N))
end

Base.get(sub::Subsets) = sub.indicator
next!(sub::Subsets) = next_indicator!(sub.indicator) 

println("Все подмножества N-элементного множества")
for sub in Subsets{4}()
    println(sub)
end

struct KSubsets{M,K} <: AbstractCombinObject
    indicator::Vector{Bool}
    KSubsets{M, K}() where{M, K} = new([zeros(Bool, length(M)-K); ones(Bool, K)])
end

Base.get(sub::KSubsets) = sub.indicator
next!(sub::KSubsets{M, K}) where{M, K} = next_indicator!(sub.indicator, K) 

for sub in KSubsets{1:6, 3}()
    sub |> println
end

mutable struct NSplit{N} <: AbstractCombinObject
    value::Vector{Int}
    num_terms::Int 
    NSplit{N}() where N = new(vec(ones(Int, N)), N)
end

Base.get(nsplit::NSplit) = nsplit.value[begin:nsplit.num_terms]
function next!(nsplit::NSplit)
    a, b = next_split!(nsplit.value, nsplit.num_terms)
    if isnothing(a) return nothing end
    nsplit.value, nsplit.num_terms = a, b
    get(nsplit)
end

println("Разбиения")
for s in NSplit{5}()
    println(s)
end


function dfs(graph::AbstractDict, vstart::T) where T <: Integer
    stack = [vstart]
    mark = zeros(Bool, length(graph)) 
    mark[vstart] = true
    while !isempty(stack)
        v = pop!(stack)
        for u in graph[v] mark[u] = (!mark[u] ? (push!(stack, u); true) : true) end
    end
    return mark
end

graph1 = Dict{Int64, Vector{Int64}}([(1, [3]), (2, [4]), (3, [1]), (4, [2, 5]), (5, [4])])
graph2 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 4]), (3, [1]), (4, [2, 5]), (5, [4])])
graph3 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 4]), (3, [1, 6]), (4, [2, 5]), (5, [4, 6]), (6, [3, 5])])
graph4 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 3, 4]), (3, [1, 2]), (4, [2, 5]), (5, [4])])
println(dfs(graph2, 1))

function isConnectedGraph(graph::AbstractDict) :: Bool
    res = dfs(graph, 1)
    return (sum(res) == length(res) ? true : false)
end

println(isConnectedGraph(graph1))

function find_connectivity_components(graph::AbstractDict, len = length(graph))
    mark = ones(Bool, len)
    ans = []
    for i in 1:len
        res = mark[i] ? (t = dfs(graph, i);  push!(ans, t); t) : Bool[0]
        mark[findall(res)] .= false
    end
    return ans
end

println(find_connectivity_components(graph1))


function  isBipartite(graph::AbstractDict, len = length(graph)) :: Bool
    color = fill(-1, len)
    queue = []
    for i in 1:len
        if color[i] != -1 continue end
        color[i] = 0
        push!(queue, i)
        while !isempty(queue)
            v = popfirst!(queue)
            if !isnothing(findfirst(isequal(color[v]), color[graph[v]])) return false end
            found = graph[v][findall(isequal(-1), color[graph[v]])]
            color[found] .= (color[v] + 1) % 2
            append!(queue, found)
        end
    end
    return true
end

println(isBipartite(graph3))
println(isBipartite(graph4))
