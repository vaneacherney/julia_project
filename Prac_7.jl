# 1
function next_repit_placement!(p::Vector{T}, n::T) where T<:Integer
    i = findlast(x->(x < n), p)
    isnothing(i) && (return nothing)
    p[i] += 1
    p[i+1:end] .= 1
    return p
end
 
println("Генерация всех размещений с повторениями из n элементов {1,2,...,n} по k")
println(next_repit_placement!([1, 1, 1], 3))

"""
n = 2; k = 3
p = ones(Int,k)
println(p)
while !isnothing(p)
    p = next_repit_placement!(p,n)
    println(p)
end
"""

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
 
println("Генерация всех размещений с повторениями из n элементов {1,2,...,n} по k")
println(next_permute!([1, 3, 4, 2]))
 

"""
p=[1,2,3,4]
println(p)
while !isnothing(p)
    p = next_permute!(p)
    println(p)
end
"""


# 3
println("Генерация всех всех подмножеств n-элементного множества {1,2,...,n} 1 способ")
 
 
indicator(i::Integer, n::Integer) = reverse(digits(Bool, i; base=2, pad=n))
 
println("1 способ")
println(indicator(12, 5))
 
# 3.2
 
function next_indicator!(indicator::AbstractVector{Bool})
    i = findlast(x->(x==0), indicator)
    isnothing(i) && return nothing
    indicator[i] = 1
    indicator[i+1:end] .= 0
    return indicator 
end
 
println("2 способ")
println(next_indicator!(indicator(12, 5)))
 
"""
n=5; A=1:n
indicator = zeros(Bool, n)
println(indicator)
while !isnothing(indicator)
    A[findall(indicator)] |> println
    indicator = next_indicator!(indicator)
    println(indicator)
end
"""
 
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
println(a[findall(next_indicator!([zeros(Bool, n-k); ones(Bool, k)], k))])
 
"""
n=6; k=3; A=1:n
indicator = [zeros(Bool,n-k); ones(Bool,k)]
A[findall(indicator)] |> println
for !isnothing(indicator)
    indicator = next_indicator!(indicator, k)
    A[findall(indicator)] |> println
end
"""
 
# 5
 
function next_split!(s ::AbstractVector{Int64}, k)
    k == 1 && return (nothing, 0)
    i = k-1 # - это потому что s[k] увеличивать нельзя
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
 
"""
n=5; s=ones(Int, n); k=n
println(s)
while !isnothing(s)
    println(s[1:k])
    s, k = next_split!(s, k)
    println(s)
end
"""
 
# 6
abstract type AbstractCombinObject
end


Base.iterate(obj::AbstractCombinObject) = (get(obj), nothing)
Base.iterate(obj::AbstractCombinObject, state) = (isnothing(next!(obj)) ? nothing : (get(obj), nothing))
 
 
# 6.1
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

 
# 6.2
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

 
# 6.3
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

 
# 6.4
struct KSubsets{M,K} <: AbstractCombinObject
    indicator::Vector{Bool}
    KSubsets{M, K}() where{M, K} = new([zeros(Bool, length(M)-K); ones(Bool, K)])
end
 
Base.get(sub::KSubsets) = sub.indicator
next!(sub::KSubsets{M, K}) where{M, K} = next_indicator!(sub.indicator, K) 
 
for sub in KSubsets{1:6, 3}()
    sub |> println
end
 
# 6.5
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
 
# 7

function dfs(graph::AbstractDict, start::T) where T <: Integer
    stack = [start]
    push!(stack, start)
    visited = falses(length(graph))
    visited[start] = true
    while !isempty(stack)
        v = pop!(stack)
        for u in graph[v] 
            if !visited[u]
                visited[u] = true
                push!(stack, u)
            end
        end
    end
    return visited
end


function bfs(graph::Dict{T, Vector{T}}, start::T) where T<:Integer
    queue = Queue{T}()
    enqueue!(queue, start)
    visited = falses(length(graph))
    visited[start] = true
    while !isempty(queue)
        v = dequeue!(queue)
        for u in graph[v] 
            visited[u] = (!visited[u] ? (enqueue!(queue, u); true) : true)
        end
    end
    return visited
end


graph1 = Dict{Int64, Vector{Int64}}([(1, [3]), (2, [4]), (3, [1]), (4, [2, 5]), (5, [4])])
graph2 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 4]), (3, [1]), (4, [2, 5]), (5, [4])])
graph3 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 4]), (3, [1, 6]), (4, [2, 5]), (5, [4, 6]), (6, [3, 5])])
graph4 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 3, 4]), (3, [1, 2]), (4, [2, 5]), (5, [4])])
println(dfs(graph2, 1))

function is_connected_graph(graph::AbstractDict) :: Bool
    res = dfs(graph, 1)
    return all(res)
end
println(is_connected_graph(graph1))


# 8
function find_connectivity_components(graph::AbstractDict, len = length(graph))
    mark = ones(Bool, len)
    ans = []
    for i in 1:len
        if mark[i]
            t = dfs(graph, i)
            push!(ans, t)
            mark[findall(t)] .= false
        else
            push!(ans, Bool[0])
        end
    end    
    return ans
end
 
println(find_connectivity_components(graph1))
 
# 9
function  isDual(graph::AbstractDict, len = length(graph)) :: Bool
    color = fill(-1, len)
    queue = []
    for i in 1:len
        if color[i] != -1 
            continue
        end
        color[i] = 0
        push!(queue, i)
        while !isempty(queue)
            v = popfirst!(queue)
            if !isnothing(findfirst(isequal(color[v]), color[graph[v]]))
                return false 
            end
            found = graph[v][findall(isequal(-1), color[graph[v]])]
            color[found] .= (color[v] + 1) % 2
            append!(queue, found)
        end
    end
    return true
end

println(isDual(graph1))
println(isDual(graph2))
println(isDual(graph3))
println(isDual(graph4))
