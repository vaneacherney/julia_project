#Practic3 
#1 Is simple? 
#2 Решето эратосфена? Все простые, не привосходящие n (1-n)
#3 n = p1^k1.... factor(n) = [p1...pm][k1...km]
#4 a1..an  Mn = 1/n Dn = SUM(i=1 to n)(ai-Mn)^2 
# Дерево

function isSimple(x::Int)
    a = x
    while(a >= sqrt(x))
        a = a - 1
        if (mod(x,a) == 0)
            return false
        end
    end
    return true
end

function findall!(x::SharedArray{Bool})
    a = Vector
end

function resheto(x::Int)
    a = SharedArray{Bool}(x)
    i = 2
    while (i < length(a))
        a[i] = 1
        i = i + 1
    end
    b = 2
    while b < x
        if (a[b] != 0)         
            i = b * b
            while (i < x)
                a[i] = 0
                i = i + b
            end
        end
        b = b + 1
    end
    return findall(a)
end

function resheto_spec(x::Int)
    a = SharedArray{Bool}(x)
    i = 2
    while (i < length(a))
        if (mod(x,i) == 0)
            a[i] = 1
        end
        i = i + 1
    end
    b = 2
    while b < x
        if (a[b] != 0)         
            i = b * b
            while (i < x)
                a[i] = 0
                i = i + b
            end
        end
        b = b + 1
    end
    return findall(a)
end

function crat(x::Int, y::Int)
    i = 0
    while (x > 1 && mod(x,y) == 0)
        x = x/y
        i = i + 1
    end
    return i
end

function factor(x::Int)
    a = resheto_spec(x)
    if (length(a) == 0)
        return ([x],[1])
    end
    b = SharedArray{Int}(length(a))
    i = 1
    while i <= length(a)
        b[i] = crat(x, a[i])
        i = i + 1
    end
    return (a,b)
end

function factorization(n)
    factors = []
    d = 2
    while n > 1
        while n % d == 0
            push!(factors, d)
            n /= d
        end
        d += 1
        if d * d > n
            if n > 1
                push!(factors, n)
                break
            end
        end
    end
    return factors
end

function standard_deviation(data)
    n = length(data)
    mean_val = sum(data) / n
    deviations = [x - mean_val for x in data]
    std_dev = sqrt(sum(deviations .^ 2) / n)
    return mean_val, std_dev
end

function meanstd(aaa)
    T = eltype(aaa)
    n = 0; s¹ = zero(T); s² = zero(T)
    for a ∈ aaa
    n += 1; s¹ += a; s² += a*a
    end
    mean = s¹ ./ n
    return mean, sqrt(s²/n - mean*mean)
end


# № 5 ----------------------------------
#5 Взаимные преобразования различных способов представления деревьев
struct Node
    index :: Int
    children :: Vector{Union{Nothing,Node}}
end
function convert!( arr :: Vector, tree :: Dict{Int,Vector})
    isempty(arr) && return

    list = []

    for subarr in arr[1:end-1]
        if isempty(subarr)
            push!(list,nothing)
            continue
        end
        if typeof(subarr) <: Int
            push!(list,subarr)
            continue
        end
        push!(list,subarr[end])
        convert!(subarr,tree)
    end

    tree[arr[end]] = list

    return tree
end

function convert!(tree :: Dict{Int,Vector}; root ::  Union{Int,Nothing}) :: Union{Vector,Int}
    arr = []
    isnothing(root) && return []
    !(root in keys(tree)) && return root
    for subroot in tree[root]
        push!(arr,convert!(tree; root = subroot))
    end
    push!(arr,root)
    return arr
end

function convert!( tree :: Dict{Int,Vector}, root :: Union{Int,Nothing}) ::Union{Node,Nothing}

    isnothing(root) && return nothing
    !(root in keys(tree)) && return Node(root,[])
    node = Node(root,[])

    for sub_root in tree[root]
        push!(node.children, convert!(tree, sub_root))
    end

    return node
end

function convert!( node :: Node) :: Union{Vector,Int}
    arr = []
    length(node.children)==0 && return node.index
    for child in node.children
        if isnothing(child)
            push!(arr, [])
            continue
        end
        push!(arr,convert!(child))
    end
    push!(arr,node.index)
    return arr

end
function convert!(node :: Node, tree :: Dict{Int, Vector}) :: Union{Dict{Int,Vector},Int}
    list = []
    for child in node.children
        if isnothing(child)
            push!(list, nothing)
            continue
        end
        push!(list,child.index)
        length(child.children) != 0 && convert!(child,tree)
    end
    tree[node.index] = list
    return tree
end

#----------------------------------------------------------------------------------------------------------#
arr = [[[[],[],6], [], 2], [[10,11,4], [[],[],5], 3],1]
tree = Dict{Int,Vector}();
tree = convert!(arr, tree)

display(tree)

_arr = convert!(tree; root = 1)
println(_arr)

node = convert!(tree, 1)
println(node)

_arr = convert!(node)
println(_arr)

tree = convert!(node,tree)
display(tree)

#--------
struct VectorTree
    data::Vector{Vector{Int}}
end

function treeHeight(tree::VectorTree, node::Int, h::Int)
    hi = h+1
    if(tree.data[node]!=[])
        for i in 1:length(tree.data[node])
            h = max(h, treeHeight(tree, tree.data[node][i], hi))
        end
        hi = h
    end
    return hi
end

function treeVal(tree::VectorTree)
    v = 0
    for i in 1:length(tree.data)
        v = max(v,length(tree.data[i]))
    end
    return v
end

function treeLeafes(tree::VectorTree, node::Int, l::T) where T
    if(tree.data[node]!=[])
        for i in 1:length(tree.data[node])
            treeLeafes(tree, tree.data[node][i], l)
        end
    else
        l.x+=1
    end
end

function treeNodes(tree::VectorTree, node::Int, n::T) where T
    if(tree.data[node]!=[])
        n.x+=1
        for i in 1:length(tree.data[node])
            treeNodes(tree, tree.data[node][i], n)
        end
    end
end

function treeMidHeight(tree::VectorTree, node::Int, h::Int, heights::T) where T
    h+=1
    if(tree.data[node]!=[])
        for i in 1:length(tree.data[node])
            treeMidHeight(tree, tree.data[node][i], h, heights)
        end
    else
        pushfirst!(heights.x, h)
    end
end


#1) проверка простоты числа, 
#2) решето Эратосфена, 
#) факторизация, 
#4) среднее квадратическое отклонение за один проход, 
#5) взаимные преобразования различных способов представления деревьев, 
#6) 5 рекурсивных алгоритмов про деревья (высота дерева, число вершин, число листьев, максимальная валентность, средняя длина пути)

function iSTest()
    tree = [
        [2, 3], # узел 1 имеет потомков 2 и 3
        [4, 5], # узел 2 имеет потомков 4 и 5
        [6, 7], # узел 3 имеет потомков 6 и 7
        [],    # узел 4 не имеет потомков
        [],    # узел 5 не имеет потомков
        [],    # узел 6 не имеет потомков
        [8, 9],  # узел 7 имеет потомков 8 и 9
        [10],   # узел 8 имеет потомка 10
        [],     # узел 9 не имеет потомков
        []     # узел 10 не имеет потомков
    ]
    
    println(treeHeight(VectorTree(tree),1,0))
    println(treeVal(VectorTree(tree)))
    l = Ref{Int}(0)
    treeLeafes(VectorTree(tree),1,l)
    println(l.x)
    l.x = 0
    treeNodes(VectorTree(tree),1,l)
    println(l.x)
    h = Ref{Vector{Int}}([])
    treeMidHeight(VectorTree(tree),1,0,h)
    println(sum(h.x)/l.x)

    println()
    tre = [
        [2],
        [3],
        [4,5,6],
        [],
        [],
        []
    ]
    println(treeHeight(VectorTree(tre),1,0))
    println(treeVal(VectorTree(tre)))
    l = Ref{Int}(0)
    treeLeafes(VectorTree(tre),1,l)
    println(l.x)
    l.x = 0
    treeNodes(VectorTree(tre),1,l)
    println(l.x)
    h.x = []
    treeMidHeight(VectorTree(tre),1,0,h)
    println(sum(h.x)/l.x)
end