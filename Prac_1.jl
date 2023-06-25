function gcd(a::T, b::T) where T<:Integer
    while b > 0
        a, b = b, a % b
    end
    return a
end

function gcd_big(a::T, b::T) where T<:Integer
    u, v = one(T), zero(T); u1, v1 = 0, 1
    #ИНВАРИАНТ:
    while b > 0
        k,r = divrem(a, b)
        a, b = b, r #a - k * b
        u, v, u1, v1 = u1, v1, u - k * u1, v - k * v1
    end
    
    return a, u, v
end

struct Z{T,N}
    a::T
    Z{T,N}(a::T) where {T<:Integer, N} = new(mod(a, N))
end

function inverse(a::Z{T,N}) where {T<:Integer, N}
    if gcd(a.a, N) != 1 
        return nothing
    else
        f, s, d = gcd_big(a.a, N)
        return Z{T,N}(s)
    end 
end

function diaphant_solve(a::T,b::T,c::T) where T<:Integer
    if mod(c,gcd(a,b))!=0
        return nothing
    end
    return gcd_big(a,b)[2:3]
end


Base. +(a::Z{T,N}, b::Z{T,N}) where {T<:Integer, N} = Z{T,N}(a.a + b.a)
Base. -(a::Z{T,N}, b::Z{T,N}) where {T<:Integer, N} = Z{T,N}(a.a - b.a)
Base. *(a::Z{T,N}, b::Z{T,N}) where {T<:Integer, N} = Z{T,N}(a.a * b.a)
Base. -(a::Z{T,N}) where {T<:Integer, N} = Z{T,N}(-a.a)
Base. display(a::Z{T,N}) where {T<:Integer, N} = println(string(a.a))



struct Polynom{T}
    a::T
    Polynom{T}(a) where T = new(a)
end
(p::Polynom)(x)=gorner2(p.a,x)

function Base. +(a::Polynom{T}, b::Polynom{T}) where {T<:Vector} 
    res = if (max(length(a.a),length(b.a))==length(a.a)) a else b  end
    b = if (min(length(a.a),length(b.a))==length(a.a)) a else b  end
    for i in eachindex(b.a) #1:length(a.a)
        res.a[i]+=b.a[i]
        
    end
    return res
end

function Base. -(a::Polynom{T}, b::Polynom{T}) where {T<:Vector} 
    res = if (max(length(a.a),length(b.a))==length(a.a)) a else b  end
    b = if (min(length(a.a),length(b.a))==length(a.a)) a else b  end
    for i in eachindex(b.a) #1:length(a.a)
        res.a[i]-=b.a[i]
        
    end
    return res
end

function Base. *(a::Polynom{T}, b::Polynom{T}) where {T<:Vector} 
    res = Polynom{T}(zeros(Int,length(a.a)+length(b.a)-1))
    for i in eachindex(a.a) #1:length(a.a)
        for j in 0:(length(b.a)-1)

            res.a[i+j]+=a.a[i]*b.a[j+1]
        end
    end
    return res
end

#####
function Base.divrem(A::Polynom{T}, B::Polynom{T}) where T
    B = copy(B.a)
    A = copy(A.a)
    D = Vector{T}(undef, length(A) - length(B) + 1)
    R = Vector{T}(undef, length(B)-1)
    for i in 1:length(A) - length(B) + 1
        D[i] = A[i] / B[1]
        A[i:end] .-= (D[i] .*B)
    end
    i = findfirst(D .!= 0) 
    D = D[i:end]
    i = findfirst(D .!= 0)
    R = A[i:end]
    return Polynom{T}(D), Polynom{T}(R)
end


function Base. /(p::Polynom{T}, q::Polynom{T})where {T<:Number}
    n, m = length(p.a), length(q.a)
    if n < m
        return Polynom([0])
    end
    coeffs = copy(p.a)
    for i in (n-m+1):-1:1
        c = coeffs[i+m-1] / q.coeffs[m]
        for j in 1:m
            coeffs[i+j-1] -= c * q.coeffs[j]
        end
        coeffs[i+m-1] = c
    end
    return Polynom(coeffs[1:(n-m+1)])
end

function Base.:mod(T::Polynom{k},M::Polynom{L})where {k,L}
    vec1 = T.a
    vec2 = M.a
    len1,len2 = length(vec1),length(vec2)
    println(len2,len1)
    if len2 > len1
        return vec1
    else
        vtemp1 = vec1
        vtemp2 = vec2
        while length(vtemp2) <= length(vtemp1)
            vtemp2 = vec2
            for j in 1:length(vtemp1) - length(vtemp2) 
                pushfirst!(vtemp2,0)
            end
            println(vtemp1,vtemp2)
            vtemp2.*= vtemp1[length(vtemp1)] / vtemp2[length(vtemp2)]
            vtemp1 = vtemp1 - vtemp2
            while vtemp1[length(vtemp1)] == 0 && length(vtemp1) != 1
                pop!(vtemp1)
            end
            println(vtemp1,vtemp2)
        end
        return Polynom{L}(vtemp1)
    end    
end 

function Base. >>(a::Polynom{T}, b::Z{M,N}) where {T<:Vector,M<:Integer, N}
    for i in eachindex(a.a) #1:length(a.a)
        
        a.a[i]=mod(a.a[i],N)
    end
    return a
end

function Base. >>(b::Z{M,N},a::Polynom{T}) where {T<:Vector,M<:Integer, N}
    a=Polynom{T}(zeros(Int, N))
    for i in 1:N
        
        a.a[i]=N-i
    end
    return a
end


function gorner(n::Int, a::AbstractVector{T}, t::T) where T
    p=zero(T)
    dp=zero(T)
    i = 1
    while (i < n)
        dp = dp * x + p
        p = p * t + a[i]
    end
    return dp, p
end

function gorner2(a::T, t) where T
    i = 2
    while (i <= length(a))
        a[i] = a[i] + a[i-1]*t 
        i = i + 1;
    end
    pop!(a)
    return a
end



#-------------------------------------
#F = Z{Int, 8}(7)
Q = Z{Int, 5}(9)

#print(inverse(Q),"\n")
#print(diaphant_solve(3,7,1))
T=Polynom{Vector}([1,2,4,5,6,8,7])
F=Polynom{Vector}([1,2,4,5])
#print(F>>T)
print(divrem(T,F))