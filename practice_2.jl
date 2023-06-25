# № 1 ----------------------------------
function fast_pow(a::T, n::Int) where T<:Any
    k=n
    p=a
    t=1
    #ИНВАРИАНТ: p^k*t=a^n
    while k>0
        if iseven(Integer(k))
            k /=2
            p *= p # - это преобразование следует из инварианта
        else
            k -= 1
            t *= p
        end
    end
    return t
end

# № 2 ----------------------------------
function fast_fib(a::Integer)
    matr = [0:1;; 1; 1]
    return ([0;;1]*fast_pow(matr, a))[1, 1]
end

# № 3 ----------------------------------
"""
z=x; t=1; y=0
while z < 1/a || z > a || t > ε 
    if z < 1/a
        z *= a # это перобразование направлено на достижения условия окончания цикла
        y -= t # тогда необходимрсть этого преобразования следует из инварианта цикла
    elseif z > a
        z /= a # это перобразование направлено на достижения условия окончания цикла
        y += t # тогда необходимрсть этого преобразования следует из инварианта цикла
    elseif t > ε
        t /= 2 # это перобразование направлено на достижения условия окончания цикла
        z *= z # тогда необходимрсть этого преобразования следует из инварианта цикла
    end
end
# y - искомое приближенное значение
"""
function log(a, x, e) # a > 1    z^t * a^y = x    
    z = x
    t = 1
    y = 0
    # Пусть y = log_a(x), y_1 = y + g, где g - погрешность. Тогда x = a^y = a^y_1-g = a^y_1 * a^-g = a^y_1 * a^t, t = -g
    while z < 1/a || z > a || t > e  # написать инвариант
        if z < 1/a # Избегаем исчезновение порядка
            z *= a 
            y -= t 
        elseif z > a # Избегаем переполнение порядка
            z /= a
            y += t
        elseif t > e 
            t /= 2 
            z *= z 
        end
    end
    return y
end

# № 4 ----------------------------------
function bisection(f::Function, a, b, epsilon)
    if f(a)*f(b) < 0 && a < b

        f_a = f(a)
        #ИНВАРИАНТ: f_a*f(b) < 0
        while b-a > epsilon
            t = (a+b)/2
            f_t = f(t)
            if f_t == 0
                return t
            elseif f_a*f_t < 0
                b=t
            else
                a, f_a = t, f_t
            end
        end  
        return (a+b)/2
    else
        @warn("Требуются разные знаки на концах отрезка")
    end
end

# № 5 ----------------------------------
cos_(x) = cos(x) - x

bisection(x -> cos(x) - x , 0, 1, 1e-8)



# № 6 ----------------------------------
function newton(f, df, x0, eps=1e-6, num_max=100)
    x = x0
    for i in 1:num_max
        fx = f(x)
        dfx = df(x)
        if abs(fx) < eps
            return x
        else
            x = x - fx / dfx
        end
    end

    error("Превышено максимальное количество итераций.")
end

# № 7 ----------------------------------

println(newton(x -> cos(x) - x, x -> -sin(x) - 1, 0.5, 1e-3, 100))

# № 8 ----------------------------------
p(x) = 6*x^5 - 23*x^4 + 12*x^2 + 86

dp(x) = 24*x - 92*x^3 + 30*x^4

println(newton(p, dp, 5, 1e-3, 1000))
# система уравнений + ???

# Найти решение x^2+y^2=1 x^3y=1
# СХема горнера