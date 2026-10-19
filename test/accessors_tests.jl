@testitem "terms of sums" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, NCAdd, terms, factors, coefficient
    NonCommutativeProducts.disable_autosort!()
    a, b, c = Fermion(:a), Fermion(:b), Fermion(:c)
    x = 2 * a * b + 3 * c + 1
    ts = collect(terms(x))
    @test length(ts) == 3
    @test all(t -> t isa NCMul, ts)
    @test sum(terms(x)) == x
    @test Set(ts) == Set([2 * a * b, 3 * c, NCMul(1, Fermion{Symbol}[])])
    # the scalar part appears as a product without factors only when nonzero
    y = 2 * a * b + 3 * c
    @test length(collect(terms(y))) == 2
    @test all(t -> !isempty(factors(t)), terms(y))
    @test sum(terms(y)) == y
    @test sum(terms(a + 1)) == a + 1
    @inferred terms(x)
    # complex coefficients
    z = (1 + 2im) * a + 2im
    @test sum(terms(z)) == z
    # a sum whose terms cancelled keeps only its scalar part
    @test collect(terms(a + 5 - a)) == [NCMul(5, Fermion{Symbol}[])]
end

@testitem "terms of products, numbers and atoms" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, terms, factors, coefficient
    NonCommutativeProducts.disable_autosort!()
    a, b = Fermion(:a), Fermion(:b)
    m = 3 * a * b
    @test collect(terms(m)) == [m]
    @test sum(terms(m)) == m
    @test collect(terms(a)) == [NCMul(a)]
    @test sum(terms(a)) == a
    @test sum(terms(a')) == a'
    @test collect(terms(5)) == [NCMul(5, Any[])]
    @test sum(terms(2.5)) == 2.5
    @test sum(terms(1 + 2im)) == 1 + 2im
    # a product without factors is a scalar term
    @test collect(terms(NCMul(4, Fermion{Symbol}[]))) == [4]
end

@testitem "terms of zeros" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, NCAdd, terms
    NonCommutativeProducts.disable_autosort!()
    a, b = Fermion(:a), Fermion(:b)
    @test isempty(terms(0))
    @test isempty(terms(0.0))
    @test isempty(terms(0 * a))
    @test isempty(terms(zero(a)))
    @test isempty(terms(zero(2 * a * b)))
    @test isempty(terms(zero(a + b)))
    @test isempty(terms(a - a))
    @test isempty(terms(a * b - a * b))
    @test isempty(terms(NCMul(0, Fermion{Symbol}[])))
    @test_throws ArgumentError sum(terms(0))
    # a zero coefficient in the sum is skipped but the other terms remain
    @test collect(terms(0 * a + b)) == [NCMul(b)]
end

@testitem "factors and coefficient" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, terms, factors, coefficient
    NonCommutativeProducts.disable_autosort!()
    a, b = Fermion(:a), Fermion(:b)
    m = 3 * a * b'
    @test factors(m) == [a, b']
    @test factors(m) isa AbstractVector
    @test coefficient(m) == 3
    @test coefficient(NCMul(a)) == 1
    @test coefficient(2.5 * a) == 2.5
    @test isempty(factors(NCMul(7, Fermion{Symbol}[])))
    @test coefficient(NCMul(7, Fermion{Symbol}[])) == 7
    # the result is a copy, so modifying it leaves the product unchanged
    fs = factors(m)
    reverse!(fs)
    @test factors(m) == [a, b']
    @test m == 3 * a * b'
    # terms reconstruct the parts of a sum
    for t in terms(2 * a * b + 5 * b)
        @test t == coefficient(t) * prod(factors(t))
    end
    # a sum with one term has the factors and coefficient of that term
    @test factors(a + b - b) == [a]
    @test coefficient(a + b - b) == 1
    @test coefficient(2 * a + 2 * b - 2 * b) == 2
    @test_throws ArgumentError factors(a + b)
    @test_throws ArgumentError coefficient(a + b)
    @test_throws "terms" NCMul(a + b)
end

@testitem "factors and coefficient of numbers and atoms" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, terms, factors, coefficient, scalar
    NonCommutativeProducts.disable_autosort!()
    a, b = Fermion(:a), Fermion(:b)
    @test coefficient(2) == 2
    @test factors(2) == []
    @test coefficient(a) == 1
    @test factors(a) == [a]
    @test factors(a) isa Vector{typeof(a)}
    @test scalar(NCMul(2, typeof(a)[])) == 2
    err = try
        scalar(a * b)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("not a scalar", err.msg)
    # the accessors agree with the terms of a sum, including the scalar term and atom terms
    x = 2 * a * b + 3 * a + 1
    for t in terms(x)
        @test factors(t) isa AbstractVector
        @test t == coefficient(t) * prod(factors(t); init=1)
    end
    @test coefficient(first(terms(a))) == 1
    @test factors(first(terms(a))) == [a]
end

@testitem "accessors on Tuple-backed products" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, NCAdd, terms, factors, coefficient
    NonCommutativeProducts.disable_autosort!()
    a, b = Fermion(:a), Fermion(:b)
    t = NCMul(2, (b, a'))
    @test t.factors isa Tuple
    @test factors(t) isa Vector
    @test factors(t) == [b, a']
    @test coefficient(t) == 2
    @test collect(terms(t)) == [t]
    @test sum(terms(t)) == t
    # a sum made from Tuple-backed products has Vector-backed terms
    s = t + NCMul(3, (a,)) + 1
    @test sum(terms(s)) == s
    @test all(u -> factors(u) isa Vector, terms(s))
    @test Set(collect(terms(s))) == Set([NCMul(2, [b, a']), NCMul(3, [a]), NCMul(1, Fermion{Symbol}[])])
end
