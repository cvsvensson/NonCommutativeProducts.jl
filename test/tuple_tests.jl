@testitem "Arithmetic and sorting on Tuple-backed products" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, NCAdd, bubble_sort
    f1 = Fermion(1)
    f2 = Fermion(2)
    # each Tuple-backed expression is compared with its Vector-backed counterpart
    t, v = NCMul(1, (f2, f1)), NCMul(1, [f2, f1])
    u, w = NCMul(2, (f1', f2)), NCMul(2, [f1', f2])
    s, sv = t + u, v + w
    for autosort in (false, true)
        autosort ? NonCommutativeProducts.enable_autosort!() : NonCommutativeProducts.disable_autosort!()
        @test t * NCMul(1, [f1']) == v * NCMul(1, [f1'])
        @test NCMul(1, [f1']) * t == NCMul(1, [f1']) * v
        @test t * u == v * w
        @test t * t == v * v
        @test t * f1' == v * f1'
        @test f1' * t == f1' * v
        @test t^2 == v^2
        @test t' == v'
        @test s * s == sv * sv
        @test s * t == sv * v
        @test t * sv == v * sv
        @test s' == sv'
        @test sort(t) == sort(v)
        @test sort!(t) == sort(v)
        @test sort(s) == sort(sv)
        @test bubble_sort(t) == bubble_sort(v)
    end
    NonCommutativeProducts.disable_autosort!()
    # the tuple is not nested as a single factor
    @test (t * NCMul(1, [f1'])).factors == [f2, f1, f1']
    @test (t * f1').factors == [f2, f1, f1']
    @test (f1' * t).factors == [f1', f2, f1]
    # sorting doesn't mutate the original, and the Tuple-backed product keeps its factors
    @test t.factors === (f2, f1)

    T = typeof(t)
    @test zero(T) == 0 && iszero(zero(T))
    @test one(T) == 1
    @test oneunit(t) == 1
    @test zero(t) + t == t
    @test one(t) * t == t
    @test copy(t) == t
end

@testitem "adjoint of mixed-species sums" setup = [Fermions, Bosons] begin
    import NonCommutativeProducts: NCMul, bubble_sort, _autosort
    using Base.ScopedValues: with
    NonCommutativeProducts.@commutative Fermion Boson
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    b = Boson()
    unsorted(f) = with(f, _autosort => false)
    s = unsorted() do
        2 * f1 * b * f2' + 1im * b' * f2 * b + f1' * f2 + 3 * b * f1 + 1
    end
    NonCommutativeProducts.enable_autosort!()
    @test s' == bubble_sort(unsorted(() -> s'))
    @test sort(s)' == bubble_sort(unsorted(() -> sort(s)'))
    @test s'' == sort(s)
    # Tuple-backed mixed products
    p = NCMul(2, (f1, b, f2'))
    @test p' == bubble_sort(unsorted(() -> p'))
    @test p' == NCMul(2, Any[f1, b, f2'])'
end

@testitem "mul_effect changing the atom type during a sort" begin
    import NonCommutativeProducts as NC
    struct Alpha
        id::Int
    end
    struct Beta
        id::Int
    end
    NC.@nc Alpha Beta
    # two equal Alphas fuse into a Beta, which commutes with everything and moves to the right
    bare_atom = Ref(true)
    function NC.mul_effect(a::Alpha, b::Alpha)
        a.id == b.id && return bare_atom[] ? Beta(a.id) : NC.NCMul(Beta(a.id))
        return a.id > b.id ? NC.Swap(-1) : nothing
    end
    NC.mul_effect(::Alpha, ::Beta) = nothing
    NC.mul_effect(::Beta, ::Alpha) = NC.Swap(1)
    NC.mul_effect(a::Beta, b::Beta) = a.id > b.id ? NC.Swap(1) : nothing

    NC.disable_autosort!()
    expected = -NC.NCMul(1, Any[Alpha(1), Alpha(3), Beta(2)])
    for bare in (true, false)
        bare_atom[] = bare
        # the factors start out as Vector{Alpha} and are widened when the Beta appears
        x = NC.NCMul(1, [Alpha(2), Alpha(1), Alpha(2), Alpha(3)])
        @test NC.bubble_sort(x) == expected
        @test NC.bubble_sort(NC.NCMul(1, (Alpha(2), Alpha(1), Alpha(2), Alpha(3)))) == expected
        @test x.factors == [Alpha(2), Alpha(1), Alpha(2), Alpha(3)]
        # in a sum, the other terms keep sorting after the element type of the terms has widened
        y = NC.bubble_sort(x + NC.NCMul(1, [Alpha(3), Alpha(1)]))
        @test y == expected - NC.NCMul(1, [Alpha(1), Alpha(3)])
        @test NC.bubble_sort(NC.NCMul(1, [Alpha(1), Alpha(1), Alpha(2), Alpha(2)])) == NC.NCMul(1, Any[Beta(1), Beta(2)])
    end
end
