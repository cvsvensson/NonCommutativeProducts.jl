@testmodule ReferenceSort begin
    # A slow reference for sorting. An expression is a Dict mapping a vector of factors to its coefficient, and it is
    # normalized by repeatedly expanding the first neighbouring pair of a product for which mul_effect is not nothing.
    using NonCommutativeProducts
    import NonCommutativeProducts: NCMul, NCAdd, AddTerms, Swap, mul_effect, with_autosort, terms, factors, coefficient
    export ref, ref_sort, ref_mul, ref_adjoint, ref_add, as_ref, rand_term, rand_sum

    const Ref_ = Dict{Vector{Any},Any}

    function addto!(d, fs, c)
        iszero(c) && return d
        v = get(d, fs, 0) + c
        iszero(v) ? delete!(d, fs) : (d[fs] = v)
        return d
    end
    # The reference form of a sum of terms (coeff, factors) and a scalar
    function ref(termlist, scalar=0)
        d = Ref_()
        addto!(d, Any[], scalar)
        for (c, fs) in termlist
            addto!(d, collect(Any, fs), c)
        end
        return d
    end

    # The value of mul_effect(a, b) as a list of (λ, replacement), meaning a*b = Σ λ * replacement
    expand(e::Number, a, b) = [(e, Any[])]
    expand(e::Swap, a, b) = [(e.λ, Any[b, a])]
    expand(e::AddTerms, a, b) = reduce(vcat, [expand(t, a, b) for t in e.terms]; init=Tuple{Any,Vector{Any}}[])
    expand(e::NCMul, a, b) = [(coefficient(e), collect(Any, factors(e)))]
    expand(e::NCAdd, a, b) = [(coefficient(t), collect(Any, factors(t))) for t in terms(e)]
    expand(e, a, b) = [(1, Any[e])]

    effect(a, b) = with_autosort(() -> mul_effect(a, b), false)

    function ref_sort(d)
        out = Ref_()
        work = collect(pairs(d))
        while !isempty(work)
            fs, c = pop!(work)
            iszero(c) && continue
            i = findfirst(i -> !isnothing(effect(fs[i], fs[i+1])), 1:length(fs)-1)
            if isnothing(i)
                addto!(out, fs, c)
            else
                for (λ, rep) in expand(effect(fs[i], fs[i+1]), fs[i], fs[i+1])
                    push!(work, Any[fs[1:i-1]; rep; fs[i+2:end]] => c * λ)
                end
            end
        end
        return out
    end
    function ref_mul(x, y)
        d = Ref_()
        for (f, c) in x, (g, e) in y
            addto!(d, Any[f; g], c * e)
        end
        return d
    end
    ref_adjoint(x) = ref([(c', reverse(map(adjoint, fs))) for (fs, c) in x])
    function ref_add(x, y)
        d = copy(x)
        for (fs, c) in y
            addto!(d, fs, c)
        end
        return d
    end

    # The package's result in reference form
    as_ref(x::Number) = ref((), x)
    as_ref(x) = ref([(coefficient(t), factors(t)) for t in terms(x)])

    randcoeff(rng) = rand(rng, -3:3)
    rand_term(rng, atoms, n) = (randcoeff(rng), rand(rng, atoms, n))
    # A random sum as the package's (unsorted) NCAdd and as its reference form
    function rand_sum(rng, atoms, maxlen)
        termlist = [rand_term(rng, atoms, rand(rng, 1:maxlen)) for _ in 1:rand(rng, 1:4)]
        s = rand(rng, Bool) ? randcoeff(rng) : 0
        x = with_autosort(false) do
            s + sum(NCMul(c, fs) for (c, fs) in termlist)
        end
        return x, ref(termlist, s)
    end
end

@testitem "Sorting agrees with a reference" setup = [ReferenceSort, Fermions, Bosons] begin
    using Random: Xoshiro
    import NonCommutativeProducts: NCMul, with_autosort, mul!!
    rng = Xoshiro(1234)
    # Fermions split into sums with a scalar part, and bosons into AddTerms
    families = [
        "fermions" => [Fermion(l, c) for l in 1:2 for c in (false, true)],
        "bosons" => [Boson(e) for e in (-1, 1)],
    ]
    # Each test lists the random inputs for which the package and the reference disagree
    for (name, atoms) in families
        @testset "$name" begin
            products = [rand_term(rng, atoms, rand(rng, 1:6)) for _ in 1:20]
            sums = [(rand_sum(rng, atoms, 3), rand_sum(rng, atoms, 3), rand_sum(rng, atoms, 3)) for _ in 1:20]
            @test [(c, fs) for (c, fs) in products if as_ref(sort(NCMul(c, fs))) != ref_sort(ref([(c, fs)]))] == []
            @test [x for ((x, rx), _, _) in sums if as_ref(sort(x)) != ref_sort(rx)] == []
            @test [(x, y) for ((x, rx), (y, ry), _) in sums
                   if as_ref(with_autosort(() -> x * y)) != ref_sort(ref_mul(rx, ry))] == []
            @test [x for ((x, rx), _, _) in sums if as_ref(with_autosort(() -> x')) != ref_sort(ref_adjoint(rx))] == []
            # mul!! mutates its first argument, so it gets a fresh copy
            @test [(z, x, y) for ((x, rx), (y, ry), (z, rz)) in sums
                   if as_ref(with_autosort(() -> mul!!(z + 0, x, y))) != ref_sort(ref_add(rz, ref_mul(rx, ry)))] == []
        end
    end
end

@testitem "Inner products agree with a reference" setup = [ReferenceSort] begin
    using Random: Xoshiro
    using LinearAlgebra: dot
    using VectorInterface: inner
    import NonCommutativeProducts: @nc, mul_effect, contraction_key, with_autosort, NCMul, AddTerms, Swap

    # Kets |n⟩ (bras when bra is true), keyed by n, and ladder operators with an unnormalized action
    struct RKet
        n::Int
        bra::Bool
    end
    Base.adjoint(k::RKet) = RKet(k.n, !k.bra)
    struct Ladder
        up::Bool
    end
    Base.adjoint(a::Ladder) = Ladder(!a.up)
    @nc RKet Ladder
    mul_effect(a::RKet, b::RKet) = a.bra && !b.bra ? Int(a.n == b.n) : nothing
    mul_effect(a::Ladder, k::RKet) = k.bra ? nothing : a.up ? RKet(k.n + 1, false) : k.n * RKet(k.n - 1, false)
    mul_effect(k::RKet, a::Ladder) = !k.bra ? nothing : a.up ? k.n * RKet(k.n - 1, true) : RKet(k.n + 1, true)
    mul_effect(a::Ladder, b::Ladder) = !a.up && b.up ? AddTerms((Swap(1), 1)) : nothing
    contraction_key(k::RKet) = k.bra ? nothing : k.n

    rng = Xoshiro(5678)
    randcoeff() = complex(rand(rng, -3:3), rand(rng, -3:3))
    # A sum of products of ladder operators applied to kets, as an NCAdd and in reference form
    function rand_state()
        termlist = [(randcoeff(), Any[rand(rng, Ladder.((false, true)), rand(rng, 0:3)); RKet(rand(rng, 0:3), false)])
                    for _ in 1:rand(rng, 1:5)]
        x = with_autosort(false) do
            sum(NCMul(c, fs) for (c, fs) in termlist) + 0
        end
        return x, ref(termlist)
    end
    # x' * y sorts to a scalar, which inner and dot must equal
    agrees((x, rx), (y, ry)) = let sorted = ref_sort(ref_mul(ref_adjoint(rx), ry))
        all(isempty, keys(sorted)) && inner(x, y) == dot(x, y) == get(sorted, Any[], 0)
    end
    cases = [(rand_state(), rand_state()) for _ in 1:30]
    @test [(x, y) for ((x, rx), (y, ry)) in cases if !agrees((x, rx), (y, ry))] == []
end
