@testitem "Collecting powers, signed swap" begin
    import NonCommutativeProducts: @nc, NCMul, mul_effect
    NonCommutativeProducts.disable_autosort!()
    struct NCInt
        n::Int
        exp::Int
    end
    NCInt(n::Int) = NCInt(n, 1)
    Base.show(io::IO, x::NCInt) = print(io, "[$(x.n)^$(x.exp)]")
    Base.:(==)(a::NCInt, b::NCInt) = a.n == b.n
    @nc NCInt
    function mul_effect(a::NCInt, b::NCInt)
        (a.n == b.n && return NCInt(a.n, a.exp + b.exp)) # Collect powers
        a.n > b.n && return -b * a # swap and multiply by -1
        return nothing # Do nothing
    end

    a, b, c, d = NCInt.(1:4)
    ab = a * b
    @test ab == NCMul(1, [a, b])
    ab2 = sort(ab)
    ab3 = sort(b * a)
    @test ab == ab2 == -ab3
    @test hash(ab) == hash(ab2)
    @test sort(ab * a) == -1 * sort((a * ab))

    op = sort(prod([a, b, c, d, a, b, c, d, a, b]))
    @test length(op.dict) == 1
    @test map(x -> x.n, only(keys(op.dict)).factors) == 1:4
    @test sum(x -> x.exp, only(keys(op.dict)).factors) == 10
end

@testitem "sort of a product returns a sum" setup = [Fermions] begin
    import NonCommutativeProducts
    import NonCommutativeProducts: NCAdd, NCMul
    NonCommutativeProducts.disable_autosort!()
    f1, f2 = Fermion(1), Fermion(2)
    for x in (NCMul(3, Any[]), 2 * f1, f1 * f2, f2 * f1, f2' * f1 * f2)
        @test sort(x) isa NCAdd
        @test sort(x) == sort(NCAdd(x))
    end
    @test sort(NCMul(3, Any[])) == 3
    @test sort(2 * f1) == 2 * f1
end

@testitem "Sorted products of atoms with a common supertype have concrete keys" begin
    import NonCommutativeProducts: @nc, mul_effect, with_autosort
    abstract type BraKet end
    struct Ket <: BraKet
        n::Int
    end
    struct Bra <: BraKet
        n::Int
    end
    @nc Ket Bra
    mul_effect(a::BraKet, b::BraKet) = a isa Bra && b isa Ket ? Int(a.n == b.n) : nothing

    H = sum(Ket(n) * Bra(n + 1) + Ket(n + 1) * Bra(n) for n in 1:4)
    ψ = sum(n * 1.0 * Ket(n) for n in 1:5)
    Hψ = with_autosort(() -> H * ψ)
    @test Hψ isa typeof(ψ)
    @test Hψ == sum((n + 1) * Ket(n) + n * Ket(n + 1) for n in 1:4)
    unsorted = with_autosort(() -> H * ψ, false)
    @test keytype(sort(unsorted).dict) == keytype(ψ.dict)
end
