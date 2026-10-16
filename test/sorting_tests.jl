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
