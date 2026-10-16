@testitem "Collecting powers, signed swap" begin
    import NonCommutativeProducts: bubble_sort, @nc, NCMul, mul_effect
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
    ab2 = bubble_sort(ab)
    ab3 = bubble_sort(b * a)
    @test ab == ab2 == -ab3
    @test hash(ab) == hash(ab2)
    @test bubble_sort(ab * a) == -1 * bubble_sort((a * ab))

    op = bubble_sort(prod([a, b, c, d, a, b, c, d, a, b]))
    @test length(op.dict) == 1
    @test map(x -> x.n, only(keys(op.dict)).factors) == 1:4
    @test sum(x -> x.exp, only(keys(op.dict)).factors) == 10
end