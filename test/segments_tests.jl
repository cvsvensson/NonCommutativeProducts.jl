@testitem "Mixed atom types sort the same in any factor container" setup = [Fermions, Bosons] begin
    import NonCommutativeProducts
    import NonCommutativeProducts: NCMul
    NonCommutativeProducts.disable_autosort!()
    NonCommutativeProducts.@commutative Fermion Boson
    f1, f2, f3 = Fermion(1), Fermion(Int32(2)), Fermion(3)
    b = Boson()

    # Fermion{Int} and Fermion{Int32} share the abstract factor type Fermion
    for F in (Fermion, Union{Fermion{Int},Fermion{Int32}}, Any)
        @test sort(NCMul(1, F[f3', f2, f1])) == -NCMul(1, Any[f3', f1, f2]) == sort(NCMul(1, Any[f3', f2, f1]))
        @test sort(NCMul(1, F[f2, f1])) == -NCMul(1, Any[f1, f2])
    end

    # @commutative orders fermions before bosons, and b b† = b† b + 1
    expected = f1 * f3 * b' * b + f1 * f3
    for F in (Any, Union{Fermion{Int},Boson})
        @test sort(NCMul(1, F[f1, b, f3, b'])) == expected
    end
    @test sort(f1 * b) == f1 * b
    @test sort(b * f1) == f1 * b
end

@testitem "Segments of mixed products" setup = [Fermions, Bosons] begin
    import NonCommutativeProducts
    import NonCommutativeProducts: NCMul, atom_group, foreach_segment, split_by_group, segment_ranges
    f1, f2, f3 = Fermion(1), Fermion(2), Fermion(3)
    b = Boson()
    segments(x) = (out = []; foreach_segment((g, seg) -> push!(out, (g, collect(seg), eltype(seg))), x); out)

    # Interleaved groups give one segment per run
    x = NCMul(2, Any[f1, f2, b, b', f3])
    @test segments(x) == [(Fermion{Int}, [f1, f2], Fermion{Int}), (Boson, [b, b'], Boson), (Fermion{Int}, [f3], Fermion{Int})]
    @test collect(segment_ranges(x)) == [Fermion{Int} => 1:2, Boson => 3:4, Fermion{Int} => 5:5]
    parts = split_by_group(x)
    @test parts == [Fermion{Int} => NCMul(2, [f1, f2]), Boson => NCMul(1, [b, b']), Fermion{Int} => NCMul(1, [f3])]
    @test eltype(last(parts[1]).factors) == Fermion{Int}
    @test NonCommutativeProducts.with_autosort(() -> prod(last, parts), false) == x

    # Concrete factors form a single segment of the same type
    @test segments(NCMul(1, [f1, f2])) == [(Fermion{Int}, [f1, f2], Fermion{Int})]
    @test (@inferred split_by_group(NCMul(2, [f1, f2]))) == [Fermion{Int} => NCMul(2, [f1, f2])]
    @test segments(NCMul(3, Fermion{Int}[])) == []
    @test isempty(split_by_group(NCMul(3, Fermion{Int}[])))

    # A custom group with two types in it
    struct Up
        n::Int
    end
    struct Down
        n::Int
    end
    NonCommutativeProducts.@nc Up Down
    atom_group(::Union{Up,Down}) = :spin
    y = NCMul(1, Any[Up(1), Down(1), b, Up(2)])
    @test segments(y) == [(:spin, [Up(1), Down(1)], Any), (Boson, [b], Boson), (:spin, [Up(2)], Up)]
    @test split_by_group(y) == [:spin => NCMul(1, Any[Up(1), Down(1)]), Boson => NCMul(1, [b]), :spin => NCMul(1, [Up(2)])]

    # Tuple factors
    t = NCMul(2, (f1, b, f2))
    @test split_by_group(t) == [Fermion{Int} => NCMul(2, [f1]), Boson => NCMul(1, [b]), Fermion{Int} => NCMul(1, [f2])]
    @test collect(segment_ranges(t)) == [Fermion{Int} => 1:1, Boson => 2:2, Fermion{Int} => 3:3]
end

@testitem "A single segment does not copy" setup = [Fermions] begin
    import NonCommutativeProducts: NCMul, split_by_group, segment_ranges
    f1, f2 = Fermion(1), Fermion(2)

    x = NCMul(2, [f1, f2])
    @test last(only(split_by_group(x))) === x
    # The first segment tells whether there are others, without allocating
    has_one_segment(x) = length(last(first(segment_ranges(x)))) == length(x.factors)
    @test has_one_segment(x)
    @test (@allocated has_one_segment(x)) == 0

    # Non-concrete factors are narrowed
    a = NCMul(2, Any[f1, f2])
    @test last(only(split_by_group(a))) == x
    @test eltype(last(only(split_by_group(a))).factors) == Fermion{Int}
end

@testitem "Segment inputs and result types" begin
    import NonCommutativeProducts
    import NonCommutativeProducts: NCMul, NCAdd, foreach_segment, split_by_group, segment_ranges, with_autosort
    struct A
        n::Int
    end
    struct B
        n::Int
    end
    NonCommutativeProducts.@nc A B
    NonCommutativeProducts.@commutative A B
    NonCommutativeProducts.mul_effect(::A, ::A) = nothing
    NonCommutativeProducts.mul_effect(::B, ::B) = nothing
    segments(x) = (out = []; foreach_segment((g, seg) -> push!(out, g => collect(seg)), x); out)

    # A sum with one term counts as that term, also from sort
    for x in (with_autosort(() -> 2 * A(1) * B(1), true), sort(NCMul(2, Any[B(1), A(1)])))
        @test x isa NCAdd
        @test collect(segment_ranges(x)) == [A => 1:1, B => 2:2]
        @test segments(x) == [A => [A(1)], B => [B(1)]]
        @test split_by_group(x) == [A => NCMul(2, [A(1)]), B => NCMul(1, [B(1)])]
    end

    # A sum with several terms has no segments of its own
    y = NCMul(1, [A(1)]) + NCMul(1, [B(1)])
    @test_throws ArgumentError segment_ranges(y)
    @test_throws ArgumentError foreach_segment((g, seg) -> nothing, y)
    @test_throws ArgumentError split_by_group(y)

    # Only atoms have a group
    @test_throws ArgumentError NonCommutativeProducts.atom_group("not an atom")
    @test_throws ArgumentError NonCommutativeProducts.atom_group(NCMul(1, [A(1)]))

    # An atom is a product with one factor
    @test collect(segment_ranges(A(1))) == [A => 1:1]
    @test segments(A(1)) == [A => [A(1)]]
    @test split_by_group(A(1)) == [A => NCMul(1, [A(1)])]
    @test_throws ArgumentError segment_ranges("not an atom")

    # Sorting does not mutate the caller's factors, also when it swaps them
    v = Any[B(1), A(1)]
    v0 = copy(v)
    @test sort(NCMul(1, v)) == NCMul(1, [A(1), B(1)])
    @test v == v0

    # The parts of a view-backed product are Vector-backed, with or without factors
    pv = [A(1), A(2), A(3)]
    @test (@inferred split_by_group(NCMul(2, view(pv, 1:2)))) == [A => NCMul(2, [A(1), A(2)])]
    @test last(only(split_by_group(NCMul(2, view(pv, 1:2))))) isa NCMul{Int,A,Vector{A}}
    @test isempty(@inferred split_by_group(NCMul(2, view(pv, 1:0))))

    @test isempty(segment_ranges(NCMul(1, A[])))
    @test isempty(collect(segment_ranges(NCMul(1, A[]))))
end
