# Multiple types

Several types can be used in the same expressions. When sorting, [`mul_effect`](@ref NonCommutativeProducts.mul_effect) is called on every pair of neighbouring factors, so each pair of types that can end up next to each other needs one. For types that commute with each other, the macro [`@commutative`](@ref NonCommutativeProducts.@commutative) defines mixed arithmetic and the mixed `mul_effect`s for you. Factors are then ordered by type, in the order the types are listed.
```jldoctest multiple
julia> import NonCommutativeProducts: @commutative, @nc

julia> struct Foo end

julia> struct Bar end

julia> @nc Foo Bar

julia> @commutative Foo Bar

julia> 2Bar()Foo() - Foo()Bar() |> sort
Foo()*Bar()
```
For types that don't commute, define `mul_effect(a::Foo, b::Bar)` and `mul_effect(a::Bar, b::Foo)` yourself.

## Groups and segments

Each atom belongs to a group, given by [`atom_group`](@ref NonCommutativeProducts.atom_group). The default group is
the type of the atom, so parametric variants such as `Fermion{Int}` and `Fermion{Int32}` are different groups.
[`foreach_segment`](@ref NonCommutativeProducts.foreach_segment) splits a product into segments, the maximal runs of
consecutive factors in the same group, [`segment_ranges`](@ref NonCommutativeProducts.segment_ranges) iterates over the
group and the index range of each segment, and [`split_by_group`](@ref NonCommutativeProducts.split_by_group) returns
the group and the factors of each segment as a product. This is useful e.g. to act with each part of a product on its
own space. Groups do not change the order of the factors, so a group can have several segments:
```jldoctest multiple
julia> import NonCommutativeProducts: NCMul, split_by_group, foreach_segment, segment_ranges, atom_group

julia> x = NCMul(2, Any[Foo(), Foo(), Bar(), Foo()]);

julia> show_parts(x) = foreach(((g, part),) -> println(g, " => ", part), split_by_group(x));

julia> show_parts(x)
Foo => 2*Foo()*Foo()
Bar => Bar()
Foo => Foo()

julia> collect(segment_ranges(x))
3-element Vector{Pair{Any, UnitRange{Int64}}}:
 Foo => 1:2
 Bar => 3:3
 Foo => 4:4

julia> foreach_segment((g, seg) -> println(g, ": ", seg, ", element type ", eltype(seg)), x)
Foo: [Foo(), Foo()], element type Foo
Bar: [Bar()], element type Bar
Foo: [Foo()], element type Foo
```
Define `atom_group` to group atoms by value, e.g. by the space they act on, or to put several types into one group:
```jldoctest multiple
julia> atom_group(::Union{Foo,Bar}) = :foobar;

julia> show_parts(x)
foobar => 2*Foo()*Foo()*Bar()*Foo()
```

## Performance of mixed products

Factors of a concrete type sort fastest. For abstract element types, use
[`foreach_segment`](@ref NonCommutativeProducts.foreach_segment) to work on each type with compiled code.
