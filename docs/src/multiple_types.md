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

Each atom belongs to a [`group`](@ref NonCommutativeProducts.group), which is its type unless you define a method for
your types. [`foreach_segment`](@ref NonCommutativeProducts.foreach_segment) splits a product into segments, the
maximal runs of consecutive factors in the same group, and
[`split_by_group`](@ref NonCommutativeProducts.split_by_group) returns the segments as products. This is useful e.g.
to act with each part of a product on its own space. Groups do not change the order of the factors, so a group appears
once for each run of its factors:
```jldoctest multiple
julia> import NonCommutativeProducts: NCMul, split_by_group, foreach_segment, group

julia> x = NCMul(2, Any[Foo(), Foo(), Bar(), Foo()]);

julia> split_by_group(x)
3-element Vector{NCMul{Int64}}:
 2*Foo()*Foo()
 Bar()
 Foo()

julia> foreach_segment((g, seg) -> println(g, ": ", seg, ", element type ", eltype(seg)), x)
Foo: [Foo(), Foo()], element type Foo
Bar: [Bar()], element type Bar
Foo: [Foo()], element type Foo
```
The factors of each segment have a concrete element type when they all have the same type, also when the factors of
the product do not, and the function given to `foreach_segment` is compiled for each type of segment.

Define `group` to group atoms by value, e.g. by the space they act on, or to put several types into one group:
```jldoctest multiple
julia> group(::Union{Foo,Bar}) = :foobar;

julia> split_by_group(x)
1-element Vector{NCMul{Int64, Any, Vector{Any}}}:
 2*Foo()*Foo()*Bar()*Foo()
```

## Performance of mixed products

Products whose factors have a concrete type, such as a single atom type, sort fastest, since each `mul_effect` call
is resolved when the sorting code is compiled. The same holds for a `Union` of two concrete types. Factors of several
types are otherwise stored in a vector with an abstract element type, e.g. `Vector{Any}` or a vector of a common
abstract supertype of the atoms, and each `mul_effect` is then found at runtime. To work on the factors of one type at
a time with compiled code, use [`foreach_segment`](@ref NonCommutativeProducts.foreach_segment).
