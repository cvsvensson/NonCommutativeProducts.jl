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
consecutive factors in the same group, [`segment_ranges`](@ref NonCommutativeProducts.segment_ranges) returns the group
and the index range of each segment, and [`split_by_group`](@ref NonCommutativeProducts.split_by_group) returns the
segments as products. This is useful e.g. to act with each part of a product on its own space. Groups do not change
the order of the factors, so interleaved groups give several runs, and a group appears once for each run of its
factors:
```jldoctest multiple
julia> import NonCommutativeProducts: NCMul, split_by_group, foreach_segment, segment_ranges, atom_group

julia> x = NCMul(2, Any[Foo(), Foo(), Bar(), Foo()]);

julia> split_by_group(x)
3-element Vector{NCMul{Int64}}:
 2*Foo()*Foo()
 Bar()
 Foo()

julia> segment_ranges(x)
3-element Vector{Pair{DataType, UnitRange{Int64}}}:
 Foo => 1:2
 Bar => 3:3
 Foo => 4:4

julia> foreach_segment((g, seg) -> println(g, ": ", seg, ", element type ", eltype(seg)), x)
Foo: [Foo(), Foo()], element type Foo
Bar: [Bar()], element type Bar
Foo: [Foo()], element type Foo
```
The factors of each segment have a concrete element type when they all have the same type, also when the factors of
the product do not, and the function given to `foreach_segment` is compiled for each type of segment. If the factors
of the product have a concrete element type, the segments are views of them. The factors of a product are never
mutated, by the convention of this package, so the views stay valid after the function returns; copy a segment only
if you intend to mutate it. `split_by_group` is inferable for factors of a concrete element type; otherwise its result
has an abstract element type, and `foreach_segment` is the better choice in hot paths. The first part of
`split_by_group` gets the coefficient of the product, also when it is zero, and the other parts get coefficient one.
Sums have no segments: apply these functions to each of their `terms`.

Define `atom_group` to group atoms by value, e.g. by the space they act on, or to put several types into one group:
```jldoctest multiple
julia> atom_group(::Union{Foo,Bar}) = :foobar;

julia> split_by_group(x)
1-element Vector{NCMul{Int64, Any, Vector{Any}}}:
 2*Foo()*Foo()*Bar()*Foo()
```
A group can also encode a kind of atom. For example, a downstream package can use `(space_id, :ket)` and
`(space_id, :op)` as groups, to keep the kets and the operators of the same subsystem in separate segments.

## Commuting groups

Sorting a product of several groups calls `mul_effect` on every pair of neighbouring factors, also to move factors of
one group past those of another. If the factors of two groups commute, or commute up to a sign, declare it with
[`group_commutation`](@ref NonCommutativeProducts.group_commutation), which returns `:commute`, `:graded` or `nothing`
(unknown, the default). For `:graded` groups, two atoms anticommute if both are odd, as given by
[`atom_parity`](@ref NonCommutativeProducts.atom_parity). Sorting then first moves the factors to the order of their
groups, with the sign of the odd atoms that cross, and sorts each group on its own, with compiled `mul_effect` calls
when the factors of the group have the same type:
```jldoctest groups
julia> import NonCommutativeProducts as NCP

julia> struct Up
           n::Int
       end

julia> struct Down
           n::Int
       end

julia> NCP.@nc Up Down

julia> NCP.mul_effect(a::T, b::T) where {T<:Union{Up,Down}} = a.n == b.n ? 0 : a.n > b.n ? NCP.Swap(-1) : nothing

julia> NCP.mul_effect(::Up, ::Down) = nothing

julia> NCP.mul_effect(::Down, ::Up) = NCP.Swap(-1)

julia> NCP.group_commutation(::Type{Up}, ::Type{Down}) = :graded;

julia> NCP.group_commutation(::Type{Down}, ::Type{Up}) = :graded;

julia> NCP.atom_parity(::Union{Up,Down}) = true;

julia> sort(NCP.NCMul(1, Any[Down(1), Up(2), Down(2), Up(1)]))
Up(1)*Up(2)*Down(1)*Down(2)
```
The declaration is a promise about `mul_effect`, which must still be defined between the groups: it must order the
factors of the two groups by group alone, with `Swap(1)`, or `Swap(-1)` for two odd atoms of graded groups, as above
and as [`@commutative`](@ref NonCommutativeProducts.@commutative) does. The order of the groups is taken from
`mul_effect`, or from [`group_isless`](@ref NonCommutativeProducts.group_isless) if it is defined, which must agree with
it. Sorting then gives exactly the same result as without the declaration, including the order of the factors, and the
declaration only makes it faster. The hooks are opt-in: without them, nothing changes.

The group order is used for products whose factors have an abstract element type, have at least 16 factors, and whose
groups all have a declared commutation with each other. Other products, e.g. those with a pair of groups whose
commutation is `nothing`, are sorted by `mul_effect` alone. The results of `group_commutation` and `group_isless` for
groups that are types are cached until methods are next defined, so they must depend only on their arguments.

## Performance of mixed products

Products whose factors have a concrete type, such as a single atom type, sort fastest, since each `mul_effect` call
is resolved when the sorting code is compiled. The same holds for a `Union` of two concrete types. Factors of several
types are otherwise stored in a vector with an abstract element type, e.g. `Vector{Any}` or a vector of a common
abstract supertype of the atoms, and each `mul_effect` is then found at runtime. To work on the factors of one type at
a time with compiled code, use [`foreach_segment`](@ref NonCommutativeProducts.foreach_segment). To sort long products of
groups that commute faster, declare their commutation, see [Commuting groups](@ref).
