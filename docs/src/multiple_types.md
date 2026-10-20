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
