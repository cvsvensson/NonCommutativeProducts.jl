# Atom types are defined downstream, so the workload covers the type-independent path of mixed-type expressions:
# products with factors::Vector{Any} and sums of them. Operations on such factors dispatch at runtime, so integers
# stand in for atoms and the workload defines no types or methods.
# Sorting runs only on products with at most one factor, which compiles the sorting loop without calling mul_effect.
# The effects that involve no downstream types are applied to a product directly, and the results are added to a sum.

@setup_workload begin
    anymul(c, xs...) = NCMul(c, Any[xs...])
    function arithmetic_workload(c)
        x, y, z = anymul(c, 1), anymul(1, 2), anymul(1, 3, 4)
        p = x * y * z
        h = sum(c * anymul(1, n) * anymul(1, n + 1) + anymul(1, n) for n in 1:2)
        h2 = h^2 + adjoint(h) - 1 + p' - p * h
        q = c * p + h2 / 2 + 2 * LinearAlgebra.I - x + (y + z)
        foreach((p, q, h2, x)) do e
            sprint(show, e)
            show(IOContext(IOBuffer(), :limit => true), MIME"text/plain"(), e)
        end
        h == h2, p == q, p == x, hash(p), hash(h), iszero(h), isscalar(h)
        anyadd(h)
        # only terms with at most one factor, see above
        sort(x), sort(c * anymul(1, 1) + anymul(1, 2) + 1)
    end
    function effect_workload(c, effect)
        pending = _push_part!!(nothing, c, Any[1, 2, 3], 1, effect)
        acc = zero(NCAdd{typeof(c),NCMul{Int,Any,Vector{Any}}})
        isnothing(pending) && return acc
        return fold_retyping!!((acc, p) -> add_term!!(acc, prefactor(p), p.factors), acc, pending)
    end
    @compile_workload begin
        with(_autosort => false) do
            foreach((1, 1.0, 1.0im)) do c
                arithmetic_workload(c)
                foreach((0, 1, c, Swap(-1), Swap(1), AddTerms((Swap(-1), 1)), anymul(-1, 2, 1), -anymul(1, 1) + 1)) do effect
                    effect_workload(c, effect)
                end
            end
        end
    end
end
