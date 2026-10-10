# Precompile workload. Atom types and their mul_effect rules are defined downstream, so the only code this package
# can cache for them is the type-independent path taken by mixed-type expressions: products with
# factors::Vector{Any}, sums keyed on them, and the sorting machinery around mul_effect. Every operation on a
# factor of such a product is a runtime dispatch, so the compiled code doesn't depend on the atom type, and plain
# integers can stand in for atoms here. That way the workload defines no types and no methods.
#
# Sorting is run on products of at most one factor, which never call mul_effect but still compile the sorting loop
# (see _mul_effect and _apply_effect!! for why that loop isn't invalidated by downstream methods). What the loop
# does with the result of mul_effect is compiled by calling _apply_effect!! directly with the effects that don't
# involve downstream types: numbers, Swap, AddTerms, and products and sums of Vector{Any} factors.

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
        # sorting, but only of terms with at most one factor, see above
        sort(x), sort(c * anymul(1, 1) + anymul(1, 2) + 1)
    end
    function effect_workload(c, effect)
        ncmul = anymul(c, 1, 2, 3)
        terms, _, _ = _apply_effect!!([ncmul], 1, ncmul, 1, effect)
        _sum_sorted_terms(terms, typeof(c))
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
