module NonCommutativeProductsKrylovKitExt

using NonCommutativeProducts
import KrylovKit

const NC = NonCommutativeProducts

# Treat NC expressions as linear operators acting by left multiplication. Results are sorted because equal operators
# left as separate terms make a zero residual look nonzero (their cancellation in inner leaves roundoff of order eps,
# not eps^2), so solvers miss the invariant subspace and normalize noise.
KrylovKit.apply(operator::NC.MulAdd, x) = NC.with_autosort(() -> operator * x)
KrylovKit.apply_normal(operator::NC.MulAdd, x) = NC.with_autosort(() -> operator * x)
KrylovKit.apply_adjoint(operator::NC.MulAdd, x) = NC.with_autosort(() -> operator' * x)

end
