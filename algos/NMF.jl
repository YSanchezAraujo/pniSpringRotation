using LinearAlgebra, Distributions, DSP;

function matrix_distance(A::Array{Float64, 2}, B::Array{Float64, 2})
    UA, sA, VA = svd(A)
    UB, sB, VB = svd(B)
    return abs(sqrt(sum(sA.^2)) - sqrt(sum(sB.^2)))
end

function avgF2norm(A::Array{Float64, 2}, B::Array{Float64, 2})
    error = 0
    K = size(A, 1)

    for k in 1:K
	error += sqrt(sum((A[k, :] - B[k, :]).^2))
    end

    return error * 1/K
end

function F2norm(A::Array{Float64, 1}, B::Array{Float64, 1})
    return sqrt(sum((A - B).^2))
end

function mu_update(A, B, C)
    return A .* ((C*B') ./ (A*(B*B')) .+ eps())
end

"""
vanilla nmf with multiplicative updates
X: input data
k: rank of the approximating matrices
d: cost function
tol: tolerance for the change in cost
maxIter: maximum iteration to go up to
"""
function nmf(X::Array{Float64, 2}, 
	     k::Int64,
	     d::Function; tol=1.0e-5, maxIter=1.0e6)
    
    N, T = size(X)
    W = rand(Uniform(0.01, 0.05), (N, k))
    H = rand(Uniform(0.01, 0.05), (k, T))
    curIter, deltaCost, F = 1.0, 1.0, 0.0

    while (deltaCost > tol) && (curIter < maxIter)
	W = mu_update(W, H, X)
	H = mu_update(H', W', X')'
        X_tild = W*H
        F = d(X_tild, X)
        deltaCost = abs(deltaCost - F)
	curIter += 1
    end

    return W, H, F, deltaCost, curIter

end

"""
indep_norm2(W[1, :], H, X, 1)
gotta think about how to compute the gradient with respect to this
ForwardDiff.gradient(W_col -> indep_norm2(W[3, :], H, X, 3), rand(2))
"""
function indep_norm2(W_col, H, X, idx)
    X_i = X[idx, :]
    W_row, X_in = W_col', norm(X_i)
    return 0.5 * (W_row*(H*H')*W_col - 2*W_row*(H*X_i) + X_in)
end


function gradH(W, H, X)
    return W*(H*H') - X*H'
end

function gradW(W, H, X)
    return (W'*W)*H - W'*X
end
#gradloss(w0, H, X, idx) = ForwardDiff.gradient(W_col -> indep_norm2(

"""
vanilla nmf as above, but this time learning is done through gradient descent updates
X: input data
k: rank of the approximating matrices
d: cost function
tol: tolerance for the change in the cost
maxIter: maximum iteration to go up to
"""
function nmf_grad(X::Array{Float64, 2},
		  k::Int64,
		  d::Function; tol=1.0e-5, maxIter=1.0e6)

    N, T = size(X)
    W = rand(Uniform(0.01, 0.05), (N, k))
    H = rand(Uniform(0.01, 0.05), (k, T))
    curIter, deltaCost, F = 1.0, 1.0, 0.0

    while (deltaCost > tol) && (curIter < maxIter)
	# unlike the above, this is not vertoized by matrix operations
        # instead we do it explicitly iterating over each row

    end

    return W, H, F, deltaCost, curIter

end

"""
X: input data
k: number of factors
l: factor duration
λ: regularization strength
"""
function fit_seqnmf(X::Array{Float64, 2}, 
		    k::Int64, 
		    l::Float64, 
		    λ::Float64; tol=1.0e-5, maxIter=1.0e6)

    N, T = size(X)
    W = rand(Uniform(0.01, 0.05), (N, k))
    H = rand(Uniform(0.01, 0.05), (k, T))
    curIter, deltaCost = 1, 1

    while (curIter < maxIter) && (deltaCost > tol)
	# update H with MU
	# shift W and H to center W's in time
	# renormalize W and H so rows of H have unit norm
	# update W using MU update
	curIter += 1
    end



end