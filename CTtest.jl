using MatrixMarket
using SparseArrays
using LinearAlgebra
# using ImageView
using Images
using MAT


include("mSet.jl")


filename = "matrices/phantom.mat"
mf = matread(filename)
A = sparse(mf["A"])
b = (mf["b"])
sol = (mf["x"])
Xkacz = (mf["Xkacz"])
Xsymk = (mf["Xsymk"])
Xrand = (mf["Xrand"])



# Look for zero rows
J = Array{Int64}(undef,0)
for row in eachrow(A)
       if ~iszero(row)
         push!(J,row.indices[1])
       end
end

M = A[J,:]
bb = b[J]
num_row, num_col = size(M)
blocksize = [1,16, 64, 256]
xzero = zeros(num_col)
X = []
println("==================")
println("Results for Bw-CRM - Matrix size: $num_row×$num_col")
println("Method  |  ",   " norm Ax-b | ", " norm x-sol | ")
for index in blocksize
    Blocks = formBlocksbySize(num_row,index)
    # elapsed_time =  0. #@elapsed
    xstar, iter, reflec, normax =  BwCRM_HyperplaneClosedForm(M,bb,Blocks,xstar = xzero,iter_max = 10)
    dist =  norm(xstar - sol)
    # elapsed_time =  @belapsed BwCRM23($M,$b,$Blocks,xstar = $xzero, EPS_VAL = $tol)
    @printf("Bw-CRM-%d     %d      %.6e          %.6e      \n", index, iter,  normax, dist)
    push!(X,xstar)
end
println("==================")
