# workspace()
# using Plots
# pgfplots()
# using LaTeXStrings
# PGFPlots.pushPGFPlotsPreamble("\\usepackage{amssymb}")
# unicodeplots()
using MatrixMarket
using SparseArrays
using LinearAlgebra
# using BenchmarkProfiles
using BenchmarkTools
include("mSet.jl")
# global const EPS_VAL = 1e-6
global const ZERO_VAL = 1e-13
# global const ITER_MAX = 1e6


FIDAP = MatrixMarket.mmread("matrices/fidap005.mtx")
# FIDAP = MatrixMarket.mmread("matrices/orsirr_1.mtx")
m,n = size(FIDAP)
spM= findnz(FIDAP)
nnzFIDAP = nnz(FIDAP)
Random.seed!(2)
M = sparse(spM[1],spM[2],randn(nnzFIDAP))
xzero = zeros(n)
numrows = 27
b = ones(numrows)
M = M[1:numrows,:]
numblocks = numrows
blocksize = [1,3,9,27]
BW = []
results=[]
tol = 1e-5
println("==================")
println("Results for Bw-CRM - Matrix size: $numrows×$n")
println("Method  |  ",  " #Blocks | ",  " #Projections | ",   " #Iterations | ", " norm Ax-b | ", " Time")
for index in blocksize
    Blocks = formBlocksbySize(numblocks,index)
    elapsed_time =  0. #@elapsed
    xstar, iter, reflec, normax =  BwCRM_HyperplaneClosedForm(M,b,Blocks,xstar = xzero,EPS_VAL = tol)
    elapsed_time =  @belapsed BwCRM_HyperplaneClosedForm($M,$b,$Blocks,xstar = $xzero, EPS_VAL = $tol)
    @printf("Bw-CRM-%d       %5d       %5d       %5d        %.4e      %.4e \n", index, div(numrows,index),reflec, iter, normax, elapsed_time)
end
println("==================")
# println("Results for MAP")
# results = MAPmSets(M,b,xstar = xzero)
# println("MAP - Projections ",results[3])
