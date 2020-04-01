# workspace()
using Plots
# using JLD
    # using HDF5
# using BenchmarkProfiles
using BenchmarkTools
using MatrixMarket
using SparseArrays
using LinearAlgebra
using ProximalOperators
using DelimitedFiles
# pgfplots()
# using LaTeXStrings
# PGFPlots.pushPGFPlotsPreamble("\\usepackage{amssymb}")
# unicodeplots()

# include("CRM.jl")
include("mSet.jl")
global const EPS_VAL = 1e-6
global const ZERO_VAL = 1e-12
global const ITER_MAX = 1e5




"""
Projection using Indicator Function form Proximal Operators
"""
function ProjectIndicator(indicator,x;γ::Float64 = 1.0)
    proj, fproj = prox(indicator,x,γ)
    return proj
end

"""
Reflection using Indicator Function form Proximal Operators
"""
function ReflectIndicator(indicator,x; γ::Float64 = 1.0)
    proj, fproj = prox(indicator,x,γ)
    reflec = 2*proj - x
    return reflec
end



function BasisPursuit(A::AbstractArray{Float64,2},b::AbstractArray{Float64,1};
     EPS_VAL = 1e-6, xzero=[],sol=[], itmax = ITER_MAX, restarts::Int64 =  1)

    num_row,num_col = size(A)
    # X = R^n
    j = 1
    firstwrite = true
    resultsDRM = Array{Float64}(undef,0,3)
    resultsCRM = Array{Float64}(undef,0,3)
    resultsMAP = Array{Float64}(undef,0,3)
    resultDRM = []
    resultCRM = []
    resultMAP = []
    Affine = IndAffine(A,b)
    Norm1 = IndBallL1(1.)
    ReflectA(x) =  ReflectIndicator(Norm1,x)
    ReflectB(x) =  ReflectIndicator(Affine,x)
    ProjectA(x) =  ProjectIndicator(Norm1,x)
    ProjectB(x) =  ProjectIndicator(Affine,x)
    Random.seed!(10)
    # Restarts

    for i = 1:restarts
        # If no initial point, compute a random one
        if isempty(xzero)
            println("xzero")
            xzero = StartingPoint(num_col)
        end
        xzero = ProjectB(xzero)
        resultDRM, resultCRM, resultMAP = TestConvex(xzero,itmax,ReflectA,ReflectB,ProjectA,ProjectB)
        resultsCRM = vcat(resultsCRM,[resultCRM.residual  resultCRM.iterations resultCRM.elapsedtime])
        resultsDRM = vcat(resultsDRM,[resultDRM.residual  resultDRM.iterations resultDRM.elapsedtime])
        resultsMAP = vcat(resultsMAP,[resultMAP.residual  resultMAP.iterations resultCRM.elapsedtime])
    end
     return resultDRM, resultCRM, resultMAP
end
####################################

A = [.5 1]
b = [1.]
xstar = [3,3.]
resultDRM, resultCRM, resultMAP = BasisPursuit(A,b,xzero = xstar,itmax = 300)
