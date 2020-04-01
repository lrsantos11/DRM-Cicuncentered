using LinearAlgebra
using DelimitedFiles
using ProximalOperators
using Random
using Printf
using BenchmarkTools
using BenchmarkProfiles, Plots
using LaTeXStrings
using DataFrames
pgfplots()
# pyplot()
# PyPlot.rc("text", usetex= true)
# PyPlot.rc("text.latex",preamble = ["\\usepackage{amsmath,amssymb}", "\\usepackage[utf8]{inputenc}"])
include("mSet.jl")
include("BLSLP.jl/src/auxfunctions.jl")

global const EPS_VAL = 1e-6
global const ZERO_VAL = 1e-13
global const ITER_MAX = 1e5




    """
    ProjectBall(x, v, r)

    Returns the reflection of x over the L2-Ball with radius r an centered at v.

    Uses the ProximalOperators toolbox

    """

    function ProjectBall(x, v, r)
            Ball = IndBallL2(r)
            proj, fproj = prox(Translate(Ball,-v),x)
            return proj
    end


    """
    ReflectBall(x, inv, r)

    Returns the reflection of x over the L2-Ball with radius r an centered at v.

    Uses the ProximalOperators toolbox

    """
    function ReflectBall(x, v, r)
            Ball = IndBallL2(r)
            proj, fproj = prox(Translate(Ball,-v),x)
            return 2. * proj - x
    end

    """
    ReflectionDiagonal(X)

    Returns proj = y = z of R^n, the orthogonal projection of X = [x1;x2] in R^{2n} onto the two dimensional subspace of tuples
    such that y=z.
    """
    function ReflectionDiagonal(X::AbstractArray)
        n = length(X)
        if !iseven(n)
            @error "Dimension of x must be even!"
        end
        n = div(n,2)
        y = X[1:n]
        z = X[n+1:end]
        proj = .5*(y+z)
        reflec = 2. *[proj;proj] - X
        return reflec
    end

    """
    ProjectionDiagonal(X)

    Returns proj = y = z of R^n, the orthogonal projection of X = [x1;x2] in R^{2n} onto the two dimensional subspace of tuples
    such that y=z.
    """
    function ProjectionDiagonal(X::AbstractArray)
        n = length(X)
        if !iseven(n)
            @error "Dimension of x must be even!"
        end
        n = div(n,2)
        y = X[1:n]
        z = X[n+1:end]
        proj = .5*(y+z)
        return [proj;proj]
    end

    """
    Projection on the Product Space of Ball × Ball
    """
    function ProjectProductBall(x,v1,v2,r1,r2)
        n = length(x)
        n = div(n,2)
        y = ProjectBall(x[1:n],v1, r1)
        z = ProjectBall(x[1:n],v2, r2)
        return [y; z]
    end

    """
    Reflection on the Product Space of Ball × Ball
    """
    function ReflectProductBall(x,v1,v2,r1,r2)
        n = length(x)
        n = div(n,2)
        y = ReflectBall(x[1:n],v1, r1)
        z = ReflectBall(x[1:n],v2, r2)
        return [y; z]
    end

    """
    Projection using Indicator Function form Proximal Operators
    """
    function ProjectIndicator(indicator,x)
        proj, fproj = prox(indicator,x)
        return proj
    end

    """
    Reflection using Indicator Function form Proximal Operators
    """
    function ReflectIndicator(indicator,x)
        proj, fproj = prox(indicator,x)
        reflec = 2*proj - x
        return reflec
    end


    # """
    # Defines ReflectA and ReflectB
    # """
    # For the Two Balls Problem
    # v1 = [.8,0]
    # r1 = 1.
    # v2 = [-.8,0]
    # r2 = 1.
    # ReflectA(x) =  ReflectProductBall(x,v1,v2,r1,r2)
    # ReflectB(x) =  ReflectionDiagonal(x)
    # ProjectA(x) =  ProjectProductBall(x,v1,v2,r1,r2)
    # ProjectB(x) =  ProjectionDiagonal(x)
    # ReflectA(x) =  ReflectBall(x,v1,r1)
    # ReflectB(x) =  ReflectBall(x,v2,r2)
    # ProjectA(x) =  ProjectBall(x,v1, r1)
    # ProjectB(x) =  ProjectBall(x,v2, r2)

    # For the Three Balls Problem
    # v1 = [.8,0]
    # r1 = 1.
    # v2 = [-.8,0]
    # r2 = 1.
    # ReflectA(x) =  ReflectProductBall(x,v1,v2,r1,r2)
    # ReflectB(x) =  ReflectionDiagonal(x)
    # ProjectA(x) =  ProjectProductBall(x,v1,v2,r1,r2)
    # ProjectB(x) =  ProjectionDiagonal(x)
    # ReflectA(x) =  ReflectBall(x,v1,r1)
    # ReflectB(x) =  ReflectBall(x,v2,r2)
    # ProjectA(x) =  ProjectBall(x,v1, r1)
    # ProjectB(x) =  ProjectBall(x,v2, r2)



# TestConvex(xzero,200)

function TestAffineSOC(;n::Int64 = 100,samples::Int64 = 2,
        EPS_VAL::Float64 = 1e-6,printfile::Bool=false,itmax::Int64 = 200, restarts::Int64 = 1)
    # Fix Random
    Random.seed!(10)
    # X = R^n
    j = 1
    firstwrite = true
    resultsDRM = Array{Float64}(undef,0,3)
    resultsCRM = Array{Float64}(undef,0,3)
    resultsMAP = Array{Float64}(undef,0,3)
    SOC = IndSOC()
    while (j <=samples)
        #########################
        # Generate Subspace
        # affine = true
        # cmax = rand(1:ceil(Integer,n/10,))
        # CC  = GenerateSamples(n,affine)
        # #  Read Files
        m = rand(1:n-1)
        A =  randn(m,n)
        b = randn(m)
        w  = (A\b)[2:end]
        b = A*[norm(w); w]
        Affine = IndAffine(A,b)
        ReflectA(x) =  ReflectIndicator(SOC,x)
        ReflectB(x) =  ReflectIndicator(Affine,x)
        ProjectA(x) =  ProjectIndicator(SOC,x)
        ProjectB(x) =  ProjectIndicator(Affine,x)
        # Restarts
        for i = 1:restarts
            xzero = StartingPoint(n)
            xzero = ProjectB(xzero)
            # @show SOC(xzero)
            while SOC(xzero) != Inf
                xzero = StartingPoint(n)
                xzero = ProjectB(xzero)
            end
            prob_name  = String("Problem$j"*"A$m"*"×$n")
            println(prob_name)
            resultDRM, resultCRM,resultMAP = TestConvex(xzero,itmax,ReflectA,ReflectB,ProjectA,ProjectB)
            # println("CRM")
            if resultDRM.iterations <= resultCRM.iterations
                println("--")
                println("--")
                println("--")
                kCRM = resultCRM.iterations
                kDRM = resultDRM.iterations
                println("CRM: $kCRM")
                println("DRM: $kDRM")
                println("--")
                println("--")
                println("--")
            end
            if printfile
                fname = @sprintf("tables/AffineSOC-N%d_Samp%d_%.1E-CRM.dat",n,samples,EPS_VAL);
                open(fname,"a") do file
                    if firstwrite
                        writedlm(file,"---\nalgname: CRM\nsuccess: conv\nfree_format: True\n---\n")
                        firstwrite = false
                    end
                    write(file,prob_name)
                    str = @sprintf(" %10.8e %10.8f %10.8f",resultCRM.residual, resultCRM.iterations,resultCRM.elapsedtime)
                    writedlm(file,str)
                end
            else
                resultsCRM = vcat(resultsCRM,[resultCRM.residual  resultCRM.iterations resultCRM.elapsedtime])
            end
            # println("MAP")
            if printfile
                fname = @sprintf("tables/AffineSOC-N%d_Samp%d_%.1E-MAP.dat",n,samples,EPS_VAL);
                open(fname,"a") do file
                    if firstwrite
                        writedlm(file,"---\nalgname: MAP\nsuccess: conv\nfree_format: True\n---\n")
                        firstwrite = false
                    end
                    write(file,prob_name)
                    str = @sprintf(" %10.8e %10.8f %10.8f",resultMAP.residual, resultMAP.iterations,resultMAP.elapsedtime)
                    writedlm(file,str)
                end
            else
                resultsMAP = vcat(resultsMAP,[resultMAP.residual  resultMAP.iterations resultMAP.elapsedtime])
            end
            # println("DRM")
            if printfile
                fname = @sprintf("tables/AffineSOC-N%d_Samp%d_%.1E-DRM.dat",n,samples,EPS_VAL);
                open(fname,"a") do file
                    if firstwrite
                        writedlm(file,"---\nalgname: DRM\nsuccess: conv\nfree_format: True\n---\n")
                        firstwrite = false
                    end
                    write(file,prob_name)
                    str = @sprintf(" %10.8e %10.8f %10.8f",resultDRM.residual, resultDRM.iterations,resultDRM.elapsedtime)
                    writedlm(file,str)
                end
            else
                resultsDRM = vcat(resultsDRM,[resultDRM.residual  resultDRM.iterations resultDRM.elapsedtime])
            end
        end
        j += 1
    end
    return resultsCRM, resultsMAP, resultsDRM
end



n = 200
samples = 100
resultsCRM, resultsMAP, resultsDRM = TestAffineSOC(n = n, samples = samples,itmax=5000, EPS_VAL = 1e-6, restarts = 10)
# # plot(s1)
perprof = performance_profile(hcat(resultsCRM[:,2],resultsDRM[:,2],resultsMAP[:,2]), ["CRM", "DRM", "MAP"],
    title=L"Performance Profile -- Gap error -- $\varepsilon = 10^{-6}$",
    legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash, :dot])
ylabel!("Percentage of problems solved", 
    )
savefig(perprof, "../../Draft/New/Two-Any-Convex-Sets/figures/AffineSoc-N$n"*"_Samp$samples"*"-bw.pdf")

df = DataFrame()
df.CRM = resultsCRM[:,2]
df.DRM = resultsDRM[:,2]
df.MAP = resultsMAP[:,2]
@show describe(df)# 
#@show norm(resultsCRM.solution - resultsMAP.solution)

 # TestAffineSOC(n = 200, samples = 1,itmax=5000, EPS_VAL = 1e-6, restarts = 1)



"""
Projection on the Product Space of Half Spaces
"""
function ProjectProductHSP(xpto,halfspacesInd)
    proj = []
    for halfspace in halfspacesInd
        push!(proj,ProjectIndicator(halfspace,xpto))
    end
    return proj
end


"""
Reflection on the Product Space of Half-spaces
"""
function ReflectProductHSP(xpto,halfspacesInd)
    reflec = []
    for halfspace in halfspacesInd
        push!(reflec,ReflectIndicator(halfspace,xpto))
    end
    return reflec
end


"""
ProjectDiagonal(X)


"""
function ProjectDiagonal(X::AbstractArray)
    num_sets = length(X)
    n = length(X[1])
    sum_vectors = zeros(n)
    for index in eachindex(X)
        sum_vectors +=X[index]
    end
    sum_vectors ./= num_sets
    proj = []
    for index in 1:num_sets
        push!(proj,sum_vectors)
    end
    return proj
end

"""
ReflectDiagonal(X)

"""
function ReflectDiagonal(X::AbstractArray)
    proj = ProjectDiagonal(X)
    return 2*proj  - X
end


"""
MAPmSet(X)

"""
function MAPmSet(halfspacesInd,xzero,itmax::Int64 = 2000;
                            printfile::Bool=true,EPS_VAL::Float64=1e-6, )
                           # Projecting onto A first
    try
         rm("tables/xMAPmSet.dat")
    catch
    end
    tolMAPmset = 1.
    k = 0
    xMAP = xzero
    m = length(halfspacesInd)
    while tolMAPmset > EPS_VAL   && k <= itmax
        # xbar = xstar #Looking for FixingPoint
        xMAPold = xMAP
        for halfspace in halfspacesInd
            xMAP = ProjectIndicator(halfspace,xMAP)
        end
        # tolMAPmset = norm(xMAP - ProjectIndicator(halfspacesInd[1],xMAP)) #Gap Error
        tolMAPmset = norm(xMAP  - xMAPold) #Fixed point Error
        println(tolMAPmset)
        k+=1
        if printfile
            open("tables/xMAPmSet.dat","a") do file
            # do stuff with the open file
              writedlm(file,[k tolMAPmset])
            end
        end
    end
    # println("xDR: $xDR")
    println("Projections of  MAPmset: $k")
    println("Norm for xMAPmSet: $tolMAPmset")
    resultMAPmset = Results(tolMAPmset,k,0.,xMAP)
    return resultMAPmset
end




"""
TestPolyhedral(;n::Int64 = 100,samples::Int64 = 2,
        EPS_VAL::Float64 = 1e-5,printfile::Bool=false,itmax::Int64 = 200)

"""
function TestPolyhedral(;n::Int64 = 100,samples::Int64 = 1,
        EPS_VAL::Float64 = 1e-5,itmax::Int64 = 1000, restarts = 1)
    # X = R^n
    firstwrite = true
    resultsDRM = Array{Float64}(undef,0,3)
    resultsCRM = Array{Float64}(undef,0,3)
    resultsMAP = Array{Float64}(undef,0,3)
    resultsMAPmset = Array{Float64}(undef,0,3)

    # Fix Random
    Random.seed!(1)
    for j in 1:samples
        # n can be graeater
        n = rand(n:2*n)
        A  = GenerateSamples(n,false)[1]
        @show m, n = size(A)
        xbar = StartingPoint(n)
        bbar = A*xbar
        # Number of semispaces
        non_slater = rand(1:m)
        indices = sort(rand(1:m,non_slater))
        b = bbar
        b[indices] = bbar[indices] + norm(bbar[indices])*rand(non_slater)
        halfspacesInd = []
        for row  in eachrow([A b])
            push!(halfspacesInd, IndHalfspace(row[1:n],row[end]))
        end



        ReflectA(x) =  ReflectProductHSP(x[1],halfspacesInd)
        ReflectB(x) =  ReflectDiagonal(x)
        ProjectA(x) =  ProjectProductHSP(x[1],halfspacesInd)
        ProjectB(x) =  ProjectDiagonal(x)
        for i = 1:restarts
            xstart = StartingPoint(n)
            xzero = []
            for index in 1:m
                push!(xzero,xstart)
            end

            xzero = ProjectB(xzero)

            xMAPmSet =  MAPmSet(halfspacesInd,xzero[1],itmax)

            resultDRM, resultCRM, resultMAP = TestConvex(xzero,itmax,ReflectA,ReflectB,ProjectA,ProjectB)
            resultsDRM = vcat(resultsDRM,[resultDRM.residual  resultDRM.iterations resultDRM.elapsedtime])
            resultsCRM = vcat(resultsCRM,[resultCRM.residual  resultCRM.iterations resultCRM.elapsedtime])
            resultsMAP = vcat(resultsMAP,[resultMAP.residual  resultMAP.iterations resultMAP.elapsedtime])
            resultsMAPmset = vcat(resultsMAPmset,[xMAPmSet.residual  xMAPmSet.iterations xMAPmSet.elapsedtime])
            # println("true error CRM")
        # println((b - A*resultCRM.solution))
        # println("true error DRM")
        # println((b - A*resultDRM.solution))
        # println("true error MAP")
        # println((b - A*resultMAP.solution))
    end
    end
    return resultsDRM, resultsCRM, resultsMAP, resultsMAPmset
end
# n = 200
# samples = 1
# restarts = 1
# resultsDRM, resultsCRM, resultsMAP, xMAPmset = TestPolyhedral(n = n, samples = samples,itmax=5000, restarts = restarts)
# xCRM = readdlm("tables/xCRM.dat")
# xDRM = readdlm("tables/xDR.dat")
# xMAP = readdlm("tables/xMAP.dat")
# xMAPmset = readdlm("tables/xMAPmSet.dat")
# p1 = plot(xCRM[:,1],xCRM[:,2],scale=:log10, label="CRM-prod",
#            title="Comparison using Product Space reformulation",
#            framestyle = :box,
#            xlabel = "Number of iterations (log scale)",
#            ylabel = "Gap error (log scale)",
#            minorticks=true);
# p1 = plot!(xDRM[:,1],xDRM[:,2],scale=:log10, label="DRM-prod",linestyle=:dash);
# p1 = plot!(xMAP[:,1],xMAP[:,2],scale=:log10, label="MAP-prod",linestyle=:dot); # linestyles=[:solid, :dash, :dot])
# # plot(p1)
# # performance_profile(hcat(resultsCRM[:,2],resultsDRM[:,2],resultsMAP[:,2]), ["CRM", "DRM", "MAP"], title="Polyhedral Intersection")
# # savefig(p1, "tables/Polyhedral"*string(now())*".pdf")
# savefig(p1, "../../Draft/New/Two-Any-Convex-Sets/figures/ComparisonCRMProd_DRM-prod_MAP-prod-bw.pdf")
# df = DataFrame()
# df.CRM = resultsCRM[:,2]
# df.DRM = resultsDRM[:,2]
# df.MAP = resultsMAP[:,2]
# describe(df)






#
# mpsfile = "BLSLP.jl/mps/25fv47.mps"
# c, A, bl, bu, xl, xu = mpstomatrix(mpsfile)
# m,n,c,A,bu,xl,xu = LPtoSTDFormat(c,A,bl,bu,xl,xu)

# Random.seed!(10)
# # n = 100
# # n = rand(n:2*n)
# # A  = GenerateSamples(n,false)[1]
# m, n = size(A)
# xbar = StartingPoint(n+m)
# # bbar = [A I]*xbar
# # Number of semispaces
# # non_slater = rand(1:m)
# # indices = sort(rand(1:m,non_slater))
# # b = bbar
# # b[indices] = bbar[indices] + norm(bbar[indices])*rand(non_slater)
# Abar = [A I]
# Affine = IndAffine(A, bu)
# NonnegativeCone = IndNonnegative()
#
# ReflectA(x) =  ReflectIndicator(NonnegativeCone, x)
# ReflectB(x) =  ReflectIndicator(Affine, x)
# ProjectA(x) =   ProjectIndicator(NonnegativeCone, x)
# ProjectB(x) =  ProjectIndicator(Affine, x)
#
# xzero = ProjectB(sparse(StartingPoint(n)))
# itmax = 1000
#
# resultDRM, resultCRM, resultMAP = TestConvex(xzero,itmax,ReflectA,ReflectB,ProjectA,ProjectB);
