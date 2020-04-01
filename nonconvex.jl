include("MAP.jl")
include("mSet.jl")
using MAP
using mSet
    global const ITER_MAX = 1e5
    global const ZERO_VAL = 1e-7

    srand(4)
A = [10 1; 1 20.]
a = rand(2)
B = [30 0; 0 5]
b = rand(2)
lbd = 1
n = 2
PA(xzero; lbd = 1, n=2)  = (lbd*A+eye(n))\(xzero + lbd*a)
PB(xzero; lbd = 1, n=2)  = (lbd*B+eye(n))\(xzero + lbd*b)
xzero = rand(n)

function Reflection(PA,xzero::Vector{Float64})
        proj_xzero = PA(xzero)
        return  2*proj_xzero - xzero

end

function DRMNonconvex(PA,PB,n::Int64,xzero::Vector{Float64};
                           printfile::Bool=false,EPS_VAL::Float64=1e-3,file::IOStream=IOStream("tables/DRM.table"))
    tol = 1e-3
    k = 1
    rate = 0.
    xstar = xzero
    while (tol >= EPS_VAL)
        # xbar = xstar #Looking for FixingPoint
        xold = xstar
        ypto = Reflection(PA,xstar)
        zpto = Reflection(PB,ypto)
        # ypto = Reflection(A,a,FA,xstar) # New method of Reflection
        # zpto = Reflection(B,b,FB,ypto) # New method of Reflection
        xstar = 0.5*(xstar + zpto)
        # println(xstar)
        # gapA = Projection(PA,aP,xstar)
        # gapB = Projection(PB,bP,xstar)
        # gapA = Projection(A,a,FA,xstar) # New method of Reflection
        # gapB = Projection(B,b,FB,xstar) # New method of Reflection
        # tol = norm(gapA - gapB,2) #Gap Distance
        # tol = norm([A;B]*xstar - [a;b],2) #True Error
        tol = norm(xold-xstar,2)
        # rate = tol/norm(xold-xbar,2)
        k +=2
        (k>ITER_MAX) && (warn("Maximum number of projections achieved in DRM"); conv=String("itmax"); break)
    end
    proj= PA(xstar)
    @show norm(A*proj - a + B*proj - b)
    
    # println(PB(xstar))
    # @show norm(PB(xstar) - xstar)


    if printfile
            str = @sprintf("    %s  %7d  %10.8e %10.8f ",conv,k,tol, rate)
            print(file,str)
            # println("Number of projections: $k")
    else
        println("DR Method")
        println("Number of projections: $k")
        @printf("Error = %s\n", tol)
        return k, tol, xstar

    end
end

function CRMNonconvex(PA,PB,n::Int64,xzero::Vector{Float64};
                           printfile::Bool=false,EPS_VAL::Float64=1e-3,file::IOStream=IOStream("tables/DRM.table"))
    tol = 1e-3
    k = 1
    rate = 0.
    xstar = xzero
    while (tol >= EPS_VAL)
        # xbar = xstar #Looking for FixingPoint
        # println(xstar)
        xold = xstar
        ypto = Reflection(PA,xstar)
        zpto = Reflection(PB,ypto)
        # norm(ypto - xstar)
        # norm(zpto - ypto)
        # norm(zpto - xstar)
        if norm(ypto - xstar)<ZERO_VAL
            xold = xstar
            xstar = 0.5*(xstar + zpto)
            # tol = norm([A;B]*xstar - [a;b],2)
            # gapA = Projection(PA,aP,xstar)
            # gapB = Projection(PB,bP,xstar)
            # tol = norm(gapA - gapB,2) #Gap Distance
            # tol = norm(xstar - xbar,2) #True Error
            # rate = tol/norm(xold - xbar,2)
            tol = norm(xold-xstar,2)
            continue
        elseif norm(zpto - ypto)<ZERO_VAL
            xold = xstar
            xstar = 0.5*(xstar + zpto)
            # gapA = Projection(PA,aP,xstar)
            # gapB = Projection(PB,bP,xstar)
            # tol = norm(gapA - gapB,2) #Gap Distance
            # tol = norm(xstar - xbar,2) #True Error
            tol = norm(xold-xstar,2)
            # rate = tol/norm(xold - xbar,2)
            continue
        elseif norm(zpto - xstar)<ZERO_VAL
            xold = xstar
            xstar = 0.5*(xstar + ypto)
            # gapA = Projection(PA,aP,xstar)
            # gapB = Projection(PB,bP,xstar)
            # tol = norm(gapA - gapB,2) #Gap Distance
            # tol = norm(xstar - xbar,2) #True Error
            # rate = tol/norm(xold - xbar,2)
            tol = norm(xold-xstar,2)
            continue    
        end


        xstar = MAP.FindCircumcenter2(xold,ypto,zpto)
        proj= PA(xstar)
        @show norm(A*proj - a + B*proj - b)
        # println(xstar)
        # gapA = Projection(PA,aP,xstar)
        # gapB = Projection(PB,bP,xstar)
        # gapA = Projection(A,a,FA,xstar) # New method of Reflection
        # gapB = Projection(B,b,FB,xstar) # New method of Reflection
        # tol = norm(gapA - gapB,2) #Gap Distance
        # tol = norm([A;B]*xstar - [a;b],2) #True Error
        tol = norm(xold-xstar,2)
        # rate = tol/norm(xold-xbar,2)
        k +=2
        (k>ITER_MAX) && (warn("Maximum number of projections achieved in CRM"); conv=String("itmax"); break)
    end
    if printfile
            str = @sprintf("    %s  %7d  %10.8e %10.8f ",conv,k,tol, rate)
            print(file,str)
            # println("Number of projections: $k")
    else
        println("Circumcenter Method")
        println("Number of projections: $k")
        @printf("Error = %s\n", tol)
        return k, tol, xstar
    end
end


k, tol, xstar = DRMNonconvex(PA,PB,n,xzero,EPS_VAL=1e-8)
# println(xstar)
# CRMNonconvex(PA,PB,n,xstar,EPS_VAL=1e-8)
CRMNonconvex(PA,PB,n,xstar+randn(2),EPS_VAL=1e-8)


