__precompile__()

module mSet
    import MAP
    export algorithmCDRMmSet, GenerateSampleApart

    global const ZERO_VAL = 1e-15
    global const ITER_MAX = 1e8


    # println("Tolerance: $EPS_VAL")
####################################
####################################
    function GenerateSampleApart(n::Int64,m::Int64)
        APartition = rand(1:n,m,1)
        R = sum(APartition)
        AA = []
        for I = 1:m
            push!(AA,randn(APartition[I],n))
        end
        return AA, R
    end
####################################




####################################
    function algorithmCDRMmSet(n::Int64,m::Int64,PA,aP,xstar::Vector{Float64},xbar::Vector{Float64},#\file::IOStream,
                            printfile=false,EPS_VAL::Float64=1e-3)
        
        tol = 1.0
        conv =  String("conv")
        # xstar =  Projection(PB,bP,xstar)  # Initial Point on span(U,V)
        k = 1
        rate = 0.
        # X is the matrix of the reflections
        X = zeros(n,m+1)
        xstar = MAP.Projection(PA[2],aP[2],xstar)
        # println(xstar)
        while (tol >= EPS_VAL)
            X[:,1] = xstar
            # xbar = xstar #Looking for FixingPoint
            # println("Entrou no loop $k")
            for J = 1:m
                X[:,J+1] = MAP.Reflection(PA[J],aP[J],X[:,J])
                # println(X[:,J+1])
            end
            # println("Criou X")
            xstar = FindCircumcentermSet(X)
            # println(xstar)
            # gapA = Projection(PA,aP,xstar)
            # gapB = Projection(PB,bP,xstar)
            # tol = norm(gapA - gapB,2) #Gap Distance
            tol = norm(xstar - xbar,2) #True Error
            # println(tol)
            # rate = tol/norm(xold - xbar,2)
            k += m
            (k>ITER_MAX) && (warn("Maximum number of projections achieved in C-DRM"); conv=String("itmax"); break)
        end
        # println(rate)
        if printfile
                str = @sprintf("    %s  %7d  %10.8e %10.8f",conv,k,tol,rate)
                print(file,str)
                # println("Number of projections: $k")
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
    end
####################################

####################################
      function algorithmMAPmSet(n::Int64,m::Int64,PA,aP,xstar::Vector{Float64},xbar::Vector{Float64},#\file::IOStream,
                                printfile=false,EPS_VAL::Float64=1e-3)
                               # Projecting onto A first
        tol = 1
        k = 0
        conv=String("itmax")
        while (tol >= EPS_VAL)
            # xbar = xstar #Looking for FixingPoint
            for J = 1:m
                xstar = MAP.Projection(PA[J],aP[J],xstar)
            end
            tol = norm(xstar - xbar,2) #True Error

            k+=m
            (k>ITER_MAX) && (warn("Maximum number of projections achieved in MAP"); conv=String("itmax"); break)
        end
        if printfile
                str = @sprintf("    conv %7d  %10.8e",k,tol)
                print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
    end
####################################


    ####################################
    function FindCircumcentermSet(X)
    # Finds the Circuncenter of poitns X = [x0, x1, x2, ...]
        Xrows = size(X,1)
        Xcols = size(X,2)
        if Xcols == 1
            return X
        end
        M = zeros(Xcols-1,Xrows)
        b = zeros(Xcols-1,1)
        x0 = X[:,1]
        for J = 1:Xcols-1
            xatual = X[:,J+1]
            v = xatual - x0
            M[J,:] = v'
            b[J] =0.5*dot(v,v)
        end 
        rankM = rank(M)
        if rankM < Xcols-1 
         warn("Rank deficient matrix")
        end
        MCond = cond(M)
        # println(M)
        # open("MatrixCond","a") do f
        #    # do stuff with the open file
        #   str = @sprintf("Matrix Cond %10.8f\n",MCond)
        #   print(f,str)
        # end
        # println(MCond)
        r = M\b
        return x0+r
    end
        ####################################
end
