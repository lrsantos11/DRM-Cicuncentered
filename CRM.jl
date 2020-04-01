__precompile__()

# module CRM
    if VERSION >=  v"0.7.0"
        using SparseArrays
        using LinearAlgebra
        using SuiteSparse
        using Printf
    end
    mutable struct Results
             residual::Float64
             iterations::Int64
             elapsedtime::Float64
             solution::AbstractArray{Float64}
    end
    # using mSet
    # export algorithmAARM, FriedrichsAngleAB, algorithmMAP, fourpointscheme, algorithmDRM
    # export algorithmCDRM, GenerateRdmPair, friedrichs, contructProjector, Projection, StartingPoint

    # global const ZERO_VAL = 1e-15
    # global const ITER_MAX = 1e5


    # println("Tolerance: $EPS_VAL")
####################################
####################################
    function GenerateSamples(n::Int64,affine::Bool=true)
        ma = rand(1:n-1) ## number of extra normals of A
        mb = rand(1:n-1-ma) ## number of extra normals of B
        A = randn(ma,n)
        B = randn(mb,n)
        if affine
            a = randn(ma)
            b = randn(mb)
        else
            a = zeros(ma)
            b = zeros(mb)
        end

        return A, a, ma, B, b, mb
    end
####################################
####################################
    function GenerateRdmPair(n::Int64,cmax::Int64,affine::Bool=true)
        n >= cmax ||  error("n must be greater then common normals")
        mcex = rand(1:cmax) ## number of common normals
        maex = rand(1:n-3*mcex) ## number of extra normals of A
        mbex = rand(1:n-2*mcex-maex) ## number of extra normals of B

        Cex = randn(mcex,n)
        while !(rank(Cex) == mcex)
            Cex = randn(mcex,n)
        end

        Aex = randn(maex,n)
        while !(rank([Aex; Cex]) == maex +  mcex)
            Aex = randn(maex,n)
        end

        Bex = randn(mbex,n)
        while !(rank([Aex; Cex; Bex]) == maex +  mcex + mbex)
            Bex = randn(mbex,n)
        end

        ma = maex + mcex
        mb = mbex + mcex

        # Space U
        A = [Aex; Cex]

        # Space V
        B = [Cex; Bex]

        if affine
            a = randn(ma)
            b = randn(mb)
        else
            a = zeros(ma)
            b = zeros(mb)
        end
        return A, a, ma, B, b, mb
    end
####################################

####################################
function friedrichs(A,B)

    cdm = rank(A)+rank(B)-rank([A;B]);
    C=[A;B[cdm+1:rank(B),:]];

    AA=pinv(A)*A;
    BB=pinv(B)*B;
    CC=pinv(C)*C

    M=CC+BB*AA-BB-AA;

    return sqrt(eigmax(M*M'));

end
####################################

####################################
    function FriedrichsAngleAB(A,B)
        # Calculating the Friederich Angle
        QA, RA = qr(A')
        QB, RB = qr(B')
        # S =  svd(QA'*QB)
        # Angle in Degrees
        angleAB = eigmax(QA'*QB)
        return acos(angleAB)
        # Angle in Radians
        # angleAB = acos(maximum(S))
        # println(S[2])
        # ind = findfirst(x -> x<(1-1e-8),S[2])
        # return maximum(S[2])
        # return acos(S[2][ind])
    end
####################################
####################################
    function contructProjector2(A::AbstractArray{Float64,2},a::Vector{Float64},n::Int64)
        MA = pinv(A)
        PA=eye(n) - MA*A
        aP = MA*a
        return PA, aP
    end
####################################
    function contructProjector(A::AbstractArray{Float64,2},a::Vector{Float64},n::Int64)
        QA, RA = qr(A')
        PA= eye(n) - QA*QA'
        aP = QA*(RA'\a)
        return PA, aP
    end
####################################
    function Projection(PA::AbstractArray{Float64,2},aP::Vector{Float64},xzero::Vector{Float64})
        # Orthogonal projection into affine subspace Ax=b if m<n.
        return PA*xzero + aP

    end
####################################
    function Reflection(PA::AbstractArray{Float64,2},aP::Vector{Float64},xzero::Vector{Float64})
        # Reflection from affine subespace Ax=b if m<n.
        # using LinearLeastSquares: min||x-xzero|| s.t. Ax=b.
        # Rbar(xzero) = R(xzero-b) + b
        # R = 2P - I where P is a projection
        proj_xzero = Projection(PA,aP,xzero);
        return  2*proj_xzero - xzero

    end
####################################
    function Projection(A,b,SF::SuiteSparse.SPQR.QRSparse{Float64,Int64},xstar)
    # Sparse Projection. SF should be computed as SF = qr(sparse(permutedims(A)))
    # A^T(AA^T)^{-1} = QR^{-T}
        if iszero(b)
            sproj = A*xstar
        else
            sproj = A*xstar - b
        end
        sproj = (permutedims(SF.R))\(sproj[SF.cpiv])
        sproj = SF.Q*(sproj)
        sproj = sproj[SF.rpivinv]
        return (xstar - sproj)
    end
####################################
    function Projection(B,b,SF::LinearAlgebra.QRCompactWY{Float64,Array{Float64,2}},xstar)
    # Dense Projection. SFB should be computed as F = qr(permutedims(A))
    # A^T(AA^T)^{-1} = QR^{-T}
        if iszero(b)
            proj = B*xstar
        else
            proj = B*xstar - b
        end
        proj = (permutedims(SF.R))\proj
        proj = SF.Q*(proj)
        # proj = pinv(B)*proj
        return (xstar - proj)
    end
####################################
    function Reflection(A,b,SF,xstar)
        # Sparse or Dense Reflection. SF should be computed as SF = qr(A')
        # Uses Projection Function either Sparse or Dense.
        # Reflection onto affine subespace Ax=b where m<n.
        proj_xstar = Projection(A,b,SF,xstar)
        return  2. *proj_xstar - xstar

    end
####################################
    function  fourpointscheme(A::AbstractArray{Float64,2},a::Vector{Float64},B::AbstractArray{Float64,2},
                                        b::Vector{Float64},n::Int64,xzero::Vector{Float64},
                                        file::IOStream,printfile=true,method::Int64=1, EPS_VAL::Float64=1e-3)
    #algoritmoRBLRS using 4 points in inner interaction
       #begining the projection onto A and onto B
        PA, aP = contructProjector(A,a,n)
        PB, bP = contructProjector(B,b,n)
        xstarA = Projection(PA,aP,xzero)
        xstarB = Projection(PB,bP,xzero)
        tol_xstarA = norm([A; B]*xstarA - [a; b],2)
        tol_xstarB = norm([A; B]*xstarB - [a; b],2)
        if tol_xstarA < tol_xstarB
            xstar = xstarA
            ind = 1
        else
            xstar = xstarB
            ind = 2
        end
    #   use ind = randi([1 2]) to implement random start
        tol = 1
        k =2
        conv=String("itmax")
        while (tol >= EPS_VAL)
            if (mod(ind,2) == 0)
                    b1 = xstar
                    #projection onto A
                    a1 = Projection(PA,aP,b1)
                    #projection onto B
                    b2 = Projection(PB,bP,a1)
                    #projection onto A
                    a2 = Projection(PA,aP,b2)

            else
                    a1 = xstar
                    #projection onto B
                    b1 = Projection(PB,bP,a1)
                    #projection onto A
                    a2 = Projection(PA,aP,b1)
                    #projection onto B
                    b2 = Projection(PB,bP,a2)
            end
            if method == 1
                # No Acceleration
                atil = a2
                btil = b2
            elseif method == 2
                # Closest Point Acceleration
                atil = accCP(a1,b1,a2)
                btil = accCP(b1,a2,b2)
            elseif method == 3
                atil = accEMAP(a1,b1,a2)
                btil = accEMAP(b1,a2,b2)
            end

            # Lembrete: Em Omega Geral, projetar atil/btil em OmA ou OmB
            tol_A = norm([A; B]*a2 - [a; b],2)
            tol_B = norm([A; B]*b2 - [a; b],2)
            tolac_A = norm([A; B]*atil - [a; b],2)
            tolac_B = norm([A; B]*btil - [a; b],2)
            vectors = [a2 b2 atil btil]
            tol, ind = findmin([tol_A tol_B tolac_A tolac_B])
            xstar = vectors[:,ind]
            k+=3
            (k>ITER_MAX) && (warn("Maximum number of projections achievied in 4point"); conv=String("itmax"); break)
        end
        if printfile
            str = @sprintf("    conv %7d  %10.8e",k,tol)
            print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
        # println(xstar)
    end
####################################
####################################
    function  algorithmMAP(A::AbstractArray{Float64,2},a::Vector{Float64},B::AbstractArray{Float64,2},
                           b::Vector{Float64},n::Int64,xzero::Vector{Float64},xbar::Vector{Float64},
                           printfile::Bool=false,method::Int64=1,EPS_VAL::Float64=1e-3,file::IOStream=IOStream("tables/MAP.table"))
        # Projecting onto A first
        # PA, aP = contructProjector(A,a,n)
        # PB, bP = contructProjector(B,b,n)
        # xstar = Projection(PA,aP,xzero)
        FA = qr(Matrix(A'))
        FB = qr(Matrix(B'))
        xstar = Projection(A,a,FA,xzero)
        toxzel = 1
        k = 1
        conv=String("itmax")
        # println("entrando no loop")
        while (tol >= EPS_VAL)
            # xbar = xstar #Looking for FixingPoint
            a1 = xstar
            #projection onto B
            # b1 = Projection(PB,bP,a1)
            b1 = Projection(B,b,FB,a1) #New Method of Projection
            #projection onto A
            # a2 = Projection(PA,aP,b1)
            a2 = Projection(A,a,FA,b1) # New Method of Projection
            if method == 1
                # RAS Acceleration
                xstar = a2
            elseif method == 2
                # Closest Point Acceleration
                xstar = accCP(a1,b1,a2)
            elseif method == 3
                # Extrapolated Acceleration
                xstar = accEMAP(a1,b1,a2)
            end
            # tol = norm(a2 - b1,2) #Gap Distance
            tol = norm(xstar - xbar,2) #True Error

            k+=2
            # println(k)
            (k>ITER_MAX) && (warn("Maximum number of projections achieved in MAP"); conv=String("itmax"); break)
        end
        if printfile
                str = @sprintf("    conv %7d  %10.8e",k,tol)
                print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
            return k, tol
        end
    end
####################################
####################################
    function algorithmDRM(A::AbstractArray{Float64,2},a::Vector{Float64},B::AbstractArray{Float64,2},
                           b::Vector{Float64},n::Int64,xzero::Vector{Float64},xbar::Vector{Float64},
                           printfile::Bool=false,EPS_VAL::Float64=1e-3,file::IOStream=IOStream("tables/DRM.table"))
        # PA, aP = contructProjector(A,a,n)
        # PB, bP = contructProjector(B,b,n)
        FA = qr(Matrix(A'))
        FB = qr(Matrix(B'))
        # tol = norm([A;B]*xstar - [a;b],2)
        tol = 1.
        conv =  String("conv")
        # Initial Point on span(U,V)
        # xstar =  Projection(PB,bP,xzero)
        xstar =  Projection(B,b,FB,xstar)  # New method of Projection
        k = 1
        rate = 0.
        while (tol >= EPS_VAL)
            # xbar = xstar #Looking for FixingPoint
            xold = xstar
            # ypto = Reflection(PA,aP,xstar)
            # zpto = Reflection(PB,bP,ypto)
            ypto = Reflection(A,a,FA,xstar) # New method of Reflection
            zpto = Reflection(B,b,FB,ypto) # New method of Reflection
            xstar = 0.5*(xstar + zpto)
            # gapA = Projection(PA,aP,xstar)
            # gapB = Projection(PB,bP,xstar)
            gapA = Projection(A,a,FA,xstar) # New method of Reflection
            gapB = Projection(B,b,FB,xstar) # New method of Reflection
            # tol = norm(gapA - gapB,2) #Gap Distance
            tol = norm([A;B]*xstar - [a;b],2) #True Error
            rate = tol/norm(xold-xbar,2)
            k +=2
            (k>ITER_MAX) && (warn("Maximum number of projections achieved in DRM"); conv=String("itmax"); break)
        end
        if printfile
                str = @sprintf("    %s  %7d  %10.8e %10.8f ",conv,k,tol, rate)
                print(file,str)
                # println("Number of projections: $k")
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
            return k, tol
        end
    end
####################################
####################################
    function algorithmCDRM(A,a::Vector{Float64},B,
                           b::Vector{Float64},n::Int64,xzero::Vector{Float64},xbar::Vector{Float64},
                           printfile::Bool=false,EPS_VAL::Float64=1e-3,file::IOStream=IOStream("tables/CDRM.table"))
        # PA, aP = contructProjector(A,a,n)
        # PB, bP = contructProjector(B,b,n)
        # println("Entrou no Algoritmo")
        FA = qr(Matrix(A'))
        FB = qr(Matrix(B'))
        tol = 1.0
        conv =  String("conv")
        # Initial Point on span(U,V)
        # xstar =  Projection(PB,bP,xzero)
        xstar =  Projection(B,b,FB,xzero)  # New method of Projection
        # Initial Point in (U+V) (Roger asked)
        # xzeroPA =  Projection(PA,aP,xzero)
        # xzeroPB =  Projection(PB,bP,xzeroPA)
        # xstar = 0.5*(xzeroPA + xzeroPB)
        k = 1
        rate = 0.
        while (tol >= EPS_VAL)
            # println("Entrou no Loop")
            # xbar = xstar #Searching for FixedPoint
            # ypto = Reflection(PA,aP,xstar)
            # zpto = Reflection(PB,bP,ypto)
            ypto = MAP.Reflection(A,a,FA,xstar) # New method of Reflection
            zpto = MAP.Reflection(B,b,FB,ypto) # New method of Reflection
            # println("Reflections done")
            if norm(ypto - xstar)<ZERO_VAL
                xold = xstar
                xstar = 0.5*(xstar + zpto)
                # tol = norm([A;B]*xstar - [a;b],2)
                # gapA = Projection(PA,aP,xstar)
                # gapB = Projection(PB,bP,xstar)
                # tol = norm(gapA - gapB,2) #Gap Distance
                tol = norm(xstar - xbar,2) #True Error
                rate = tol/norm(xold - xbar,2)
                continue
            elseif norm(zpto - ypto)<ZERO_VAL
                xold = xstar
                xstar = 0.5*(xstar + zpto)
                # gapA = Projection(PA,aP,xstar)
                # gapB = Projection(PB,bP,xstar)
                # tol = norm(gapA - gapB,2) #Gap Distance
                tol = norm(xstar - xbar,2) #True Error
                rate = tol/norm(xold - xbar,2)
                continue
            end
            xold = xstar
            xstar = mSet.FindCircumcentermSet([xold ypto zpto])
            # xstar = FindCircumcenter2(xold,ypto,zpto)

            # gapA = MAP.Projection(A,a,FA,xstar)
            # gapB = MAP.Projection(B,b,FB,xstar)
            # tol = norm(gapA - gapB,2) #Gap Distanc
            # tol = norm([A;B]*xstar - [a;b],2) #True Error
            tol = norm(xstar - xbar,Inf) #True Error
            # println("Estou aqui")
            # rate = tol/norm(xold - xbar,2)
            @show k += 2
            @show tol
            # println(k)
            # println(xstar)
            (k>ITER_MAX) && (warn("Maximum number of projections achieved in C-DRM"); println(tol); conv=String("itmax"); break)
        end
        # println(rate)
        if printfile
                str = @sprintf("    %s  %7d  %10.8e %10.8f",conv,k,tol,rate)
                print(file,str)
                # println("Number of projections: $k")
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
            return k, tol
        end
    end
####################################
####################################
    function algorithmAARM(A::AbstractArray{Float64,2},a::Vector{Float64},B::AbstractArray{Float64,2},
                           b::Vector{Float64},n::Int64,xzero::Vector{Float64},xbar::Vector{Float64},
                           printfile::Bool=false,EPS_VAL::Float64=1e-3,file::IOStream=IOStream("tables/AARM.table"))
        # PA, aP = contructProjector(A,a,n)
        # PB, bP = contructProjector(B,b,n)
        FA = qr(Matrix(A'))
        FB = qr(Matrix(B'))
        AngleAB = FriedrichsAngleAB(A,B)
        alpha = 1.0
        beta_opt=1.0/(1.0+sin(AngleAB))
        # Initial Point on span(U,V)
        # xstar =  Projection(PB,bP,xzero)
        xstar =  Projection(B,b,FB,xzero)  # New method of Projection
        # Initial Point in (U+V) (Roger)
        # xzeroPA =  Projection(PA,aP,xzero)
        # xzeroPB =  Projection(PB,bP,xzeroPA)
        # xstar = 0.5*(xzeroPA + xzeroPB)
        k = 1
        rate = 0.
        err=10
        q = xstar
        while err>EPS_VAL
            proj_A = Projection(A,a,FA,(xstar+q))
            QA=2*beta_opt*(proj_A-q)-xstar
            proj_B = Projection(B,b,FB,(QA+q))
            QB=2*beta_opt*(proj_B-q)-QA
            xnew=(1-alpha)*xstar+alpha*QB
            zn=Projection(A,a,FA,(xnew+q))
            err=norm(zn-xbar)
            k+=3nso
            xstar=xnew
        end
        if printfile
                str = @sprintf("    %s  %7d  %10.8e %10.8f",conv,k,err,rate)
                print(file,str)
                # println("Number of projections: $k")
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", err)
            return k, err
        end
    end
####################################

####################################
    function accEMAP(x1,x2,x3)

        # x1 and x3 are in subspace A
        # x3 is in subspace B.
        # Finds the Extrapolation VonNewmann-Accelerated
        # point xacc in subspace given by span(x3-x1)
        t = norm(x1-x2)^2/norm(x3-x1)^2
        return x1 + (x3 - x1)*t
    end
####################################
####################################
    function accCP(x1,x2,x3)
    # Finds the Closest Point Acceleration xac
    # x2 = PB*x1 in B if x1 in A.
    # x3 = PA*PB*x1, if x1 in A.

        # if (x1 ~= x3)
            t = dot(x1,x1-x3)/((norm(x1 - x3,2))^2)
        # else
            # t = 1
        # end
        return x1 + t*(x3 -x1)

    end
####################################
    function FindCircumcenter(x0,x1,x2)
    # Finds the Circumcenter of points x0, x1, x2
        v1 = x1 - x0
        v2 = x2 - x0
        M = [v1'; v2']
        # MCond = cond(M)
        # open("MatrixCond","a") do f
        #    # do stuff with the open file
        #   str = @sprintf("Matrix Cond %10.8f\n",MCond)
        #   print(f,str)
        # end
        # println("Method 1: ",MCond)
        r = M\([0.5*dot(v1,v1); 0.5*dot(v2,v2)])
        return x0+r
    end
####################################

    function FindCircumcenter2(x,y,z)
    # Finds the Circumcenter of points x, y, z
        su = y - x
        sv = z - x
        nsu2 = norm(su)^2.0
        nsv2 = norm(sv)^2.0
        sudotsv = dot(su,sv)
        M = [nsu2 sudotsv; sudotsv nsv2]
        @show rank(M)
        # MCond = cond(M)
        # open("MatrixCond","a") do f
        #    # do stuff with the open file
        #   str = @sprintf("Matrix Cond %10.8f\n",MCond)
        #   print(f,str)
        # end
        # println("Method 2: ",MCond)
        r = M\([0.5*nsu2; 0.5*nsv2])
        return x+r[1]*su + r[2]*sv
    end

####################################
    function StartingPoint(n::Int64)
        ## Creates a random point in R^n

        x=zeros(n);
        while norm(x)<2
            x = randn(n);
        end
        # norm between 5 and 15
        foonorm = (15-5)*rand() + 5
        return foonorm*x/norm(x);
    end
####################################

    """
    TestConvex(xzero,itmax)
    """

function TestConvex(xzero,itmax,ReflectA,ReflectB,ProjectA,ProjectB, EPS_VAL = 1e-6, printfile::Bool = true, xbar = [])
        if printfile
            try
                rm("tables/xDR.dat")
            catch
            end
            try
                rm("tables/xCRM.dat")
            catch
            end
            try
                rm("tables/xMAP.dat")
            catch
            end
        end
        n_col = length(xzero)
        xDR = xzero
        # Douglas-Rachford
        k = 1
        tolDR = 1.
        while tolDR > EPS_VAL   && k <= itmax
            # DR
            xDRold = xDR
            # println(xDR)
            yDR = ReflectA(xDR)
            # println(yDR)
            zDR = ReflectB(yDR)
            # println(zDR)
            xDR = 0.5*(xDR + zDR)
            # tolDR = norm(xDR - xbarDR,2)
            tolDR = norm(xDR - xDRold,2)
            # tolDR = norm(xDR  - ProjectA(xDR))
            # tolDR = norm(min.(0,ProjectB(xDR)))
            open("tables/xDR.dat","a") do file
            #    # do stuff with the open file
                    writedlm(file,[k tolDR])
            end
            k += 1
        end
        elapsedtime = 0.0
        # println("xDR: $xDR")
        println("Projections of  DR: $k")
        println("Norm for xDR: $tolDR")
        # println(xDR)
        resultDR = Results(tolDR,k,elapsedtime,xDR)
        # CRM
        k = 1
        tolCRM = 1.
        xCRM = xzero


        while tolCRM > EPS_VAL && k <= itmax
            # CRM
            xCRMold = xCRM
            ypto = ReflectA(xCRM)
            zpto = ReflectB(ypto)
            if norm(ypto - xCRM)<ZERO_VAL
                xCRM = 0.5*(xCRM + zpto)
                # tolCRM = norm(xCRM - xbar,2)
                tolCRM = norm(xCRM - xCRMold,2)
                # tolCRM = norm(xCRM - ProjectA(xCRM))
                # tolCRM = norm(min.(0,ProjectB(xCRM)))
                open("tables/xCRM.dat","a") do file
                # do stuff with the open file
                        writedlm(file,[k tolCRM])
                end
                k += 1
                continue
            elseif norm(zpto - ypto)<ZERO_VAL
                xCRM = 0.5*(xCRM + zpto)
                # tolCRM = norm(xCRM - xbar,2)
                tolCRM = norm(xCRM - xCRMold,2)
                # tolCRM = norm(xCRM - ProjectA(xCRM))
                open("tables/xCRM.dat","a") do file
                #    # do stuff with the open file
                        writedlm(file,[k tolCRM])
                end
                k += 1
                continue
            end
            xCRM = FindCircumcentermSet([xCRM, ypto, zpto])
            # println("xpto: $xCRM")
            # tolCRM = norm(xCRM - xbar,2)
            tolCRM = norm(xCRM - xCRMold,2)
            # tolCRM = norm(xCRM - ProjectA(xCRM))
            # tolCRM = norm(min.(0,ProjectB(xCRM)))
            open("tables/xCRM.dat","a") do file
            #    # do stuff with the open file
                    writedlm(file,[k tolCRM])
            end
            k += 1
        end
        # println("xC: $xCRM")
        println("Projections of CRM: $k")
        println("Norm for xC: $tolCRM")
        resultCRM = Results(tolCRM,k,elapsedtime,xCRM)

        # MAP
        k = 1
        tolMAP = 1.
        xMAP = xzero
        # open("tables/xMAP.dat","a") do file
        #    # do stuff with the open file
        #   writedlm(file,xMAP')
        # end

        while tolMAP > EPS_VAL && k <= itmax
            # MAP
            xMAPold = xMAP
            # open("tables/xMAP.dat","a") do file
                ypto = ProjectA(xMAP)
            #     writedlm(file,ypto')
                xMAP = ProjectB(ypto)
            #     writedlm(file,xMAP')
            # end
                # tolMAP = norm(xMAP - xbar,2)
                tolMAP = norm(xMAP - xMAPold,2)
                # tolMAP = norm(xMAP - ProjectA(xMAP))
                # tolMAP = norm(min.(0,ProjectB(xMAP)))
                open("tables/xMAP.dat","a") do file
                    # do stuff with the open file
                            writedlm(file,[k tolMAP])
                end
                # println(tolMAP)
                k += 1
               # do stuff with the open file

        end
        # println("xC: $xMAP")
        println("Projections of MAP: $k")
        println("Norm for MAP: $tolMAP")
        resultMAP = Results(tolMAP,k,elapsedtime,xMAP)

        return resultDR, resultCRM, resultMAP
    end
####################################



# end
