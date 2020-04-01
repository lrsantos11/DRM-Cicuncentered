# workspace()
# using JLD
# using HDF5
# using BenchmarkProfiles
using BenchmarkTools
using MatrixMarket
using SparseArrays
using LinearAlgebra
using DelimitedFiles
using ProximalOperators
using Plots
using LaTeXStrings
# pgfplots()
# PGFPlots.pushPGFPlotsPreamble("\\usepackage{amssymb}")
# unicodeplots()
pyplot()
PyPlot.rc("text", usetex= true)
PyPlot.rc("text.latex",preamble = ["\\usepackage{amsmath,amssymb}", "\\usepackage[utf8]{inputenc}"])

include("MAP.jl")
include("mSet.jl")
global const EPS_VAL = 1e-5
global const ZERO_VAL = 1e-13
global const ITER_MAX = 1e5

    function GeneralTest(n::Int64 = 100,samples::Int64 = 2,restarts::Int64=1,
        affine::Bool=false,EPS_VAL::Float64 = 1e-3,printfile::Bool=false)
        # Fix Random
        Random.seed!(2)
        # X = R^n
        j = 1
        firstwrite = true
        resultsMAP = Array{Float64}(0,5)
        resultsDRM = Array{Float64}(0,5)
        resultsCDRM = Array{Float64}(0,5)
        resultsAARM = Array{Float64}(0,5)
        results4pts = Array{Float64}(0,5)
        while (j <=samples)
            #########################
            # Generate Random Pair of Subspaces
            cmax = rand(1:ceil(Integer,n/10,))
            CC  = GenerateRdmPair(n,cmax,affine)
            #########################
            #  Include Random Generated Pairs saved in JLD
            # Dict =  load("tables/highcF.jld")

            #########################
            #  Read Files
            A = CC[1]
            a = CC[2]
            ma = CC[3]
            B = CC[4]
            b= CC[5]
            mb = CC[6]
            #########################
            # Include Bauschke et all Pairs
            # fld = @sprintf("../Experiment_Julia/pairs2/pair%03d.mat",j)
            # include(fld) ## read A,B,AngleAB
            # println(A)
            # ma, n = size(A)
            # a = zeros(ma)
            # mb, n = size(B)
            # b = zeros(mb)

            # println(rank([A;B]))
            # println(size(A))
            # println(rank(A))
            # println(size(B))
            # println(rank(B))
            AngleAB = FriedrichsAngleAB(A,B)
            cF = cos(AngleAB)
            # println(AngleAB)
            # if AngleAB > 1e-2
                # println("Find better AngleAB")
                # continue
             # end
            # P, p = contructProjector([A; B],[a; b],n)
            QR = qrfact([A;B]')
            for i = 1:restarts
                xzero = StartingPoint(n)
                xbar = Projection([A;B],[a;b],QR,xzero)
                prob_name  = String("Problem$j"*"?-$i"*"A$ma"*"B$mb")
                println(prob_name)
                #########################
                # Saving data to file
                # jldopen("tables/highcF.jld", "r+") do file
                #     write(file, prob_name,CC)
                # end
                #########################
                # println("MAP with no acceleration")
                # if printfile
                #     fname = @sprintf("tables/N%d_Samp%d_%.1E-MAP-true.table",n,samples,EPS_VAL);
                #     open(fname,"a") do file
                #         if firstwrite
                #             str = @sprintf("---\nalgname: MAP\nsuccess: conv\nfree_format: True\n---\n")
                #             print(file,str)
                #         end
                #         write(file,prob_name)
                #         time = @elapsed algorithmMAP(A,a,B,b,n,xzero,xbar,printfile,1,EPS_VAL,file)
                #         str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                #     # print(str)
                #         print(file,str)
                #     end
                # else

                #     time = @elapsed k, err = algorithmMAP(A,a,B,b,n,xzero,xbar,printfile,1,EPS_VAL)
                #     resultsMAP = vcat(resultsMAP,[k err time AngleAB cF])
                # end
                #########################
                # println("MAP with CP acceleration")
                # @time MAP.algorithmMAP(A,a,B,b,n,xzero,2)
                #######################
                # println("MAP with EMAP acceleration")
                # fname = @sprintf("tables/N%d_Samp%d_%.1E-EMAP-true.table",n,samples,EPS_VAL);
                # open(fname,"a") do file
                #     if firstwrite
                #         str = @sprintf("---\nalgname: E-MAP\nsuccess: conv\nfree_format: True\n---\n")
                #         print(file,str)
                #     end
                #     write(file,prob_name)
                #     time = @elapsed algorithmMAP(A,a,B,b,n,xzero,xbar,file,true,3,EPS_VAL)
                #     str = @sprintf(" %10.8e %10.8f\n",time,AngleAB)
                #     # print(str)
                #     print(file,str)
                # end
                ########################
                # println("Four Point Scheme with no acceleration")
                # fname = @sprintf("tables/%d_%.1E-4pointMAP-true.table",samples,EPS_VAL);
                # open(fname,"a") do file
                #     write(file,prob_name)
                #     time = @elapsed fourpointscheme(A,a,B,b,n,xzero,file,true,1,EPS_VAL)
                #     str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                #     # print(str)
                #     print(file,str)
                # end
                #########################
                # println("Four Point Scheme with CP acceleration")
                # if printfile
                #     fname = @sprintf("tables/N%d_Samp%d_%.1E-C-DRM-true.table",n,samples,EPS_VAL);
                #     open(fname,"a") do file
                #         if firstwrite
                #             str = @sprintf("---\nalgname: FourPointScheme\nsuccess: conv\nfree_format: True\n---\n")
                #             print(file,str)
                #         end
                #         write(file,prob_name)
                #         time = @elapsed fourpointscheme(A,a,B,b,n,xzero,xbar,2,printfile,EPS_VAL,file)
                #         str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)
                #         # print(str)
                #         print(file,str)
                #     end
                # else
                #     time = @elapsed k, err = fourpointscheme(A,a,B,b,n,xzero,xbar,2,printfile,EPS_VAL)
                #     results4pts = vcat(results4pts,[k err time AngleAB cF])
                # end
                #########################
                # println("Four Point Scheme with EMAP acceleration")
                # fname = @sprintf("tables/%d_%.1E-4pointEMAP-true.table",samples,EPS_VAL);
                # open(fname,"a") do file
                #     write(file,prob_name)
                #     time = @elapsed fourpointscheme(A,a,B,b,n,xzero,xbar,file,true,3,EPS_VAL)
                #     str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                #     # print(str)
                #     print(file,str)
                # end

                # #######################
                println("C-DR Method")
                if printfile
                    fname = @sprintf("tables/N%d_Samp%d_%.1E-C-DRM-true.table",n,samples,EPS_VAL);
                    open(fname,"a") do file
                        if firstwrite
                            str = @sprintf("---\nalgname: C-DRM\nsuccess: conv\nfree_format: True\n---\n")
                            print(file,str)
                        end
                        write(file,prob_name)
                        time = @elapsed algorithmCDRM(A,a,B,b,n,xzero,xbar,printfile,EPS_VAL,file)
                        str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)
                        # print(str)
                        print(file,str)
                    end
                else
                    time = @elapsed k, err = algorithmCDRM(A,a,B,b,n,xzero,xbar,printfile,EPS_VAL)
                    resultsCDRM = vcat(resultsCDRM,[k err time AngleAB cF])
                end
                ####################
                # println("DR Method")
                # if printfile
                #     fname = @sprintf("tables/N%d_Samp%d_%.1E-DRM-true.table",n,samples,EPS_VAL);
                #     open(fname,"a") do file
                #         if firstwrite
                #             str = @sprintf("---\nalgname: DRM\nsuccess: conv\nfree_format: True\n---\n")
                #             print(file,str)
                #         end
                #         write(file,prob_name)
                #         time = @elapsed algorithmDRM(A,a,B,b,n,xzero,xbar,printfile,EPS_VAL,file)
                #         str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                    # print(str)
                #         print(file,str)
                #     end
                # else
                #     time = @elapsed  k, err = algorithmDRM(A,a,B,b,n,xzero,xbar,printfile,EPS_VAL)
                #     resultsDRM = vcat(resultsDRM,[k err time AngleAB cF])
                # end
                # ###################
                # println("AARM Method")
                # if printfile
                #     fname = @sprintf("tables/N%d_Samp%d_%.1E-AARM-true.table",n,samples,EPS_VAL);
                #     open(fname,"a") do file
                #         if firstwrite
                #             str = @sprintf("---\nalgname: AARM\nsuccess: conv\nfree_format: True\n---\n")
                #             print(file,str)
                #         end
                #         write(file,prob_name)
                #         time = @elapsed algorithmAARM(A,a,B,b,n,xzero,xbar,printfile,EPS_VAL,file)
                #         str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                    # print(str)
                #         print(file,str)
                #     end
                # else
                #     time = @elapsed  k, err = algorithmAARM(A,a,B,b,n,xzero,xbar,printfile,EPS_VAL)
                #     resultsAARM = vcat(resultsAARM,[k err time AngleAB cF])
                # end
                # #  ###################
                # firstwrite = false
            end
            j += 1
        end
        # jldopen("tables/resultshighcF.jld", "w") do file
        #     @write file resultsCDRM
        #     @write file resultsAARM
        # end
        # println(resultsCDRM)
        # scatter(resultsMAP[:,4],resultsMAP[:,1],label="MAP")
        # scatter!(resultsDRM[:,4],resultsDRM[:,1],label="DRM")
       s1 = scatter(resultsCDRM[:,4],resultsCDRM[:,1],label="CDRM")
       # s1 = scatter!(resultsAARM[:,4],resultsAARM[:,1],label="AARM",xlabel="Friederich Angle (in radians)",ylabel = "Number of projections", yscale = :log10, title= L"100 instances in $\mathbb{R}^{200}$ (True error -- $\varepsilon_2 = 10^{-6}$)")
       plot(s1)
       # performance_profile(hcat(resultsMAP[:,1],resultsDRM[:,1],resultsCDRM[:,1],resultsAARM[:,1]), ["MAP", "DRM", "CDRM", "AARM"], title="Comparing Number of Projections")
       # performance_profile(hcat(resultsMAP[:,3],resultsDRM[:,3],resultsCDRM[:,3],resultsAARM[:,3]), ["MAP", "DRM", "CDRM", "AARM"], title="Comparing Number of Projections")
       # p1 = performance_profile(hcat(resultsCDRM[:,1],resultsAARM[:,1]), ["CDRM", "AARM"], title="Performance Profile", ylabel = "Number of projections",title="Performance Profile")
       # plt = plot(s1,p1,layout = (1,2),legend = true)
       # save("tables/N$n""_Samp$samples""_$EPS_VAL""-true.pdf", plt)

    end


# GeneralTest(10,3,2,false)

# GeneralTest(50,100,1,false)

# GeneralTest(200,25,1,false)
# GeneralTestmSet(3,3,false)

    # function ReflecBall(x, v, r)
        # A is the Ball centered in v and ray r
            # if norm(x)<=1
            #     return x
            # end
            # v=[-sqrt(2)/2;sqrt(2)/2]
            # v=[-.5,0]
            # r = 1.0
            # proj = v + (x-v)/(max(1,norm(x-v)/r))
            # proj = v + (x-v)/(norm(x-v)/r)
            # return 2. *proj - x
    # end

    function ReflecBall(x, v, r)
        # A is the Ball centered in v
        #     if norm(x)<=1
        #         return x
        #     end
            # v=[.5,0.]
            # r = 1.
            proj = v + (x-v)/(max(1,norm(x-v)/r))
            return 2. * proj - x
    end
    function ProjectBall(x, v, r)
        # A is the Ball centered in v
        #     if norm(x)<=1
        #         return x
        #     end
            # v=[.5,0.]
            # r = 1.
            proj = v + (x-v)/(max(1,norm(x-v)/r))
            return proj
    end
    function ReflecAffine(x)
    #     # B is the  (affine) subspace that line that crosses [1;1]
            n = length(x)
            OneVec = ones(n)
            proj = dot(OneVec,x)/dot(OneVec,OneVec)*OneVec
            return 2. *proj - x
    end

    # ReflecA(x) =  ReflecBall(x,[-.5,0], 1.)

    # ReflecB(x) =  ReflecAffine(x)


    """
    ProjectionYequalZ(X)

    Returns proj = y = z of R^n, the orthogonal projection of X = [x1;x2] in R^{2n} onto the two dimensional subspace of tuples
    such that y=z.
    """
    function ProjectionYequalZ(X::AbstractArray)
        n = length(X)
        if !iseven(n)
            @error "Dimension of x must be even!"
        end
        n = div(n,2)
        x = X[1:n]
        y = X[n+1:end]
        proj = .5*(x+y)
        return [proj;proj]
    end

    """
    Reflection on the Product Space

    """
    function ReflecA(x)
        n = length(x)
        n = div(n,2)
        y = ReflecBall(x[1:n],[-1.,0], 1.)
        z = ReflecBall(x[1:n],[1.,0], 1.)
        return [y; z]
    end

    """
    Reflection on the Diagonal. Returns vector on

    """
    function ReflecB(x)
        proj = ProjectionYequalZ(x)
        reflec = 2. *proj - x
        return reflec
    end

    function testconvex(itmax; xzero = [])
        try
            rm("tables/xDR.dat")
        catch
        end
        try
            rm("tables/xCRM.dat")
        catch
        end
        if isempty(xzero)
            xDRC = [-3.,   -3.] #5.*randn(2) #[1.233205883326263 ; .9317100309585165]
        else
            xDRC = xzero
        end
        n = length(xDRC)
        OneVec = ones(2)
        xDR = xDRC
        # xbar = [0.; 0.]
        # Douglas-Rachford
        k = 0
        tolDR = 1.
        tolC = 1.
        open("tables/xDR.dat","a") do f
            # do stuff with the open file
            writedlm(f, xDR')
        end
        while tolDR > EPS_VAL   && k <= itmax
            # DR
            xDRold = xDR
            # println(xDR)
            yDR = ReflecA(xDR)
            # println(yDR)
            zDR = ReflecB(yDR)
            # println(zDR)
            xDR = 0.5*(xDR + zDR)
            # xDRT = dot(OneVec,xDR)/dot(OneVec,OneVec)*OneVec
            # tolDR = norm(xDRT - xbar,2)
            tolDR = norm(xDR - xDRold,2)
            k += 1
            # println(xDR)
            open("tables/xDR.dat","a") do f
               # do stuff with the open file
              writedlm(f,xDR')
            end
        end
        println("xDR: $xDR")
        println("Projections of  DR: $k")
        println("Norm for xDR: $tolDR")

        # CRM
        k = 0
        open("tables/xCRM.dat","a") do f
           # do stuff with the open file
          writedlm(f,xDRC')
        end

        while tolC > EPS_VAL && k <= itmax
            # CRM
            xDRCold = xDRC
            ypto = ReflecA(xDRC)
            zpto = ReflecB(ypto)
            # if norm(ypto - xDRC)<ZERO_VAL
            #     xDRC = 0.5*(xDRC + zpto)
            #     # tolC = norm(xDRC - xbar,2)
            #     tolC = norm(xDRC - xDRCold,2)
            #     k += 1
            #     open("tables/xCRM.dat","a") do f
            #        # do stuff with the open file
            #       writedlm(f,xDRC')
            #     end
            #     continue
            # elseif norm(zpto - ypto)<ZERO_VAL
            #     xDRC = 0.5*(xDRC + zpto)
            #     # tolC = norm(xDRC - xbar,2)
            #     tolC = norm(xDRC - xDRCold,2)
            #     k += 1
            #     open("tables/xCRM.dat","a") do f
            #        # do stuff with the open file
            #       writedlm(f,xDRC')
            #     end
            #     continue
            # end
            xDRC = FindCircumcentermSet([xDRC ypto zpto])
            # println("xpto: $xDRC")
            # tolC = norm(xDRC - xbar,2)
            tolC = norm(xDRC - xDRCold,2)
            k += 1
            open("tables/xCRM.dat","a") do f
               # do stuff with the open file
              writedlm(f,xDRC')
            end
        end
        println("xC: $xDRC")
        println("Projections of CRM: $k")
        println("Norm for xC: $tolC")

        # MAP
        k = 0
        open("tables/xMAP.dat","a") do f
           # do stuff with the open file
          writedlm(f,xDRC')
        end

        while tolC > EPS_VAL && k <= itmax
            # MAP
            xDRCold = xDRC
            ypto = ReflecA(xDRC)
            zpto = ReflecB(ypto)
            # if norm(ypto - xDRC)<ZERO_VAL
            #     xDRC = 0.5*(xDRC + zpto)
            #     # tolC = norm(xDRC - xbar,2)
            #     tolC = norm(xDRC - xDRCold,2)
            #     k += 1
            #     open("tables/xMAP.dat","a") do f
            #        # do stuff with the open file
            #       writedlm(f,xDRC')
            #     end
            #     continue
            # elseif norm(zpto - ypto)<ZERO_VAL
            #     xDRC = 0.5*(xDRC + zpto)
            #     # tolC = norm(xDRC - xbar,2)
            #     tolC = norm(xDRC - xDRCold,2)
            #     k += 1
            #     open("tables/xMAP.dat","a") do f
            #        # do stuff with the open file
            #       writedlm(f,xDRC')
            #     end
            #     continue
            # end
            xDRC = FindCircumcentermSet([xDRC ypto zpto])
            # println("xpto: $xDRC")
            # tolC = norm(xDRC - xbar,2)
            tolC = norm(xDRC - xDRCold,2)
            k += 1
            open("tables/xMAP.dat","a") do f
               # do stuff with the open file
              writedlm(f,xDRC')
            end
        end
        # println("xC: $xDRC")
        println("Projections of MAP: $k")
        println("Norm for xC: $tolC")
    end



    function GeneralTestmSet(n::Int64 = 100,m::Int64 = 2,affine::Bool=false)
        # Fix Random
        # srand(4)
        xstar = MAP.StartingPoint(n)
        # println(xstar)
        AA, R = GenerateSampleApart(n,m)
        println(AA)
        println(R)
        PA = []
        aP = []
        A = []
        for J = 1:m
            a  = zeros(size(AA[J],1))
            PA_J, aP_J = MAP.contructProjector(AA[J],a,n)
            push!(PA,PA_J)
            push!(aP,aP_J)
            if J == 1
                A = AA[J]
            else
                A = vcat(A,AA[J])
            end
        end
        # println(A)
        P, p = MAP.contructProjector(A,zeros(size(A,1)),n)
        xbar = MAP.Projection(P,p,xstar)
        println("C-DRM")
        time = @elapsed mSet.algorithmCDRMmSet(n,m,PA,aP,xstar,xbar,false,EPS_VAL)
        println("MAP")
        time = @elapsed mSet.algorithmMAPmSet(n,m,PA,aP,xstar,xbar,false,EPS_VAL)

    end



    function BlocksModelMM(file)
        M = MatrixMarket.mmread(file)
        m,n = size(M)
        sol = ones(n)
        b = M*sol
        xzero = StartingPoint(n)
        return M, b, xzero, sol
    end

    """
    circleShape(h::Float64,k::Float64,r::Float64)

    Returns the parametric function for a circle with center `(h,k)` and radius `r`.
    """
    function circleShape(h::Float64,k::Float64,r::Float64)
        θ = LinRange(0,2*π,500)
        return h .+ r*sin.(θ), k .+ r*cos.(θ)
    end





testconvex(5)

plot(circleShape(-.5,0.,1.), seriestype=:shape, lw = 1, c=:blue, linecolor=:black, legend = false, fillalpha = .4, aspectratio = 1,framestyle=:none)
# plot!(x->x,-1.5,2,lw=1,color=:black)
plot!(circleShape(.5,.0,1.), seriestype=:shape, lw = 1, c=:green, linecolor=:black, legend = false, fillalpha = .4, aspectratio = 1)
# X = readdlm("tables/xDR.dat")
# scatter!(X[:,1],X[:,2], color=:blue,markersize=3,lw = .5)
 # X = readdlm("tables/xCRM.dat")
 # scatter!(X[:,1],X[:,2], color=:red,markersize=3,lw = .5)


# MAPdata = readdlm("tables/N50_Samp100_1.0E-03-MAP.table");

# DRMdata = readdlm("tables/N50_Samp100_1.0E-03-DRM.table");

# DRMCdata = readdlm("tables/N50_Samp100_1.0E-03-DRM-C.table");
