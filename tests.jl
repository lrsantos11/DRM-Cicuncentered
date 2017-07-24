workspace()
include("MAP.jl")
include("mSet.jl")
# using JLD
importall MAP
importall mSet
global const EPS_VAL = 1e-6
global const ZERO_VAL = 1e-15

    function GeneralTest(n::Int64 = 100,samples::Int64 = 2,restarts::Int64=1, affine::Bool=true)
        # Fix Random
        srand(2)
        # X = R^n
        j = 1
        firstwrite = true
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
            cF = FriedrichsAngleAB(A,B)
            AngleAB = acos(cF)
            # println(AngleAB)
            if AngleAB > 1e-2
                # println("Find better AngleAB")
               continue
            end
            P, p = contructProjector([A; B],[a; b],n)
            for i = 1:restarts
                xzero = StartingPoint(n)
                xbar = Projection(P,p,xzero)
                prob_name  = String("Problem$j""-$i""A$ma""B$mb")
                #########################
                # Saving data to file
                # jldopen("tables/highcF.jld", "r+") do file
                #     write(file, prob_name,CC)  
                # end

                #########################
                println("MAP with no acceleration")
                fname = @sprintf("tables/N%d_Samp%d_%.1E-MAP-true.table",n,samples,EPS_VAL);
                open(fname,"a") do file
                    if firstwrite
                        str = @sprintf("---\nalgname: MAP\nsuccess: conv\nfree_format: True\n---\n")
                        print(file,str)
                    end
                    write(file,prob_name)
                    time = @elapsed algorithmMAP(A,a,B,b,n,xzero,xbar,file,true,1,EPS_VAL)
                    str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                #     # print(str)
                    print(file,str)
                end
                #########################
                # println("MAP with CP acceleration")
                # @time MAP.algorithmMAP(A,a,B,b,n,xzero,2)
                ########################
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
                # @time MAP.fourpointscheme(A,a,B,b,n,xzero,2)
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
                println("DRM Circumcentered algorithm")
                fname = @sprintf("tables/N%d_Samp%d_%.1E-C-DRM-true.table",n,samples,EPS_VAL);
                open(fname,"a") do file
                    if firstwrite
                        str = @sprintf("---\nalgname: C-DRM\nsuccess: conv\nfree_format: True\n---\n")
                        print(file,str)
                    end
                    write(file,prob_name)
                    time = @elapsed algorithmDRM_C(A,a,B,b,n,xzero,xbar,file,true,EPS_VAL)
                    str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)
                    # print(str)
                    print(file,str)
                end
                #####################
                println("DRM algorithm")
                fname = @sprintf("tables/N%d_Samp%d_%.1E-DRM-true.table",n,samples,EPS_VAL);
                open(fname,"a") do file
                    if firstwrite
                        str = @sprintf("---\nalgname: DRM\nsuccess: conv\nfree_format: True\n---\n")
                        print(file,str)
                    end
                    write(file,prob_name)
                    time = @elapsed algorithmDRM(A,a,B,b,n,xzero,xbar,file,true,EPS_VAL)
                    str = @sprintf(" %10.8e %10.8f %10.8f\n",time,AngleAB,cF)                    # print(str)
                    print(file,str)
                end
                firstwrite = false
            end
            j += 1
        end
        
    end


    function ReflecA(x)
        # A is the Ball centered in v and ray r
            # if norm(x)<=1 
            #     return x
            # end
            # v=[-sqrt(2)/2;sqrt(2)/2]
            v=[-.5;.5]
            r = 1.0
            # proj = v + (x-v)/(max(1,norm(x-v)/r))
            proj = v + (x-v)/(norm(x-v)/r)
            return 2.*proj - x
    end
    # function ReflecB(x)
    #     # A is the Ball centered in v
    #     #     if norm(x)<=1 
    #     #         return x
    #     #     end
    #         v=[1.5;0.]
    #         r = 1.
    #         proj = v + (x-v)/(max(1,norm(x-v)/r))
    #         return 2.*proj - x
    # end
    function ReflecB(x)
        # B is the  (affine) subspace that line that crosses [1;1]
            n = length(x)
            OneVec = ones(n)
            proj = dot(OneVec,x)/dot(OneVec,OneVec)*OneVec
            return 2.*proj - x
    end
    function testeconvex(itmax)
        if isfile("tables/xDR.dat")
            rm("tables/xDR.dat")
            rm("tables/xDR-C.dat")
        end    
        xDRC = [-0.06578485818175905 ;   2.1374332764804542] #5.*randn(2) #[1.233205883326263 ; .9317100309585165]
        OneVec = ones(2)    
        xDR = xDRC
        xbar = [0.; 0.]
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

        println("Projections of  DR: $k")
        println("Norm for xDR: $tolDR")
        k = 0
        open("tables/xDR-C.dat","a") do f
           # do stuff with the open file
          writedlm(f,xDRC')
        end
        while tolC > EPS_VAL && k <= itmax
            # DR-C
            xDRCold = xDRC    
            ypto = ReflecA(xDRC)
            zpto = ReflecB(ypto)
            if norm(ypto - xDRC)<ZERO_VAL
                xDRC = 0.5*(xDRC + zpto)
                # tolC = norm(xDRC - xbar,2)
                tolC = norm(xDRC - xDRCold,2)
                k += 1
                open("tables/xDR-C.dat","a") do f
                   # do stuff with the open file
                  writedlm(f,xDRC')
                end 
                continue
            elseif norm(zpto - ypto)<ZERO_VAL
                xDRC = 0.5*(xDRC + zpto)
                # tolC = norm(xDRC - xbar,2)
                tolC = norm(xDRC - xDRCold,2)
                k += 1
                open("tables/xDR-C.dat","a") do f
                   # do stuff with the open file
                  writedlm(f,xDRC')
                end 
                continue
            end
            xDRC = MAP.FindCircumcenter(xDRC,ypto,zpto)
            # println("xpto: $xDRC")
            # tolC = norm(xDRC - xbar,2)
            tolC = norm(xDRC - xDRCold,2)
            k += 1
            open("tables/xDR-C.dat","a") do f
               # do stuff with the open file
              writedlm(f,xDRC')
            end  
        end
        # println("xC: $xDRC")
        println("Projections of DR-C: $k")
        println("Norm for xC: $tolC")
    end


    function GeneralTestmSet(n::Int64 = 100,m::Int64 = 2,affine::Bool=false)
        # Fix Random
        # srand(30)
        xstar = MAP.StartingPoint(n)
        # println(xstar)
        AA, R = GenerateSampleApart(n,m)
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

        time = @elapsed mSet.algorithmCDRMmSet(n,m,PA,aP,xstar,xbar,false,EPS_VAL) 
        time = @elapsed mSet.algorithmMAPmSet(n,m,PA,aP,xstar,xbar,false,EPS_VAL) 

    end



# testeconvex(1000)

# GeneralTest(200,25,1,false)

  GeneralTestmSet(50,10,false)

# MAPdata = readdlm("tables/N50_Samp100_1.0E-03-MAP.table");

# DRMdata = readdlm("tables/N50_Samp100_1.0E-03-DRM.table");

# DRMCdata = readdlm("tables/N50_Samp100_1.0E-03-DRM-C.table");
