using DelimitedFiles
using Plots
using LaTeXStrings
pgfplots()
# PGFPlots.pushPGFPlotsPreamble("\\usepackage{amssymb}")
# unicodeplots()
# pyplot()
# PyPlot.rc("text", usetex= true)
# PyPlot.rc("text.latex",preamble = ["\\usepackage{amsmath,amssymb}", "\\usepackage[utf8]{inputenc}"])
# gr()


"""
circleShape(h::Float64,k::Float64,r::Float64)

Returns the parametric function for a circle with center `v=(h,k)` and radius `r`.
"""
function circleShape(v,r::Float64)
    h = v[1]
    k = v[2]
    θ = range(0,2*π,length=500)
    return h .+ r*sin.(θ), k .+ r*cos.(θ)
end

v1 = [.8,0]
r1 = 1.
v2 = [-.8,0]
r2 = 1.

function plotBalls()
    plot(circleShape(v1,r1), seriestype=:shape, lw = 1, c=:blue, linecolor=:blue,  fillalpha = .4, aspectratio = 1,framestyle=:none,label="")
    # plot!(x->x,-1.5,2,lw=1,color=:black)
    plot!(circleShape(v2,r2), seriestype=:shape, lw = 1, c=:green, linecolor=:green,  fillalpha = .4, aspectratio = 1,label="")
    xMAP = readdlm("tables/xMAP-FigBalls1.dat")
    k = div(size(xMAP,1)+1,2) - 1
    scatter!(xMAP[:,1],xMAP[:,2], color=:green, line=(:path, .7), marker=(:square, Plots.stroke(.4, :black)), label="MAP - $k it")
    xDR = readdlm("tables/xDR-FigBalls1.dat")
    k = size(xDR,1)-1
    scatter!(xDR[:,1],xDR[:,2], color=:red,line=(:path, .7), marker=(5,:circle),label="DRM - $k it")
    XCRM = readdlm("tables/xCRM-FigBalls1.dat")
    k = size(XCRM,1) - 1
    scatter!(XCRM[:,1],XCRM[:,2], color=:blue, line=(:path, .7), marker=(5,:diamond),label="CRM - $k it")
    xCRMProd = readdlm("tables/xCRM-prod-FigBalls1.dat")
    k = size(xCRMProd,1) - 1
    scatter!(xCRMProd[:,1],xCRMProd[:,2], color=:black, line=(:path, .7), marker=(5,:utriangle),label="CRM-prod - $k it")

    annotate!([(-2.2,   -3.0,text(L"x^0",14))])
    annotate!([(.8,  0,text(L"X_2",16,:bottom, color=:blue))])
    annotate!([(-.8,  0,text(L"X_1",16,:bottom, color=:green))])

    # annotate!([(-4.814493449496354,  -6.713723121438827,text(L"x_2",14,:top))])
    # annotate!([(-14.199862805565246, 7.627976817750572,text(L"x_3",14,:bottom))])
    savefig("../../Draft/New/Two-Any-Convex-Sets/figures/Figure1_ComparisonDR-MAP-CRM-CRMProd1.pdf")
end
plotBalls()
# A = [1/4 1.]
# b = [-5.]
# xzero = [-16., -1]
# f(x) = (-A[1,1]*x + b[1])/A[1,2]
# plot(0:10,x->x,c=:red)
# plot!(0:10,x->-x,c=:red,aspect_ratio=:equal)
# scatter!(xzero[1:1],xzero[2:2])
# plot!(-21:20,f,framestyle=:origin,leg=false,aspect_ratio=:equal)
