using Distributed
using DifferentialEquations
using DelimitedFiles
using BenchmarkTools
using LinearAlgebra
using LaTeXStrings  
using Plots; pyplot()
using Statistics




function Energy(th1,th2,mom1,mom2)
    m1, m2, l1, l2, g = 1,1,1,1,10
    return (m2*(l2^2)*(mom1^2) + (m1+m2)*(l1^2)*(mom2^2)-2*m2*l1*l2*mom1*mom2*cos(th1-th2))/(2*m2*(l1^2)*(l2^2)*(m1+m2*(sin(th1-th2)^2))) - (m1+m2)*g*l1*cos(th1) - m2*g*l2*cos(th2)
end

function compr(one, two)
    if one > two 
        return 1
    else
        return 0
    end
end


function tline(size,th1s,th1e,th2)
    yuan = Energy((pi/180)*0, (pi/180)*180, 0,0)
    Ua = zeros(1, size)
    for i in 1:size
        Aj = [(pi / 180) * (th1s + ((th1e-th1s)/size)*i), (pi / 180) * th2, 0, 0]
        Ua[1, i] = compr(yuan, Energy(Aj[1],Aj[2],Aj[3],Aj[4]))
    end
    return Ua
end

function tplot(size,th1s,th1e,th2s,th2e)
    Aa = zeros(size, size)
    for j in 1:size
        Aa[j, :] = tline(size,th1s,th1e,th2s + ((th2e-th2s)/size)*j)
        # print(100*(1/size)*j)
    end
    writedlm("EngBound.csv", Aa, ',')
    return range(th1s,stop=th1e,length=size), range(th2s,stop=th2e,length=size), Aa
end

F = tplot(3000,-180,180,-180,180)
# plot(F,st=:surface,camera=(-30,30))
heatmap(F[1],F[2],F[3])