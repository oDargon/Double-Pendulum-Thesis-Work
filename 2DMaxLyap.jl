using Distributed
using DifferentialEquations
using DelimitedFiles
using BenchmarkTools
using LinearAlgebra
using LaTeXStrings  
using Plots; pyplot()
using Statistics


function Vern(A, B, C)
    tspan = (C[1], C[2])
    Tt = C[1]:C[3]:C[2]
    l1, l2, m1, m2, g = B[1], B[2], B[3], B[4], 9.81
    u0 = [A[1]; A[2]; A[3]; A[4]]
    function Hamil!(du, u, p, t)
        th1, th2, mom1, mom2 = u
        h1 = (mom1 * mom2 * sin(th1 - th2)) / (l1 * l2 * (m1 + m2 * (sin(th1 - th2)^2)))
        h2 = (m2 * (l2^2) * (mom1^2) + (m1 + m2) * (l1^2) * (mom2^2) - 2 * m2 * l1 * l2 * mom1 * mom2 * cos(th1 - th2)) / (2 * (l1^2) * (l2^2) * (m1 + m2 * (sin(th1 - th2)^2))^2)
        du[1] = (l2 * mom1 - l1 * mom2 * cos(th1 - th2)) / ((l1^2) * l2 * (m1 + m2 * (sin(th1 - th2)^2)))
        du[2] = (-m2 * l2 * mom1 * cos(th1 - th2) + (m1 + m2) * l1 * mom2) / (m2 * l1 * (l2^2) * (m1 + m2 * sin(th1 - th2)^2))
        du[3] = -(m1 + m2) * g * l1 * sin(th1) - h1 + h2 * sin(2 * (th1 - th2))
        du[4] = -m2 * g * l2 * sin(th2) + h1 - h2 * sin(2 * (th1 - th2))
    end
    prob = ODEProblem(Hamil!, u0, tspan, maxiters=10^10)
    sol = solve(prob, Vern9(),reltol=1e-12,abstol=1e-12)
    Bb = sol(Tt)
    return [Bb[1,:],Bb[2,:],Bb[3,:],Bb[4,:],LinRange(C[1],C[2],trunc(Int,(1/C[3])*(C[2]-C[1])) +1)]
end

function Vern2(A, B, C)
    tspan = (C[1], C[2])
    l1, l2, m1, m2, g = B[1], B[2], B[3], B[4], 9.81
    u0 = [A[1]; A[2]; A[3]; A[4]]
    function Hamil!(du, u, p, t)
        th1, th2, mom1, mom2 = u
        h1 = (mom1 * mom2 * sin(th1 - th2)) / (l1 * l2 * (m1 + m2 * (sin(th1 - th2)^2)))
        h2 = (m2 * (l2^2) * (mom1^2) + (m1 + m2) * (l1^2) * (mom2^2) - 2 * m2 * l1 * l2 * mom1 * mom2 * cos(th1 - th2)) / (2 * (l1^2) * (l2^2) * (m1 + m2 * (sin(th1 - th2)^2))^2)
        du[1] = (l2 * mom1 - l1 * mom2 * cos(th1 - th2)) / ((l1^2) * l2 * (m1 + m2 * (sin(th1 - th2)^2)))
        du[2] = (-m2 * l2 * mom1 * cos(th1 - th2) + (m1 + m2) * l1 * mom2) / (m2 * l1 * (l2^2) * (m1 + m2 * sin(th1 - th2)^2))
        du[3] = -(m1 + m2) * g * l1 * sin(th1) - h1 + h2 * sin(2 * (th1 - th2))
        du[4] = -m2 * g * l2 * sin(th2) + h1 - h2 * sin(2 * (th1 - th2))
    end
    prob = ODEProblem(Hamil!, u0, tspan, maxiters=10^10)
    sol = solve(prob, Vern9(),reltol=1e-9,abstol=1e-9)
    return [sol[1, :], sol[2, :], sol[3, :], sol[4, :], sol.t]
end

function Ojo(Z,Intcord, Intconst, tspan)
    O = Vern(Intcord, Intconst, tspan)
    th1, th2, mom1, mom2, tims = O[1], O[2], O[3], O[4], O[5]
    qu = 1e-5
    vec = Z*(qu/norm(Z))
    for i in 2:size(th1)[1]
        F = Vern2([th1[i-1], th2[i-1], mom1[i-1], mom2[i-1]]+vec,Intconst,[tims[i-1],tims[i]])
        nvec = [last(F[1])-th1[i],last(F[2])-th2[i],last(F[3])-mom1[i],last(F[4])-mom2[i]]
        vec = nvec*(qu/norm(nvec))
    end 
    return vec
end

function NwLyp(Z,Intcord, Intconst, tspan)
    O = Vern(Intcord, Intconst, tspan)
    th1, th2, mom1, mom2, tims = O[1], O[2], O[3], O[4], O[5]
    qu = 1e-5
    vectt = Z*(qu/norm(Z))
    vec = Ojo(vectt,Intcord,Intconst,[0,100,0.25])
    d0 = Vector{Float64}[]
    d1 = Vector{Float64}[]
    sol = []
    for i in 2:size(th1)[1]
        push!(d0, vec)
        F = Vern2([th1[i-1], th2[i-1], mom1[i-1], mom2[i-1]]+vec,Intconst,[tims[i-1],tims[i]])
        nvec = [last(F[1])-th1[i],last(F[2])-th2[i],last(F[3])-mom1[i],last(F[4])-mom2[i]]
        push!(d1, nvec)
        push!(sol,(log(norm(d1[i-1]/d0[i-1])))/(tims[i]-tims[i-1]))
        vec = nvec*(qu/norm(nvec))
    end
    return mean(sol)
end

function tline(size,th1s,th1e,th2)
    Ua = zeros(1, size)
    angs = []
    for i in 1:size
        Aj = [(pi / 180) * (th1s + ((th1e-th1s)/size)*i), (pi / 180) * th2, 0, 0]
        Bj = [1, 1, 1, 1]
        Ua[1, i] = NwLyp([0,0,10,0],Aj,Bj,[0,500,0.25])
        push!(angs,th1s + ((th1e-th1s)/size)*i)
    end
    return Ua
end

function tline2(size,th1s,th1e)
    Ua = zeros(1, size)
    angs = []
    for i in 1:size
        Aj = [(pi / 180) * (th1s + ((th1e-th1s)/size)*i), (pi / 180) * (th1s + ((th1e-th1s)/size)*i), 0, 0]
        Bj = [1, 1, 1, 1]
        Ua[1, i] = NwLyp([0,0,10,0],Aj,Bj,[0,500,0.25])
        push!(angs,th1s + ((th1e-th1s)/size)*i)
    end
    return angs,Ua'
end

function tline3(size,th2s,th2e,th1)
    Ua = zeros(1, size)
    angs = []
    for i in 1:size
        Aj = [(pi / 180) * (th1), (pi / 180) * (th2s + ((th2e-th2s)/size)*i), 0, 0]
        Bj = [1, 1, 1, 1]
        Ua[1, i] = NwLyp([0,0,10,0],Aj,Bj,[0,500,0.25])
        push!(angs,th2s + ((th2e-th2s)/size)*i)
    end
    return angs,Ua'
end

function tplot(size,th1s,th1e,th2s,th2e)
    Aa = zeros(size, size)
    for j in 1:size
        Aa[j, :] = tline(size,th1s,th1e,th2s + ((th2e-th2s)/size)*j)
        print(100*(1/size)*j)
    end
    print("done")
    writedlm("3DLyap.csv", Aa, ',')
    return range(th1s,stop=th1e,length=size), range(th2s,stop=th2e,length=size), Aa
end

# F = tline(10,100,130,90.5)
# plot(F)   

# @btime NwLyp([0,0,10,0],[(pi/180)*30,(pi/180)*30,0,0],[1,1,1,1],[0,100,0.25])

F = tplot(500,100,130,70,100)
# plot(F,st=:surface,camera=(-30,30))
heatmap(F[1],F[2],F[3])


# F = tline(1000,0,180,0)
# F2 = tline2(1000,0,180) 
# plot(F[1],[F[2],F2[2]],xlab = L"\theta_1(0)",ylab = L"\lambda_{max} (s^{-1})",labels= [L"\theta_1 = \theta_2 " L"\theta_2 = 0"],legend=:topleft,grid = false,legendfontsize=10,yguidefontsize=13,xguidefontsize=13)

# G = tline3(1000,0,180,100)

# plot(G,xlab = L"\theta_2(0)",ylab = L"\lambda_{max} (s^{-1})",labels=  L"\theta_1 = 100",legend=:topleft,grid = false,legendfontsize=10,yguidefontsize=13,xguidefontsize=13)