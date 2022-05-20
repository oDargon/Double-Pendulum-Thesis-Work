using DifferentialEquations
using BenchmarkTools
using Plots
using LaTeXStrings  
pyplot()


function Vern(A, B, C)
    tspan = (C[1],C[2])
    Tt = C[1]:C[3]:C[2]
    l1, l2, m1, m2, g = B[1], B[2], B[3], B[4], 10
    u0 = [A[1];A[2];A[3];A[4]]
    function Hamil!(du,u,p,t)
        th1, th2, mom1, mom2 = u

        h1 = (mom1*mom2*sin(th1-th2))/(l1*l2*(m1+m2*(sin(th1-th2)^2)))
        h2 = (m2*(l2^2)*(mom1^2) +(m1+m2)*(l1^2)*(mom2^2)-2*m2*l1*l2*mom1*mom2*cos(th1-th2))/(2*(l1^2)*(l2^2)*(m1+m2*(sin(th1-th2)^2))^2)

        du[1] = (l2*mom1 -l1*mom2*cos(th1-th2))/((l1^2)*l2*(m1 + m2*(sin(th1-th2)^2)))

        du[2] = (-m2*l2*mom1*cos(th1-th2) + (m1+m2)*l1*mom2)/(m2*l1*(l2^2)*(m1 +m2*sin(th1 -th2)^2))

        du[3] = -(m1 +m2)*g*l1*sin(th1) -h1 + h2*sin(2*(th1-th2))

        du[4] = -m2*g*l2*sin(th2) +h1 -h2*sin(2*(th1-th2))
    end
    prob = ODEProblem(Hamil!,u0,tspan)
    sol = solve(prob,RK4();dt = 0.0001,adaptive=false)
    # sol = solve(prob,RK4();dt = 0.0001,adaptive=false)
    # [sol(Tt)[1,:],sol(Tt)[2,:],sol(Tt)[3,:],sol(Tt)[4,:],Tt]
    # [Bb[1,:],Bb[2,:],Bb[3,:],Bb[4,:]]
    return [sol[1,:],sol[2,:],sol[3,:],sol[4,:],sol.t]
end



function EPlot(sol,A,B)
    l1, l2, m1, m2, g = B[1], B[2], B[3], B[4], 10
    th1i, th2i, mom1i, mom2i = A[1], A[2], A[3], A[4]
    function Energy(th1,th2,mom1,mom2)
        I = (m2*(l2^2)*(mom1i^2) + (m1+m2)*(l1^2)*(mom2i^2)-2*m2*l1*l2*mom1i*mom2i*cos(th1i-th2i))/(2*m2*(l1^2)*(l2^2)*(m1+m2*(sin(th1i-th2i)^2))) - (m1+m2)*g*l1*cos(th1i) - m2*g*l2*cos(th2i)
        return 100*((m2*(l2^2)*(mom1^2) + (m1+m2)*(l1^2)*(mom2^2)-2*m2*l1*l2*mom1*mom2*cos(th1-th2))/(2*m2*(l1^2)*(l2^2)*(m1+m2*(sin(th1-th2)^2))) - (m1+m2)*g*l1*cos(th1) - m2*g*l2*cos(th2) -I)/I
    end
    time = sol[5]
    energies = abs.(Energy.(sol[1],sol[2],sol[3],sol[4]))
    t2 = []
    en2 = []
    for i in 1:size(sol[5])[1]
        if energies[i] == 0
            nothing
        else
            push!(t2,time[i])
            push!(en2,energies[i])
        end
    end
    return t2,en2
end

# function EnergySpec(sol,cun)
#     H = []
#     for i in 1:size(sol[1])[1]
#         push!(H,Energy((sol[1])[i],(sol[2])[i],(sol[3])[i],(sol[4])[i],cun))
#     end
#     return sol[5],H
# end

Ytick = [1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3]

A1 = Vern([(pi/180)*50,(pi/180)*50,0,0],[1,1,1,1],[0,1000,0.01])
Q1 = EPlot(A1,[(pi/180)*50,(pi/180)*50,0,0],[1,1,1,1])
print("hi")
A2 = Vern([(pi/180)*80,(pi/180)*80,0,0],[1,1,1,1],[0,1000,0.01])
Q2 = EPlot(A2,[(pi/180)*80,(pi/180)*80,0,0],[1,1,1,1])
print("hi")
A3 = Vern([(pi/180)*110,(pi/180)*110,0,0],[1,1,1,1],[0,1000,0.01])
Q3 = EPlot(A3,[(pi/180)*110,(pi/180)*110,0,0],[1,1,1,1])
print("hi")
A4 = Vern([(pi/180)*140,(pi/180)*140,0,0],[1,1,1,1],[0,1000,0.01])
Q4 = EPlot(A4,[(pi/180)*140,(pi/180)*140,0,0],[1,1,1,1])
print("hi")
A5 = Vern([(pi/180)*170,(pi/180)*170,0,0],[1,1,1,1],[0,1000,0.01])
Q5 = EPlot(A5,[(pi/180)*170,(pi/180)*170,0,0],[1,1,1,1])
print("hi")
A6 = Vern([(pi/180)*20,(pi/180)*20,0,0],[1,1,1,1],[0,1000,0.01])
Q6 = EPlot(A6,[(pi/180)*20,(pi/180)*20,0,0],[1,1,1,1])
print("hi")
p1 = plot([Q6[1],Q1[1],Q2[1],Q3[1],Q4[1],Q5[1]],[Q6[2],Q1[2],Q2[2],Q3[2],Q4[2],Q5[2]],yaxis=:log10, xlabel="time (s)", ylabel="Energy Drift (%)", labels=[ L"\theta_{1} = \theta_{2} = 20" L"\theta_{1} = \theta_{2} = 50" L"\theta_{1} = \theta_{2} = 80" L"\theta_{1} = \theta_{2} = 110" L"\theta_{1} = \theta_{2} = 140" L"\theta_{1} = \theta_{2} = 170"], yticks=Ytick,legend=:bottomright,legendfontsize =  9, reuse = false, linewidth=1, grid = false, tickfontsize = 10, guidefontsize = 12)
display(p1)