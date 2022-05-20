using DifferentialEquations
using BenchmarkTools


function Vern(A, B, C)
    tspan = (C[1],C[2])
    Tt = C[1]:C[3]:C[2]
    l1, l2, m1, m2, g = B[1], B[2], B[3], B[4], 9.8
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
    #= sol = solve(prob,Vern9(), reltol=1e-9,abstol=1e-9) =#
    sol = solve(prob,Vern9(),reltol=1e-9,abstol=1e-9)
    # [sol(Tt)[1,:],sol(Tt)[2,:],sol(Tt)[3,:],sol(Tt)[4,:],Tt]
    # [Bb[1,:],Bb[2,:],Bb[3,:],Bb[4,:]]
    Bb = sol(Tt)
    return [Bb[1,:],Bb[2,:],Bb[3,:],Bb[4,:]]
end



# @profview Vern([(pi/180)*140,(pi/180)*146,0,0],[1,1,1,1],[0,1,0.01])  # run once to trigger compilation (ignore this one)
