using Distributed
using DifferentialEquations
using DelimitedFiles
using BenchmarkTools

function Vern(A, B, C)
    tspan = (C[1], C[2])
    Tt = C[1]:C[3]:C[2]
    l1, l2, m1, m2, g = B[1], B[2], B[3], B[4], 10
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

function DIF(data)
    Td = [0.0]
    vals = data[2]
    timequt = data[5]
    Hu = [1]
    for i in 1:size(vals)[1]
        if vals[i] > vals[1] + 2.0 * pi || vals[i] < vals[1] - 2.0 * pi
            push!(Td, timequt[i])
            break
        end
    end
    return last(Td)
end

function fractline((th1), (th2), size, time, Jay)
    size = size - 1
    j = Jay
    Ua = zeros(1, size + 1)
    for i in 1:size+1
        Aj = [(pi / 180) * (th1[1] + ((th1[2] - th1[1]) / size) * (i - 1)), (pi / 180) * (th2[2] - ((th2[2] - th2[1]) / size) * (j - 1)), 0, 0]
        Bj = [1, 1, 1, 1]
        # Aa[j,i] = DIF(Vern(Aj,Bj,time),time)
        # Ua[1, i] = th1[1] + ((th1[2] - th1[1]) / size) * (i - 1) + (th2[2] - ((th2[2] - th2[1]) / size) * (j - 1))
        Ua[1, i] = DIF(Vern(Aj,Bj,time))
        # if i % 79 == 0 && j == 6
        #     print((pi / 180) * (th1[1] + ((th1[2] - th1[1]) / size) * (i)), (pi / 180) * (th2[2] - ((th2[2] - th2[1]) / size) * (j - 1)))
        # end
    end
    return Ua
end

function abc(amount)
    load,prec = 0, 100/amount

    Aa = zeros(amount, amount)
    Threads.@threads for i in 1:amount
        Aa[i, :] = fractline([-180,180], [-180,180], amount, [0, 100, 0.01], i)
        if i % Threads.nthreads() == 0
            load = load +1
            println(';',load*prec*Threads.nthreads())
        end
    end
    writedlm("Parfractdat.csv", Aa, ',')
    print("done")
end

@time abc(100)