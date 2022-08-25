using Plots
const N = 4000
const N_E = 3200
const N_I = N - N_E

const T = 1000
const DT = 1.0
const NT = 1000

const TAU_M = 20.0
const TAU_E = 5.0
const TAU_I = 10.0

const V_REST = -49.0
const V_INIT = -60.0
const V_RESET = -60.0
const THETA = -50.0

const G_E = 1.62 / TAU_E
const G_I = -9.0 / TAU_I
const P = 0.02

v = Array{Float64}(undef, N)
ge = Array{Float64}(undef, N)
gi = Array{Float64}(undef, N)
w = Array{Float64}(undef, N * N)
s = Array{Bool}(undef, N)

anst = []
ansi = []

function initialize()
    for i in 1:N
        v[i] = V_INIT + 10.0 * rand(Float16)
    end

    for i in 1:N
        for j in 1:N
            w[j+N*(i-1)] = rand(Float16) < P ? 1.0 : 0.0
        end
    end
end

function calculate_synptic_inputs()
    for i in 1:N
        re = 0
        ri = 0
        for j in 1:N
            r = w[j+N*(i-1)] * s[j]
            if j <= N_E
                re += r
            else
                ri += r
            end
        end
        ge[i] = exp.(-DT / TAU_E) * ge[i] + re
        gi[i] = exp.(-DT / TAU_I) * gi[i] + ri
    end
end

function update_cell_parameters()
    for i in 1:N
        v[i] += DT * (-(v[i] - V_REST) + G_E * ge[i] + G_I * gi[i]) / TAU_M
        s[i] = (v[i] > THETA)
        v[i] = s[i] * V_RESET + !s[i] * v[i]
    end
end

function push_spike(nt)
    for i in 1:N
        if s[i]
            push!(anst, DT * (nt + 1))
            push!(ansi, i)
        end
    end

end

function loop()
    for nt in 0:NT-1
        calculate_synptic_inputs()
        update_cell_parameters()
        push_spike(nt)
        println("$nt")
    end
end

function main()
    initialize()
    loop()
    plot(anst, ansi, fmt=:png, seriestype=:scatter, markersize=0.05, title="random network raster", xlabel="time[m/s]", ylabel="neuron number", label=false)
    png("figures/lif/random.png")
end

main()