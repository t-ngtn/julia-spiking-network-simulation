using Plots

const TAU = 20.0
const V_REST = -65.0
const V_RESET = -65.0
const THETA = -55.0
const R_M = 1.0
const DT = 1.0
const T = 1000.0
const NT = 1000
const I_EXT = 12.0
const N = 2

function main()
    v = [V_REST, V_REST - 15.0]
    s = [false, false]

    tl = [[], []]
    vl = [[], []]

    for nt in 0:NT-1
        t = DT * nt

        for i in 1:N
            push!(tl[i], t)
            push!(vl[i], v[i])
        end

        for i in 1:N
            v[i] += DT * (-(v[i] - V_REST) + R_M * I_EXT) / TAU
            s[i] = (v[i] > THETA)

            for j in 1:N
                if s[i]
                    push!(tl[i], t + DT)
                    push!(vl[i], v[i])
                    push!(tl[i], t + DT)
                    push!(vl[i], 0.0)
                end
            end
        end

        for i in 1:N
            v[i] = s[i] * V_RESET + (!s[i]) * v[i]
        end
    end

    plot(tl[1], vl[1], fmt=:png, title="lif model(2 neuron)", xlabel="time[m/s]", ylabel="membrane potential[mV]", label="v1", ylims=(-70, 0))
    plot!(tl[2], vl[2], label="v2")
    png("figures/lif/lif2.png")

end

main()