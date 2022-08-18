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

function main()
    v = V_REST

    tl = []
    vl = []

    for nt in 0:NT-1
        t = DT * nt

        push!(tl, t)
        push!(vl, v)

        v += DT * (-(v - V_REST) + R_M * I_EXT) / TAU
        s = (v > THETA)

        if s
            push!(tl, t + DT)
            push!(vl, v)
            push!(tl, t + DT)
            push!(vl, 0.0)
        end

        v = s * V_RESET + (!s) * v

    end

    plot(tl, vl, fmt=:png, title=:"lif model", xlabel=:"time[m/s]", ylabel="membrane potential[mV]", label=:"v", ylims=(-70, 0))
    png("figures/lif/lif.png")

end

main()