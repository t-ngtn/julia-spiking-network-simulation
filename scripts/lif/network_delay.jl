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
const TAU_SYN = 5.0
const R_SYN = 1.0
const W = 2.0  # 同期状態　逆相同期にするには符号を反転
const T_REFR = 5.0
const NT_REFR = 5
const DELAY_SYN = 2.0
const NDELAY_SYN = 2

function main()
    v = [V_REST, V_REST - 15.0]  # 同期がわかりやすくしている 逆相同期の場合は0.1とか2とかにする
    i_syn = [0.0, 0.0]
    s = [false, false]
    ts = [0, 0]
    refr = [0, 0]

    tl = [[], []]
    vl = [[], []]

    for nt in 0:NT-1
        t = DT * nt

        for i in 1:N
            push!(tl[i], t)
            push!(vl[i], v[i])
        end

        for i in 1:N
            i_syn[i] = exp.(-DT / TAU_SYN) * i_syn[i] + W * (ts[i == 1 ? 2 : 1] + NDELAY_SYN == nt)
        end
        for i in 1:N
            v[i] += DT * (-(v[i] - V_REST) + R_SYN * i_syn[i] + R_M * I_EXT) / TAU
            s[i] = (v[i] > THETA)
            ts[i] = s[i] * (nt + 1) + !s[i] * ts[i]
        end

        for i in 1:N
            if s[i]
                push!(tl[i], t + DT)
                push!(vl[i], v[i])
                push!(tl[i], t + DT)
                push!(vl[i], 0.0)
            end
        end

        for i in 1:N
            refr[i] = s[i] * NT_REFR + !s[i] * (refr[i] - 1)
            v[i] = refr[i] > 0 ? V_RESET : v[i]
        end
    end

    plot(tl[1], vl[1], fmt=:png, title="lif model(2 neuron network delay)", xlabel="time[m/s]", ylabel="membrane potential[mV]", label="v1", ylims=(-70, 0))
    plot!(tl[2], vl[2], label="v2")
    png("figures/lif/network_delay.png")

end

main()