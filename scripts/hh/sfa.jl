using Plots

const E_REST = -65.0
const C = 1.0
const G_LEAK = 0.3
const E_LEAK = 10.6 + E_REST
const G_NA = 120.0
const E_NA = 115.0 + E_REST
const G_K = 36.0
const E_K = -12.0 + E_REST

const TAU_AHP = 200
const G_AHP = 1400

const DT = 0.01
const T = 1000
const NT = 100000

alpha_m(v) = (2.5 - 0.1 * (v - E_REST)) / (exp.(2.5 - 0.1 * (v - E_REST)) - 1.0)
beta_m(v) = 4.0 * exp.(-(v - E_REST) / 18.0)
alpha_h(v) = 0.07 * exp.(-(v - E_REST) / 20.0)
beta_h(v) = 1.0 / (exp.(3.0 - 0.1 * (v - E_REST)) + 1.0)
alpha_n(v) = (0.1 - 0.01 * (v - E_REST)) / (exp.(1 - 0.1 * (v - E_REST)) - 1.0)
beta_n(v) = 0.125 * exp.(-(v - E_REST) / 80.0)

m0(v) = alpha_m(v) / (alpha_m(v) + beta_m(v))
h0(v) = alpha_h(v) / (alpha_h(v) + beta_h(v))
n0(v) = alpha_n(v) / (alpha_n(v) + beta_n(v))
tau_m(v) = 1.0 / (alpha_m(v) + beta_m(v))
tau_h(v) = 1.0 / (alpha_h(v) + beta_h(v))
tau_n(v) = 1.0 / (alpha_n(v) + beta_n(v))

dmdt(v, m) = (1.0 / tau_m(v)) * (-m + m0(v))
dhdt(v, h) = (1.0 / tau_h(v)) * (-h + h0(v))
dndt(v, n) = (1.0 / tau_n(v)) * (-n + n0(v))
dadt(s, a) = (1.0 / TAU_AHP) * (-a + s)
dvdt(v, m, h, n, a, i_ext) = (-G_LEAK * (v - E_LEAK)
                              -
                              G_NA * m * m * m * h * (v - E_NA)
                              -
                              G_K * n * n * n * n * (v - E_K)
                              -
                              G_AHP * a * (v - E_K)
                              +
                              i_ext) / C

function main()
    v = E_REST
    m = m0(v)
    h = h0(v)
    n = n0(v)
    a = 0
    v_old = E_REST
    i_ext = 40.0

    tl = []
    vl = []

    for nt in 0:NT-1
        t = DT * nt
        #println("$t $v $m $h $n")
        push!(tl, t)
        push!(vl, v)

        s = (v > 0 && v_old < 0)
        v_old = v

        dmdt1 = dmdt(v, m)
        dhdt1 = dhdt(v, h)
        dndt1 = dndt(v, n)
        dadt1 = dadt(s, a)
        dvdt1 = dvdt(v, m, h, n, a, i_ext)

        dmdt2 = dmdt(v + 0.5 * DT * dvdt1, m + 0.5 * DT * dmdt1)
        dhdt2 = dhdt(v + 0.5 * DT * dvdt1, h + 0.5 * DT * dhdt1)
        dndt2 = dndt(v + 0.5 * DT * dvdt1, n + 0.5 * DT * dndt1)
        dadt2 = dadt(s, a + 0.5 * DT * dadt1)
        dvdt2 = dvdt(v + 0.5 * DT * dvdt1, m + 0.5 * DT * dmdt1, h + 0.5 * DT * dhdt1, n + 0.5 * DT * dndt1, a + 0.5 * DT * dadt1, i_ext)

        dmdt3 = dmdt(v + 0.5 * DT * dvdt2, m + 0.5 * DT * dmdt2)
        dhdt3 = dhdt(v + 0.5 * DT * dvdt2, h + 0.5 * DT * dhdt2)
        dndt3 = dndt(v + 0.5 * DT * dvdt2, n + 0.5 * DT * dndt2)
        dadt3 = dadt(s, a + 0.5 * DT * dadt2)
        dvdt3 = dvdt(v + 0.5 * DT * dvdt2, m + 0.5 * DT * dmdt2, h + 0.5 * DT * dhdt2, n + 0.5 * DT * dndt2, a + 0.5 * DT * dadt2, i_ext)

        dmdt4 = dmdt(v + 0.5 * DT * dvdt3, m + 0.5 * DT * dmdt3)
        dhdt4 = dhdt(v + 0.5 * DT * dvdt3, h + 0.5 * DT * dhdt3)
        dndt4 = dndt(v + 0.5 * DT * dvdt3, n + 0.5 * DT * dndt3)
        dadt4 = dadt(s, a + 0.5 * DT * dadt3)
        dvdt4 = dvdt(v + 0.5 * DT * dvdt3, m + 0.5 * DT * dmdt3, h + 0.5 * DT * dhdt3, n + 0.5 * DT * dndt3, a + 0.5 * DT * dadt3, i_ext)

        m += DT * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6.0
        h += DT * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6.0
        n += DT * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6.0
        a += DT * (dadt1 + 2 * dadt2 + 2 * dadt3 + dadt4) / 6.0
        v += DT * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6.0

    end

    plot(tl, vl, fmt=:png, title=:"hh model(sfa 0~100ms)", xlabel=:"time[m/s]", ylabel="membrane potential[mV]", label=:"v", ylims=(-80, 60), xlims=(0, 100))
    png("figures/hh/sfa_0~100ms.png")
    plot(tl, vl, fmt=:png, title=:"hh model(sfa 900~1000ms)", xlabel=:"time[m/s]", ylabel="membrane potential[mV]", label=:"v", ylims=(-80, 60), xlims=(900, 1000))
    png("figures/hh/sfa_900~1000ms.png")


end

main()