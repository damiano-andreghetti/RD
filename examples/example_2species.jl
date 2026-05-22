using MembraneRD, Random, Colors, MembraneRD.Filters, JLD2

species = @species A B EA EB

function build_model(L)
    G,layout = MembraneRD.gen_hex_lattice(L)
    N = nv(G)
    unit_length = 1.0
    theta = 0.2126944621086619 #Stot/V
    Ac = unit_length^2 #area associated to each cell, set unit lengthscale
    Stot = N*Ac #area of each cell is assumed to be 1, setting unit lengthscale
    V = Stot/theta
    dA = dB = dEA = dEB = 0.01 #this sets unit timescale
    #rates in the theory (do not correspond excatly to those used in the simulations), some slight dimensional changes are needed
    kAa_th = 1.0; kAd_th = 1.0*kAa_th; kAc_th = 1.0
    kBa_th = 1.0; kBd_th = 1.0*kBa_th; kBc_th = 1.3*kAc_th
    KMMA_th = 1.0; KMMB_th = 1.0
    #rates to implement
    kAc, kAd, kAa = kAc_th, kAd_th, kAa_th / V
    kBc, kBd, kBa = kBc_th, kBd_th, kBa_th / V
    KMMA, KMMB = KMMA_th*Ac, KMMB_th*Ac
    kAb, kBb = 1e-4, 1e-4


    cat = ((@catalytic (kAc,KMMA) EA+B => EA+A), 
           (@catalytic (kBc,KMMB) EB+A => EB+B))
    att = ((EA,A,kAa), (EB,B,kBa))
    det = ((EA,kAd),(EB,kBd))
    dif = ((A,dA),(B,dB),(EA,dEA),(EB,dEB))
    rea = ((A,),(B,),kAb), ((B,),(A,),kBb)

    M = Model(; species, G, cat, att, det, dif, rea, rho_0 = 0.0)

    # create initial state distribution
    totmol = 10N
    totA, totB = floor(Int, 0.5totmol), floor(Int, 0.5totmol)
    totEA, totEB = floor(Int, 0.1N), floor(Int, 0.1N)
    memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+theta/(kAd_th/kAa_th)*totA/Stot))
    memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+theta/(kBd_th/kBa_th)*totB/Stot))
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    (; M, mem = [totA, totB, memEA, memEB], cyto = [0.0, 0.0, cytoEA, cytoEB], layout)
end



# structure to hold measures
struct Measure
	time::Float64
    phi_av::Float64
    bc::Float64
    cytoEA::Float64
    cytoEB::Float64
    rho::Float64
    function Measure(M::Model, s::State, t)
        Nc = nsites(M)
        #calculate phi at each cell (normalized, goes from -1, to 1)
        ϕ(i) = (s.membrane[i,B] - s.membrane[i,A])/(s.membrane[i,A] + s.membrane[i,B] + 1e-10)
        @views totA, totB = sum(s.membrane[:,A]), sum(s.membrane[:,B])
        #measure phi
        ϕav = (totB - totA) / (totA + totB)
        #in some denominators I added + 10^-10 in order to avoid NaNs
        m2 = sum((ϕ(i)-ϕav)^2 for i in 1:Nc) / Nc
        m4 = sum((ϕ(i)-ϕav)^4 for i in 1:Nc) / Nc
        #measure Binder cumulant
        bc = 1 - (m4 / (m2 ^ 2 + 1e-10)) / 3			
        #the following is to check for convergence of rho to 1
        rho = M.rho_0 * ((M.det[1][2] / M.att[1][3]) + totA) / ((M.det[2][2] / M.att[2][3]) + totB)
        new(t, ϕav, bc, s.cytosol[EA], s.cytosol[EB], rho)
    end
end

# save measures and print summary
function Saver(M::Model; name)
    function save_measure(t, s, _...)
        m = Measure(M, s, t)
        println()
        println("T = $t and <ϕ>/c = $(m.phi_av)")
        println("Binder cumulant: $(m.bc)")
        save("$name/config_T=$(round(t,digits=3)).jld",compress=true, "state", s)
        save("$name/measure_T=$(round(t,digits=3)).jld",compress=true, "measure", m)
    end
end



T = 20000.0
L = 100

# for reproducibility
rng = Random.Xoshiro(22)
(; M, mem, cyto, layout) = build_model(L)
s = State(M, mem, cyto; rng)

colors = ([colorant"yellow",colorant"blue",colorant"black",colorant"black"])/30
# push measures into a stack, access them with measurer.stack
measurer = Pusher((t,s,_...)->Measure(M,s,t), Measure)
# push states into a stack, access them with states.stack
states = Pusher((t,s,_...)->deepcopy(s), State)
# save every 500 time units into a file (and show a summary)
saver = TimeFilter(Saver(M; name="test_example_1"); times = 0:500:T)
# display the membrane once every 1s if the output is graphics-capable, e.g. in VSCode.
displayer = StopWatchFilter(display ∘ Plotter(layout; colors); seconds=2)
# for a text-only terminal, look at Sixel and use something like: displayer = StopWatchFilter((x->(println(); println(sixel_encode(x)))) ∘ raster ∘ Plotter(layout; colors); seconds=1.0)

# put everything together
stats = TimeFilter(ProgressShower(T), 
#    measurer,
    states,
    displayer,
#    saver, 
    ; times = 0:100:T)

# run simulation
@time run_RD!(s, M, T; stats, rng)


# to generate video afterwards, use this:
# savevideo("video.mp4", Iterators.map(raster∘Plotter(layout; colors), states.stack))
