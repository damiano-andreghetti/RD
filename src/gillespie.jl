@enum Event evdif evrea evatt evdet evcat evend

function build_queues(M)
    N = nsites(M)
    Qn = Tuple(StaticExponentialQueue(N) for _ in M.species)
    Qcat = Tuple(StaticExponentialQueue(N) for _ in M.cat)
    Qrea = Tuple(StaticExponentialQueue(N) for r in M.rea)
    Qatt = Tuple(Qn[m1]*Ref(0.0) for (m,m1,ka) in M.att)
    Q = NestedQueue((
        ((evdif,m) => Qn[m] * d for (m,d) in M.dif)...,
        ((evatt,m) => q for ((m,_,_),q) in zip(M.att, Qatt))...,
        ((evdet,m) => Qn[m] * kd for (m,kd) in M.det)...,
        ((evcat,r) => q*kc for (r,(_,_,_,kc,_),q) in zip(eachindex(M.cat),M.cat, Qcat))...,
        # reactions with only 1 substrate are treated differently because they don't need an extra queue
        ((evrea,r) => (length(s) == 1 ? Qn[only(s)] : q)*k for (r,(s,p,k),q) in zip(eachindex(M.rea), M.rea, Qrea))...,
       ))
    (Qn, Qcat, Qrea, Qatt, Q)
end


function run_RD!(state::State, M::Model, T; 
        stats = (_...,)->nothing, 
        rng = Random.default_rng(),
        queues = build_queues(M))

    Qn, Qcat, Qrea, Qatt, Q = queues

    function update(i::Int)
        for ((e,s,_,_,km),q) in zip(M.cat, Qcat)
            @inbounds q[i] = state.membrane[i,e]  / (1 + km/state.membrane[i,s])
        end
        for ((s,), q) in zip(M.rea, Qrea)
            if length(s) > 1
                @inbounds q[i] = prod(state.membrane[i,m] for m in s)
            end
        end
        for ((m,_,ka),q) in zip(M.att, Qatt)
            @inbounds q.f[] = state.cytosol[m]*ka
        end
        for (m,q) in pairs(Qn)
            @inbounds q[i] = state.membrane[i,m]
        end
    end

    foreach(update, 1:nsites(M))
            
    println("starting simulation, $(length(Q)) events in the queue")

    t::Float64 = 0.0
    while !isempty(Q)
        ((ev,m),i),dt = peek(Q; rng)
        t += dt
        t > T && break # reached end time of simulation
        stats(t, state, ev, m, i)
        @inbounds if ev == evdif # diffusion
            #arrival is chosen uniformly between its neighbours
            j = rand(rng, neighbors(M.G, i))
            state.membrane[i,m] -= 1
            state.membrane[j,m] += 1
            update(i)
            update(j)
        elseif ev == evcat # catalytic reaction
            _, s, p = M.cat[m]
            state.membrane[i,s] -= 1
            state.membrane[i,p] += 1
            update(i)
        elseif ev == evatt # attachment to membrane
            state.cytosol[m] -= 1
            state.membrane[i,m] += 1
            update(i)
        elseif ev == evdet # detachment from membrane
            state.membrane[i,m] -= 1
            state.cytosol[m] += 1
            update(i)
        else # ev == evrea # reaction
            s, p = M.rea[m]
            for m in s
                state.membrane[i,m] -= 1
            end
            for m in p
                state.membrane[i,m] += 1
            end
            update(i)
        end
    end
    stats(T, state, evend, 0, 0)
end
