struct State{M,V}
    membrane::M
    cytosol::V
end

function State(M::Model, totmembrane::Vector, cytosol::Vector; rng = Random.default_rng())
    Nspecies, Nsites = nspecies(M), nsites(M)
    membrane = fill!(FixedSizeMatrixDefault{Int}(undef, Nspecies, Nsites), 0) |> transpose
    for m in 1:Nspecies
        for _ in 1:totmembrane[m]
            membrane[rand(rng, 1:Nsites), m] += 1
        end
    end
    State(membrane, cytosol)
end
