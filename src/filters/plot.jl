using ColorVectorSpace, Images, Colors
import Compose: polygon, context, fill, mm, set_default_graphic_size, compose, draw, PNG, PDF
import Cairo, Fontconfig

function hexagon(x, y, r)
    polygon([[(x[i]+r*cos(π/3*(k+3/2)),
               y[i]+r*sin(π/3*(k+3/2))) for k in 0:6] 
                    for i in eachindex(x,y)])
end

square(x, y, r) = rectangle(x.-r,y.-r,fill(2r,length(x)),fill(2r,length(x)))

"""
`composed(s::State; layout, colors, Δ)`

Generates a `Compose` image of a lattice of a given `layout` from `s::State` 
"""
function composed(s::State; layout, colors, Δ)
    (cellshape, posx, posy) = layout
    X, Y = Δ .* posx, Δ .* posy
    (x0,x1),(y0,y1) = extrema(X), extrema(Y)
    set_default_graphic_size(x1-x0+3Δ, y1-y0+3Δ)
    if cellshape == :hexagon
        compose(context(), 
            hexagon(X,Y,Δ),
            fill(s.membrane * colors)
        )
    elseif cellshape == :square
        compose(context(),
            square(X,Y,Δ),
            fill(s.membrane * colors)
        )
    else
        throw(ArgumentError("Unrecognized shape $cellshape"))
    end 
end

"""
`Plotter(posx, posy; colors, Δ=mm)`

A filter to generate a `Compose` image from each state `s`. Use e.g. the
`display ∘ Plotter(posx, posy; colors)` filter for iterative display, etc.
"""
function Plotter(layout; colors, Δ=mm)
    f(s::State, args...) = composed(s; layout, colors, Δ)
    f(_, s::State, args...) = f(s)
    f((_,s,args...)) = f(s)
    f
end

"""
`raster(c::Compose.Context) -> Matrix{<:RGB}`

Converts a `Context` object into a raster image.
"""
function raster(c)
    io = IOBuffer()
    draw(PNG(io; emit_on_finish=false), c)
    A = Images.load(io)
    m, n = 2 .* (size(A) .÷ 2)
    #crop to even dimensions and N0f8 colorspace
    RGB{N0f8}.(@view A[1:m,1:n])
end


function savepdf(filename, p)
    open(filename, "w") do io 
        draw(PDF(io), p); 
        close(io)
    end
end