using RecipesBase
using ColorTypes
import PlotUtils: cgrad

const mygreen = RGBA{Float64}(151/255,180/255,118/255,1)
const mygreen2 = RGBA{Float64}(113/255,161/255,103/255,1)
const myblue = RGBA{Float64}(74/255,144/255,226/255,1)

@userplot Streamlines

@recipe function f(s::Streamlines)


    #if (isa(elements, Array) || isa(elements, Tuple)) &&
    #    ConformalBody in typeof.(elements)
    if (length(s.args)==4)

        ζ = s.args[3]
        elements = s.args[4]

        ψ = streamfunction(ζ,elements)

        #ζ = [r*exp(im*θ) for θ in s.args[2], r in s.args[1]]
        #b = elements[findfirst(isa.(elements,PotentialFlow.ConformalBody))]
        #Z = conftransform(ζ,b)

        @series begin
          seriestype --> :contour
            grid --> :none
            aspect_ratio --> 1
            linewidth --> 1
            legend --> :none
            seriescolor --> [:black, :black]
            #levels --> v
            #real.(Z), imag.(Z), ψ
            s.args[1], s.args[2], ψ
        end

    else
      elements = s.args[3]

      z = [x + im*y for y in s.args[2], x in s.args[1]]
      ψ = streamfunction(z, elements)

      @series begin
          seriestype --> :contour
          grid --> :none
            seriescolor --> cgrad([:grey, :grey])

            s.args[1], s.args[2], ψ
      end
    end
end

@recipe function plot(points::Array{P};
                      source_marker = :xcross,
                      vortex_marker = :circle) where {P <: Union{Points.Point, Blobs.Blob}}
    z = Elements.position(points)
    x = real.(z)
    y = imag.(z)

    sources = Int[]
    vortices = Int[]

    for (i, p) in enumerate(points)
        if abs2(real(p.S)) > 0
            push!(vortices, i)
        end

        if abs2(imag(p.S)) > 0
            push!(sources, i)
        end
    end

    if !isempty(vortices)
        @series begin
            seriestype --> :scatter
            markershape := vortex_marker
            marker_z --> Elements.circulation.(points[vortices])
            x[vortices], y[vortices]
        end
    end

    if !isempty(sources)
        @series begin
            seriestype --> :scatter
            markershape := source_marker
            marker_z --> Elements.flux.(points[sources])
            x[sources], y[sources]
        end
    end
end

@recipe function plot(s::Vortex.Sheet)
    z = s.zs
    Γ = Elements.circulation.(s.blobs)
    line_z --> 0.5(Γ[1:end-1] + Γ[2:end])./abs.(diff(z))
    x := real.(z)
    y := imag.(z)
    ()
end

@recipe function plot(p::Plate)
    z = [p.zs[1], p.zs[end]]
    linecolor --> :black
    x := real.(z)
    y := imag.(z)
    ()
end

@recipe function plot(b::ConformalBody)
    z = [b.zs; b.zs[1]]
    linecolor --> mygreen
    fillrange --> 0
    fillcolor --> mygreen
    aspect_ratio --> 1
    legend --> :none
    x := real.(z)
    y := imag.(z)
    ()
end

@recipe function plot(c::Corner{N};radius=2.0) where {N}

    θ0 = angle(c) - 0.5*N*π*c.ν
    θ1 = angle(c) + 0.5*N*π*c.ν
    θ = range(θ1-2π,θ0,length=200)
    xv = [0;2*radius*cos.(θ);0]
    yv = [0;2*radius*sin.(θ);0]
    z = xv+im*yv

    linecolor --> mygreen
    fillrange --> 0
    fillcolor --> mygreen
    aspect_ratio --> 1
    legend --> :none
    x := real.(z)
    y := imag.(z)
    ()
end

const VortexSystem = NTuple{N, Union{Element, Tuple, Array{V} where {V <: Element}}} where N
function RecipesBase.RecipesBase.apply_recipe(plotattributes::Dict{Symbol, Any}, sys::VortexSystem)
    series_list = RecipesBase.RecipeData[]
    for s in sys
        append!(series_list, RecipesBase.RecipesBase.apply_recipe(copy(plotattributes), s) )
    end
    series_list
end
