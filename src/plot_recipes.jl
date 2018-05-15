using RecipesBase
using ColorTypes
import PlotUtils: cgrad

@userplot Streamlines

@recipe function f(s::Streamlines)
    elements = s.args[3]

    if ConformalBody in typeof.(elements)

        b = elements[findfirst(isa.(elements,PotentialFlow.ConformalBody))]
        ζ = [r*exp(im*θ) for θ in s.args[2], r in s.args[1]]
        ψ = streamfunction(ζ,elements)
        Z = conftransform(ζ,b)
        vmin = minimum(ψ)
        vmax = maximum(ψ)
        vmean = vmax-vmin
        vfact = 1
        n = 31
        vmin = vmean+vfact*(vmin-vmean)
        vmax = vmean+vfact*(vmax-vmean)
        v = linspace(vmin,vmax,n)
        if (sign(vmin) != sign(vmax))
            v = v .- v[abs.(v).==minimum(abs.(v))];
        end
        @series begin
          seriestype --> :contour
            grid --> :none
            ratio --> 1
            linewidth --> 1
            seriescolor --> [:black, :black]
            levels --> v
            real.(Z), imag.(Z), ψ
        end

        z = [b.zs; b.zs[1]]
        @series begin
          linecolor --> RGBA{Float64}(151/255,180/255,118/255,1.0)
          fillrange --> 0
          fillcolor --> RGBA{Float64}(151/255,180/255,118/255,1.0)
          ratio --> 1
          legend --> :none
          real.(z), imag.(z)
        end

    else
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
                      vortex_marker = :circle) where {P <: Union{Points.Point{T} where T, Blobs.Blob{T} where T}}
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
    linecolor --> :black
    fillrange --> 0
    fillcolor --> :black
    ratio --> 1
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
