using RecipesBase

@userplot Streamlines

@recipe function f(s::Streamlines)
    elements = s.args[3]

    z = [x + im*y for y in s.args[2], x in s.args[1]]
    ψ = streamfunction(z, elements)

    @series begin
        seriestype --> :contour
        grid --> :none

        s.args[1], s.args[2], ψ
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

const VortexSystem = NTuple{N, Union{Element, Tuple, Array{V} where {V <: Element}}} where N
function RecipesBase.RecipesBase.apply_recipe(plotattributes::Dict{Symbol, Any}, sys::VortexSystem)
    series_list = RecipesBase.RecipeData[]
    for s in sys
        append!(series_list, RecipesBase.RecipesBase.apply_recipe(copy(plotattributes), s) )
    end
    series_list
end
