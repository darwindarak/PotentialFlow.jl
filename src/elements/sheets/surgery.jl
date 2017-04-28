using Dierckx

function truncate!(sheet, n::Int)
    @assert 2 ≤ n < length(sheet)-1
    ΔΓ = sheet.Γs[n] - sheet.Γs[1]
    splice!(sheet.Γs, 1:n-1)
    splice!(sheet.blobs, 1:n-1)
    sheet.blobs[1] = Vortex.Blob(sheet.blobs[1].z, 0.5(sheet.Γs[2] - sheet.Γs[1]), sheet.δ)
    return ΔΓ
end

function remesh!(sheet, zs, Γs)
    sheet.Γs = Γs
    sheet.blobs = Vortex.Blob.(zs, compute_trapezoidal_weights(Γs), sheet.δ)
    sheet
end

function append_segment!(sheet::Sheet, z, Γ)
    b₋ = sheet.blobs[end]
    sheet.blobs[end] = Vortex.Blob(b₋.z, b₋.Γ + 0.5Γ, b₋.δ)
    push!(sheet.blobs, Vortex.Blob(z, 0.5Γ, sheet.δ))
    push!(sheet.Γs, sheet.Γs[end] + Γ)
    nothing
end

function filter_by_arclength(z, Δs, Δf, properties...)
    @assert Δs < Δf

    l = vcat(0, accumulate(+, abs.(diff(z))))

    n = round(Int, l[end]/Δs)
    if n ≤ 1
        warn("Sampling interval too large for arc length")
        return z, properties
    end
    
    x = Spline1D(l, real.(z); k=1, bc="nearest", s=0.0)
    y = Spline1D(l, imag.(z); k=1, bc="nearest", s=0.0)

    splines = map(properties) do p
        Spline1D(l, p; k=1, bc="nearest", s=0.0)
    end

    L = linspace(0, l[end], n)
    ΔL = step(L)

    z₌ = @. evaluate(x, L) + im*evaluate(y, L)
    splines₌ = map(splines) do spline
        evaluate.(spline, L)
    end

    ẑ = dct(z₌)
    Ŝ = dct.(splines₌)

    cutoff = ceil(Int, 2L[end]/Δf) + 1
    ẑ[cutoff:end] = zero(Complex128)
    foreach(Ŝ) do S
        S[cutoff:end] = 0.0
    end

    idct(ẑ), idct.(Ŝ)
end
