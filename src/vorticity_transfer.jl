export impulse_matching_correction, transfer_circulation!, sheettip_to_blob!, newest_element,
        shed_new_vorticity!


"""
    shed_new_vorticity!(sys,motion,δ,t,lesp=0.0,tesp=0.0)

Returns the plate/vortex system after adding new vortices of the appropriate type
at the plate edges, in order to satisfy the respective regularization conditions at those edges.
"""
function shed_new_vorticity!(sys, motion, δ, t, lesp = 0.0, tesp = 0.0)
    plate, ambient_sys = sys
    le_sys, te_sys = ambient_sys


    nlev = new_vortex(le_sys..., plate.zs[end], plate.α, δ, plate, motion)
    ntev = new_vortex(te_sys..., plate.zs[1], plate.α, δ, plate, motion)

    Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, nlev, ntev, t, lesp, tesp)

    plate, (add_new_vortex!(le_sys..., nlev, Γ₊), add_new_vortex!(te_sys..., ntev, Γ₋))

end

"""
    impulse_matching_correction(vs, vt, plate) -> ComplexF64

Compute the impulse matching position correction for a target vortex `vt` when the circulation
of a source vortex `vs` is to be fully transferred to it. Equation (A26) of Darakananda and Eldredge JFM 2019
with Γ̇ = ΓₛΔt. Assumes a flat plate body.
"""
function impulse_matching_correction(vs, vt, plate::Plate)
    @unpack c, α, L = plate

    c⁺ = plate.zs[end]
    c⁻ = plate.zs[1]

    z̃t = 2(vt.z - c)*exp(-im*α)/L
    z̃s = 2(vs.z - c)*exp(-im*α)/L

    pt = Plates.unit_impulse(z̃t)
    ps = Plates.unit_impulse(z̃s)

    η = z̃t/(√(z̃t + 1)*√(z̃t - 1))

    Δp = ps - pt
    Δz = im*0.5L*exp(im*α)*(Δp*(1 + η') + Δp'*(η' - 1))/(η + η')*(circulation(vs)/circulation(vt))
#     if isnan(Δz) || isinf(Δz)
#         error("Δz is $Δz")
#     end
    return Δz
end


"""
    transfer_circulation!(sheet,point,Δz,plate,ϵ=1e-3;[max_segments=10]) -> Vortex.Blob

Return a vortex blob that is initially identical to `point` and, where appropriate,
accrues the circulation of the last several segments of the vortex sheet `sheet`.
In these cases, it also removes the last several segments of `sheet`. The position
correction `Δz`, which accounts for the impulse matching correction due to the
transfer of circulation from each segment, is added to the position of `point`
to create the output vortex blob. `ϵ` is the maximum allowable discrepancy in
the impulse due to the circulation transfer (compared to the impulse if no circulation
is transferred). `max_segments` is the maximum number of sheet segments that can be considered
for transfer; it defaults to 10.
"""
function transfer_circulation!(sheet, point, Δz, plate, ϵ = 1e-3; max_segments=10)
    z₀  = Elements.position(point)
    Γ₀  = circulation(point)

    ΣΔz = zero(ComplexF64)
    ΣΓ = 0.0

    Σp = Plates.vortex_impulse(plate,Vortex.Point(z₀, Γ₀)) #total_impulse(Vortex.Point(z₀, Γ₀), plate)

    segments = 0

    # If the sheet isn't long enough to transfer, then return with no transfer
    if length(sheet) < max_segments+2
        return point
    end

    for n in 1:length(sheet)-max_segments  # Don't eat into the last `max_segments` segments
        if isnan(Δz[n]) || isinf(Δz[n])
            if n ≤ 2
                return point
            else
                segments = n-1
                break;
            end
        end

        ΣΔz += Δz[n]
        ΣΓ  += circulation(sheet.blobs[n])

        p₊ = Plates.vortex_impulse(plate,Vortex.Point(z₀ + ΣΔz, Γ₀ + ΣΓ)) # impulse of active vortex with this segment transferred
        Σp += Plates.vortex_impulse(plate,sheet.blobs[n])  # impulse of the sheet if this segment is not transferred


        if abs2(p₊ - Σp) > ϵ^2
            if n ≤ 2
                return point
            else
                segments = n-1
                break
            end
        end
        segments += 1
    end

    # Change the position and strength of `point`, using the impulse correction
    # and strength of each segment
    z = z₀ + sum(Δz[1:segments-1]) + 0.5Δz[segments]
    Γ = Γ₀ + sum(circulation(sheet.blobs[1:segments-1])) +
                 0.5(sheet.Ss[segments] -sheet.Ss[segments-1])
    point₊ = Vortex.Blob(z, Γ, point.δ)

    # Remove the transferred segments from the sheet
    splice!(sheet.blobs, 1:segments-1)
    splice!(sheet.Ss, 1:segments-1)
    sheet.blobs[1] = sheet.blobs[1](Γ = 0.5(sheet.Ss[2] - sheet.Ss[1]))

    return point₊
end

"""
    sheettip_to_blob!(blobs,sheet)

Transform the tip of a vortex sheet into a vortex blob. Add it to the beginning of the list
of blobs (and remove it from the sheet).
"""
function sheettip_to_blob!(blobs::Vector{T},sheet::Vortex.Sheet) where T <: Vortex.Blob
    tips = sheet.blobs[1:2]

    # factor of 2 in circulation is due to the trapezoidal weight on blob strength of first blob
    new_blob = Vortex.Blob(0.5*(tips[1].z + tips[2].z), 2circulation(tips[1]), tips[1].δ)

    # remove the first entry from the sheet, since this has been transferred to the active blob
    splice!(sheet.blobs, 1)
    splice!(sheet.Ss, 1)

    # now the second sheet element is the first element and needs to have its strength adjusted for trapezoidal weight
    sheet.blobs[1] = Vortex.Blob(tips[2].z, 0.5*(sheet.Ss[2] - sheet.Ss[1]), tips[2].δ)
    pushfirst!(blobs, new_blob)
    nothing
end


"""
    new_vortex(edge,angle,δ,plate,motion) -> Vortex.Blob

Create a new vortex blob near the edge `edge` of the plate, with unit strength
"""
new_vortex(edge, α, δ, plate, motion) = Vortex.Blob(edge - 1e-2im*exp(im*α)*sign(Plates.normal(motion.ċ+im*motion.α̇*(edge-plate.c),α)), 1.0, δ)


"""
    new_vortex(blobs::Vector{Vortex.Blob},edge,angle,δ,plate,motion) -> Vortex.Sheet

Create a new vortex sheet near the edge `edge` of the plate with unit strength, using the vortex blobs in `blobs`
to place the sheet segments. (The blobs are placed at the sheet segment centroids.)
"""
function new_vortex(blobs::Vector{T}, edge, α, δ, plate, motion) where T <: Vortex.Blob
    Δz = (blobs[end].z - edge)/2
    #Vortex.Sheet([0.5, 1.5] .* Δz .+ edge, [0.0, 1.0], δ)
    Vortex.Sheet([1.5, 0.5] .* Δz .+ edge, [1.0, 0.0], δ)
end

"""
    new_vortex(blobs::Vector{Vortex.Blob},sheet::Vortex.Sheet,edge,angle,δ,plate,motion) -> Vortex.Blob

Create vortex blobs near the edge `edge` of the plate with unit strength, using the last vortex sheet segment in `sheet`
to place the blobs.
"""
function new_vortex(::Vector{T}, sheet::Vortex.Sheet, edge, α, δ, plate, motion) where T <: Vortex.Blob
    z = (edge + 2sheet.zs[end])/3
    Vortex.Blob.([sheet.zs[end], z], [0.5, 0.5], δ)
end




"""
    add_new_vortex!(blob::Vortex.Blob,Γ) -> (blob,)

Create a new vortex blob with the position and radius of `blob`, but with strength `Γ`,
unless `Γ` is zero, in which case it returns an empty tuple.
"""
function add_new_vortex!(b::Vortex.Blob, Γ)
    if Γ != 0
        ([Vortex.Blob(b.z, Γ, b.δ)],)
    else
        ()
    end
end

"""
    add_new_vortex!(blobs::Vector{Vortex.Blob},sheet,Γ) -> (blob,)

Return blobs and a sheet with the same positions and radii as `sheet`, but with the strengths
set to 0 and `Γ`. This is meant for cases in which the sheet is initially created.
"""
function add_new_vortex!(blobs::Vector{T}, sheet::Vortex.Sheet, Γ) where T <: Vortex.Blob
    blobs, Vortex.Sheet(sheet.zs, [0.0, Γ], sheet.blobs[1].δ)
end

"""
    add_new_vortex!(blobs::Vector{Vortex.Blob},sheet,segment::Vector{Vortex.Blob},Γ) -> (blob,)

Return blobs and a sheet that appends the last vortex in `segment` to the end of `sheet`, with strength `Γ`.
This is the newly-created entry from the edge.
"""
function add_new_vortex!(blobs::Vector{T}, sheet::Vortex.Sheet, segment::Vector{T}, Γ) where T <: Vortex.Blob
    Vortex.Sheets.append_segment!(sheet, segment[2].z, Γ)
    blobs, sheet
end


"""
    newest_element() -> ()

Returns an empty tuple
"""
newest_element() = ()

"""
    newest_element(blobs::Vector{Vortex.Blob}) -> (Vortex.Blob,)

Returns the most recently added vortex in `blobs`
"""
newest_element(blobs::Vector{T}) where T <: Vortex.Blob = (blobs[1],)

"""
    newest_element(blobs::Vector{Vortex.Blob},sheet::Vortex.Sheet) -> (Vortex.Blob,)

Returns the most recently added vortex in `sheet`
"""
function newest_element(blobs::Vector{T}, sheet::Vortex.Sheet) where T <: Vortex.Blob
    (sheet.blobs[end],)
end
