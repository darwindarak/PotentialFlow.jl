module Source

import ..Points
import ..Blobs
import ..Elements: circulation

#== Wrapper for a point source ==#

const Point = Points.Point{Complex128}
Point(s::Point; z = s.z, S = imag(s.S)) = Point(z, S)
function Base.show(io::IO, s::Point)
    if iszero(real(s.S))
        print(io, "Source.Point($(s.z), $(imag(s.S))")
    else
        print(io, "Points.Point($(s.z), $(s.S)")
    end
end
circulation(::Point) = 0.0

#== Wrapper for a blob source ==#

const Blob = Blobs.Blob{Complex128}
Blob(s::Blob; z = s.z, S = imag(s.S), δ = s.δ) = Blob(z, S, δ)
function Base.show(io::IO, s::Blob)
    if iszero(real(s.S))
        print(io, "Source.Blob($(s.z), $(imag(s.S)), $(s.δ)")
    else
        print(io, "Blobs.Blob($(s.z), $(imag(s.S)), $(s.δ)")
    end
end
circulation(::Blob) = 0.0

end
