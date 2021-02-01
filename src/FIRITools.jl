"""
    FIRITools

`FIRITools` is a collection of convenience functions for working with
Faraday-International Reference Ionosphere (FIRI) model profiles.

# References

Friedrich, M., Pock, C., & Torkar, K. (2018).
FIRI-2018, an updated empirical model of the lower ionosphere.
Journal of Geophysical Research: Space Physics, 123, 6737-6751.
[doi: 10.1029/2018JA025437](https://doi.org/10.1029/2018JA025437)
"""
module FIRITools

using Artifacts, Statistics
using CSV, DataFrames
using LsqFit

export header, profile, buildmask

# TODO: Just read in the file once and then split up components appropriately
HEADER = CSV.read(joinpath(artifact"firi", "firi2018.csv"), DataFrame,
    transpose=true, select=1:7)
DATA = CSV.read(joinpath(artifact"firi", "firi2018.csv"), DataFrame,
    skipto=9, footerskip=2, select=2:1981)
ALTITUDE = CSV.read(joinpath(artifact"firi", "firi2018.csv"), DataFrame,
    skipto=9, footerskip=2, select=[1], header=["altitude, km"], threaded=false, silencewarnings=true)
@assert nrow(ALTITUDE) == nrow(DATA) "ALTITUDE does not match DATA"


"""
    profile()

Return a copy of profile data as a `DataFrame`.
"""
function profile()
    return copy(DATA)
end

"""
    profile(mask)

Return a copy of profile data with columns masked by `mask`.

See also: [`buildmask`](@ref)
"""
function profile(mask)
    return DATA[:,mask]
end

"""
    header()

Return a copy of the data header as a `DataFrame`.
"""
function header()
    return copy(HEADER)
end

"""
    header(mask)

Return a copy of the data header as a `DataFrame` masked by `mask`.
"""
function header(mask)
    return HEADER[mask,:]
end

"""
    altitude()

Return FIRI2018 profile altitudes.
"""
function altitude()
    return ALTITUDE[:,1]
end

"""
    select(h, x::Tuple)

Return `x[1] .<= h .<= x[2]`.
"""
function select(h, x::T) where T<:Tuple
    first(x) .<= h .<= last(x)
end

"""
    select(h, x)

Return `x .== h`.
"""
function select(h, x::T) where T<:Number
    x .== h
end
# TODO: Array-like or range type

"""
    buildmask(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))

Construct a `BitArray` mask to be applied across the columns of `FIRITools.DATA`.

With `Tuple` arguments, the first and last elements of the tuple are treated as the extents
of the mask field, inclusive. If arguments are of `Number` type, an exact match is required.

By default, no masking occurs.

# Examples

To look at daytime ionospheres only:
```julia-repl
julia> mask = buildmask(chi=(0, 90))
```
"""
function buildmask(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    mask = trues(nrow(HEADER))

    mask .&= select(HEADER[!,"Chi, deg"], chi) .&
             select(HEADER[!,"Lat, deg"], lat) .&
             select(HEADER[!,:F10_7], f10_7)   .&
             select(HEADER[!,:DOY], doy)       .&
             select(HEADER[!,:Month], month)

    return mask
end

"""
    quantile(mask, p)

Compute the quantile(s) `p` at each `ALTITUDE` with data columns masked by `mask`.
"""
function quantile(mask, p::T) where T<:Number
    profile = Vector{Float64}(undef, nrow(ALTITUDE))
    i = 1
    for r in eachrow(DATA[!,mask])
        profile[i] = Statistics.quantile(r, p)
        i += 1
    end
    return profile
end

function quantile(mask, p)
    profile = Array{Float64, 2}(undef, nrow(ALTITUDE), length(p))
    for j in eachindex(p)
        profile[:,j] = quantile(mask, p[j])
    end
    return profile
end

expmodel(z, p) = p[1]*exp.(z*p[2])

function jacobian_expmodel(z, p)
    J = Matrix{eltype(p)}(undef, length(z), length(p))
    J[:,1] = exp.(p[2] .* z)       # df/dp[1]
    J[:,2] = z .* p[1] .* J[:,1]   # df/dp[2]
    return J
end

"""
    extrapolate([z,] profile, newz; kwargs...)

Return `profile`, originally defined at sorted altitudes `z`, exponentially extrapolated
at altitudes of `newz` below `z`.

If `z` is not specified, it is assumed that `profile` is defined at `FIRITools.altitude()`.

The extrapolation is performed using profile samples from the bottom altitude of `z` up to
and including `max_altitude`.

If `N` is provided, then `max_altitude` is automatically determined as the `N`th element of
`z`.
"""
function extrapolate(z, profile::AbstractVector, newz; max_altitude=nothing, N=nothing)
    isnothing(max_altitude) && isnothing(N) &&
        throw(ArgumentError("At least one of `max_altitude` or `N` are required."))
    issorted(z) || throw(ArgumentError("`z` must be sorted."))
    length(z) == length(profile) ||
        throw(ArgumentError("`z` and `profile` must be equal lengths."))
    maximum(newz) <= maximum(z) ||
        @warn "`newz` extends above `maximum(z)`. `profile` will only be extrapolated below `max_altitude`."
    
    if !isnothing(N)
        max_altitude = z[N]
    end

    p0 = [1.0, 0.3]
    mask = z .<= max_altitude
    fit = curve_fit(expmodel, jacobian_expmodel, z[mask], profile[mask], p0)

    extrapz = filter(x -> x <= max_altitude, newz)
    reducedmask = max_altitude > minimum(newz) ? (max_altitude .< z .<= maximum(newz)) :
        (minimum(newz) .<= z .<= maximum(newz))

    return [expmodel(extrapz, fit.param); profile[reducedmask]]
end
extrapolate(profile, newz; max_altitude=nothing, N=nothing) =
    extrapolate(altitude(), profile, newz; max_altitude=max_altitude, N=N)    

function extrapolate(z, profile::AbstractMatrix, newz; max_altitude=nothing, N=nothing)
    expprofile = Vector{eltype(profile)}()
    nrows = 0
    for i in axes(profile, 2)
        p = extrapolate(z, view(profile,:,i), newz, max_altitude=max_altitude, N=N)
        append!(expprofile, p)
        nrows = length(p)
    end
    return reshape(expprofile, nrows, size(profile, 2))
end


# NOTE: at this time, models() isn't very important because only 1 model is supported
"""
    models()

Return the available models.
"""
models() = first.(splitext.(filter(endswith(".csv"), readdir(artifact"firi"))))

end  # module