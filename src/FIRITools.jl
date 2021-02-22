"""
    FIRITools

`FIRITools` is a collection of convenience functions for working with
Faraday-International Reference Ionosphere (FIRI) model profiles.

The underlying FIRI-2018 model `DATA`, `HEADER`, and `ALTITUDE`'s can be accessed
as e.g. `FIRITools.HEADER`.

# References

Friedrich, M., Pock, C., & Torkar, K. (2018).
FIRI-2018, an updated empirical model of the lower ionosphere.
Journal of Geophysical Research: Space Physics, 123, 6737-6751.
[doi: 10.1029/2018JA025437](https://doi.org/10.1029/2018JA025437)
"""
module FIRITools

using Artifacts, Statistics
using CSV, DataFrames
using LsqFit, Interpolations

export firi, buildmask

const DF = CSV.read(joinpath(artifact"firi", "firi2018.csv"), DataFrame)

function parseheader()
    types = Dict(:Code=>String, :Month=>Int, :DOY=>Int, "Chi, deg"=>Int, "Lat, deg"=>Int,
                 :F10_7=>Int)

    tmpheader = permutedims(first(DF, 6), 1)

    df = DataFrame()
    for (s, t) in types
        if t == String
            df[!,s] = tmpheader[!,s]
        else
            df[!,s] = parse.(t, tmpheader[!,s])
        end
    end
    return df
end

const HEADER = parseheader()
const DATA = convert(Matrix, parse.(Float64, DF[8:end-2, 2:end]))  # last 2 rows have Missing
const ALTITUDE = parse.(Int, DF[8:end-2, 1])
@assert length(ALTITUDE) == size(DATA, 1) "ALTITUDE does not match DATA"


"""
    firi(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))

Return a copy of profile with columns masked by `mask` averaged at each altitude.
"""
function firi(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    mask = trues(nrow(HEADER))

    cols = Dict("Chi, deg"=>chi, "Lat, deg"=>lat, "F10_7"=>f10_7, "DOY"=>doy, "Month"=>month)
    for (s,v) in cols
        m = select(HEADER[!,s], v)
        if !isnothing(m)
            mask .&= m
            delete!(cols, s)
        else
            @info "Interpolating $s"
        end
    end

    # The fields left in `cols` need to be interpolated
    newcols = Matrix{Float64}(undef, length(ALTITUDE), 0)
    newmask = trues(nrow(HEADER))
    for (s,v) in cols
        uvals = unique(HEADER[!,s])
        I = sortperm(abs.(v .- uvals))
        bounds = (uvals[I[1]], uvals[I[2]])

        # The lower and upper bound appear many times for every unmasked combination of the
        # other model parameters
        mask_lower = (HEADER[!,s] .== bounds[1]) .& mask
        mask_upper = (HEADER[!,s] .== bounds[2]) .& mask
        N = count(mask_lower)
        @assert N == count(mask_upper) "Unequal number of instances of bound values in header"

        mask_lower_idxs = findall(mask_lower)
        mask_upper_idxs = findall(mask_upper)

        newcols = hcat(newcols, Matrix{Float64}(undef, length(ALTITUDE), N))

        # Generate a single interpolated column of data for every combination of the other
        # acceptable parameters
        for i = 1:N
            d = [DATA[:,mask_lower_idxs[i]] DATA[:,mask_upper_idxs[i]]]
            for a in eachindex(ALTITUDE)
                itp = LinearInterpolation(bounds, d[a,:], extrapolation_bc=Line())
                newcols[a,i] = itp(v)
            end

            newmask[mask_lower_idxs[i]] = 0
            newmask[mask_upper_idxs[i]] = 0
        end
    end

    mask .&= newmask  # mask out the columns that were used for interpolation
    profile = DATA[:,mask]
    profile = hcat(profile, newcols)  # append the columns resulting from interpolation

    return dropdims(mean(profile, dims=2), dims=2)
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

Return `x .== h` if `h` is in `x`, otherwise `nothing`.
"""
function select(h, x::T) where T<:Number
    if h in x
        return x .== h
    else
        return nothing
    end
end

"""
    quantile(mask, p)

Compute the quantile(s) `p` at each `ALTITUDE` with data columns masked by `mask`.
"""
function quantile(mask, p::T) where T<:Number
    profile = Vector{Float64}(undef, length(ALTITUDE))
    i = 1
    for r in eachrow(DATA[!,mask])
        profile[i] = Statistics.quantile(r, p)
        i += 1
    end
    return profile
end

function quantile(mask, p)
    profile = Array{Float64, 2}(undef, length(ALTITUDE), length(p))
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