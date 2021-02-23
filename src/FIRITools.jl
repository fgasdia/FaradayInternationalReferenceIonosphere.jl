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

export firi, quantile

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

Return the average FIRI profile across model parameters selected by the keyword arguments.
"""
function firi(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    profiles = selectprofiles(chi=chi, lat=lat, f10_7=f10_7, doy=doy, month=month)
    return dropdims(mean(profiles, dims=2), dims=2)
end

"""
    firi(chi, lat; f10_7=(75, 200), doy=(15, 350), month=(1, 12))

Return the average FIRI profile across the model parameters selected by the keyword
arguments and interpolated at solar zenith angle `chi` and latitude `lat`.
"""
function firi(chi, lat; f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    profiles = selectprofiles(chi, lat, f10_7=f10_7, doy=doy, month=month)
    return dropdims(mean(profiles, dims=2), dims=2)
end

"""
    selectprofiles(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))

Return matrix of selected model profiles.
"""
function selectprofiles(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    mask = trues(nrow(HEADER))

    mask .&= select(HEADER[!,"Chi, deg"], chi) .&
             select(HEADER[!,"Lat, deg"], lat) .&
             select(HEADER[!,:F10_7], f10_7)   .&
             select(HEADER[!,:DOY], doy)       .&
             select(HEADER[!,:Month], month)

    return DATA[:,mask]
end

"""
Return matrix of selected model profiles interpolated at solar zenith angle `chi` and
latitude `lat`.
"""
function selectprofiles(chi, lat; f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    mask = trues(nrow(HEADER))

    mask .&= select(HEADER[!,:F10_7], f10_7) .&
             select(HEADER[!,:DOY], doy)     .&
             select(HEADER[!,:Month], month)

    function twoclosest(s, v)
        uvals = unique(HEADER[!,s])
        I = sortperm(abs.(v .- uvals))
        return [uvals[I[1]], uvals[I[2]]]  # 2 closest matches
    end

    chibnds = twoclosest("Chi, deg", chi)
    latbnds = twoclosest("Lat, deg", lat)

    # Mask all but the 2 closest values of chi and lat to the target chi and lat
    mask .&= .!isnothing.(indexin(HEADER[!,"Chi, deg"], chibnds)) .&
             .!isnothing.(indexin(HEADER[!,"Lat, deg"], latbnds))

    N = count(mask)
    maskidxs = findall(mask)

    profiles = Matrix{Float64}(undef, length(ALTITUDE), N÷4)

    # Although not a general solution, we know the format of Chi guarantees that occurences
    # of the two closest values will be in pairs and that both Chi and Lat are sorted
    # so every 4 indices in the mask form a group that covers
    # [(chi_low, lat_low) (chi_high, lat_low);
    #  (chi_low, lat_high) (chi_high, lat_high)]
    j = 1
    for i = 1:4:N
        chis = HEADER[[maskidxs[i], maskidxs[i+1]], "Chi, deg"]
        lats = HEADER[[maskidxs[i], maskidxs[i+2]], "Lat, deg"]
        for a in eachindex(ALTITUDE)
            d = [DATA[a,maskidxs[i]] DATA[a,maskidxs[i+2]];
                 DATA[a,maskidxs[i+1]] DATA[a,maskidxs[i+3]]]
            itp = LinearInterpolation((chis, lats), d, extrapolation_bc=Line())
            profiles[a,j] = itp(chi, lat)
        end
        j += 1
    end

    return profiles
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
    x in h || @warn "$x is not in model parameters. Looking for interpolating form of `selectprofiles`?"
    
    return x .== h
end

"""
    quantile(p; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))

Compute the quantile(s) `p` at each `ALTITUDE` with data columns selected by the keyword
arguments.
"""
function quantile(p; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    profiles = selectprofiles(chi=chi, lat=lat, f10_7=f10_7, doy=doy, month=month)
    profile = _quantile(profiles, p)

    if p isa Number
        # make `profile` a Vector for consistency with `firi`
        profile = dropdims(profile, dims=2)
    end

    return profile
end

"""
    quantile(p; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), doy=(15, 350), month=(1, 12))

Compute the quantile(s) `p` at each `ALTITUDE` with data columns selected by keyword
arguments and interpolated to solar zenith angle `chi` and latitude `lat`.
"""
function quantile(chi, lat, p; f10_7=(75, 200), doy=(15, 350), month=(1, 12))
    profiles = selectprofiles(chi, lat, f10_7=f10_7, doy=doy, month=month)
    profile = _quantile(profiles, p)

    if p isa Number
        # make `profile` a Vector for consistency with `firi`
        profile = dropdims(profile, dims=2)
    end

    return profile
end

function _quantile(profiles, p)
    profile = Array{Float64, 2}(undef, length(ALTITUDE), length(p))
    for j in eachindex(p)
        i = 1
        for r in eachrow(profiles)
            profile[i,j] = Statistics.quantile(r, p[j])
            i += 1
        end
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