"""
    FIRITools

`FIRITools` is a collection of convenience functions for working with
Faraday-International Reference Ionosphere (FIRI) model profiles.

The underlying FIRI-2018 model `DATA`, `HEADER`, and `ALTITUDE`'s can be accessed
as e.g. `FIRITools.HEADER`. The unit of `ALTITUDE` is meters and `DATA` is electrons per
cubic meter.

Only the function `firi` is exported, but also see `FIRITools.quantile` and
`FIRITools.extrapolate`.

# References

Friedrich, M., Pock, C., & Torkar, K. (2018).
FIRI-2018, an updated empirical model of the lower ionosphere.
Journal of Geophysical Research: Space Physics, 123, 6737-6751.
[doi: 10.1029/2018JA025437](https://doi.org/10.1029/2018JA025437)
"""
module FIRITools

using Artifacts, Statistics
using CSV, DataFrames
using Interpolations

export firi

const MIN_DENSITY = 1e-4

const DF = CSV.read(joinpath(artifact"firi", "firi2018.csv"), DataFrame)

function parseheader()
    types = Dict(:Code=>String, :Month=>Int, "Chi, deg"=>Int, "Lat, deg"=>Int, :F10_7=>Int)

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

# Start altitude at 60 km even though table has down to 55 km (60 km is stated lower limit)
const DATA = convert(Matrix, parse.(Float64, DF[13:end-2, 2:end]))  # last 2 rows have Missing
const ALTITUDE = parse.(Int, DF[13:end-2, 1])*1000
@assert length(ALTITUDE) == size(DATA, 1) "ALTITUDE does not match DATA"

"""
    values(s)

Return the FIRI model values of parameter `s`.
"""
function values(s)
    if s in ("f10_7", :f10_7)
        return [75, 130, 200]  # == unique(FIRITools.HEADER[!,:F10_7])
    elseif s in ("chi", :chi)
        return [0, 30, 45, 60, 75, 80, 85, 90, 95, 100, 130]
    elseif s in ("lat", :lat)
        return [0, 15, 30, 45, 60]
    elseif s in ("month", :month)
        return [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    else
        throw(ArgumentError("$s is not an FIRI model parameter"))
    end
end


"""
    firi(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))

Return the average FIRI profile across model parameters selected by the keyword arguments.
"""
function firi(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))
    profiles = selectprofiles(chi=chi, lat=lat, f10_7=f10_7, month=month)
    profile = mean(profiles, dims=2)
    replace!(x -> x<MIN_DENSITY ? MIN_DENSITY : x, profile)
    return dropdims(profile, dims=2)
end

"""
    firi(chi, lat; f10_7=(75, 200), month=(1, 12))

Return the average FIRI profile across the model parameters selected by the keyword
arguments and interpolated at solar zenith angle `chi` and latitude `lat`.
"""
function firi(chi, lat; f10_7=(75, 200), month=(1, 12))
    profiles = selectprofiles(chi, lat, f10_7=f10_7, month=month)
    profile = mean(profiles, dims=2)
    replace!(x -> x<MIN_DENSITY ? MIN_DENSITY : x, profile)
    return dropdims(profile, dims=2)
end

"""
    selectprofiles(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))

Return matrix of selected model profiles.
"""
function selectprofiles(;chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))
    mask = trues(nrow(HEADER))

    mask .&= select(HEADER[!,"Chi, deg"], chi) .&
             select(HEADER[!,"Lat, deg"], lat) .&
             select(HEADER[!,:F10_7], f10_7)   .&
             select(HEADER[!,:Month], month)

    return DATA[:,mask]
end

"""
    selectprofiles(chi, lat; f10_7=(75, 200), month=(1, 12))

Return matrix of selected model profiles interpolated at solar zenith angle `chi` and
latitude `lat`.
"""
function selectprofiles(chi, lat; f10_7=(75, 200), month=(1, 12))
    issorted(unique(HEADER[!, "Chi, deg"])) || throw(DomainError("Chi column of HEADER is not sorted"))
    issorted(unique(HEADER[!, "Lat, deg"])) || throw(DomainError("Lat column of HEADER is not sorted"))

    if chi > 130
        @warn "`chi` greater than 130° uses `chi = 130°`" maxlog=3
        chi = oftype(chi, 130)
    end

    if chi < 0 || lat < 0
        throw(ArgumentError("`chi` or `lat` below 0° is not supported. See Friedrich et al., “FIRI-2018”."))
    end

    mask = trues(nrow(HEADER))

    mask .&= select(HEADER[!,:F10_7], f10_7) .&
             select(HEADER[!,:Month], month)

    chibnds = twoclosest(HEADER[!,"Chi, deg"], chi)
    latbnds = twoclosest(HEADER[!,"Lat, deg"], lat)

    # Mask all but the 2 closest values of chi and lat to the target chi and lat
    mask .&= .!isnothing.(indexin(HEADER[!,"Chi, deg"], chibnds)) .&
             .!isnothing.(indexin(HEADER[!,"Lat, deg"], latbnds))

    N = count(mask)
    maskidxs = findall(mask)

    # We need to explicitly handle 1D vs 2D interpolation (if lat, or chi don't need interp)
    if isequal(chibnds...) && isequal(latbnds...)
        profiles = DATA[:, maskidxs]
    elseif isequal(chibnds...)
        profiles = Matrix{Float64}(undef, length(ALTITUDE), N÷2)

        j = 1
        for i = 1:2:N
            for a in axes(DATA, 1)  # altitude
                d = DATA[a,maskidxs[[i, i+1]]]
                itp = LinearInterpolation(latbnds, d, extrapolation_bc=Line())
                profiles[a,j] = itp(lat)
            end
            j += 1
        end
    elseif isequal(latbnds...)
        profiles = Matrix{Float64}(undef, length(ALTITUDE), N÷2)

        j = 1
        for i = 1:2:N
            for a in axes(DATA, 1)  # altitude
                d = DATA[a,maskidxs[[i, i+1]]]
                itp = LinearInterpolation(chibnds, d, extrapolation_bc=Line())
                profiles[a,j] = itp(chi)
            end
            j += 1
        end
    else
        # Although not a general solution, we know the format of Chi guarantees that occurences
        # of the two closest values will be in pairs and that both Chi and Lat are sorted
        # so every 4 indices in the mask form a group that covers
        # [(chi_low, lat_low) (chi_high, lat_low);
        #  (chi_low, lat_high) (chi_high, lat_high)]
        profiles = Matrix{Float64}(undef, length(ALTITUDE), N÷4)

        j = 1
        for i = 1:4:N
            for a in axes(DATA, 1)  # altitude
                d = [DATA[a,maskidxs[i]] DATA[a,maskidxs[i+2]];
                     DATA[a,maskidxs[i+1]] DATA[a,maskidxs[i+3]]]
                itp = LinearInterpolation((chibnds, latbnds), d, extrapolation_bc=Line())
                profiles[a,j] = itp(chi, lat)
            end
            j += 1
        end
    end

    return profiles
end

"""
    twoclosest(s, v)

Return:
    - the two closest values of vector `s` to `v` if `v` is outside of the range of
        values in `s` (extrapolating)
    - or the closest value on each side of `v` in column `s` if interpolating
    - or `[v, v]` if `v` is exactly in `s`.

The two values are sorted.
"""
function twoclosest(s, v)
    uvals = unique(s)

    if v in uvals
        return [v, v]
    elseif v >= maximum(uvals) || v <= minimum(uvals)
        # We're extrapolating or at a min/max
        I = sortperm(abs.(v .- uvals))
        return sort([uvals[I[1]], uvals[I[2]]])  # 2 closest matches
    else
        # We're interpolating, make sure we return the closest on each side of `v`
        udiffs = v .- uvals
        v1 = minimum(x->x > 0 ? x : Inf, udiffs)  # closest val with positive diff
        v2 = -minimum(x->x < 0 ? abs(x) : Inf, udiffs)  # closest val with negative diff
        
        I1 = findfirst(isequal(v1), udiffs)
        I2 = findfirst(isequal(v2), udiffs)
        return sort([uvals[I1], uvals[I2]])
    end
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
    quantile(p; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))

Compute the quantile(s) `p` at each `ALTITUDE` with data columns selected by the keyword
arguments.
"""
function quantile(p; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))
    profiles = selectprofiles(chi=chi, lat=lat, f10_7=f10_7, month=month)
    profile = _quantile(profiles, p)

    replace!(x -> x<MIN_DENSITY ? MIN_DENSITY : x, profile)

    if p isa Number
        # make `profile` a Vector for consistency with `firi`
        profile = dropdims(profile, dims=2)
    end

    return profile
end

"""
    quantile(chi, lat, p; f10_7=(75, 200), month=(1, 12))

Compute the quantile(s) `p` at each `ALTITUDE` with data columns selected by keyword
arguments and interpolated to solar zenith angle `chi` and latitude `lat`.
"""
function quantile(chi, lat, p; f10_7=(75, 200), month=(1, 12))
    profiles = selectprofiles(chi, lat, f10_7=f10_7, month=month)
    profile = _quantile(profiles, p)

    replace!(x -> x<MIN_DENSITY ? MIN_DENSITY : x, profile)

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


"""
    extrapolate([z,] profile, newz)

Return `profile`, originally defined at sorted altitudes `z`, exponentially interpolated
and extrapolated to altitudes of `newz`.

If `z` is not specified, it is assumed that `profile` is defined at `FIRITools.ALTITUDE`.
"""
function extrapolate(z::AbstractRange, profile::AbstractVector, newz::AbstractRange)
    length(z) == length(profile) ||
        throw(ArgumentError("`z` and `profile` must be equal lengths."))
    
    itp = LinearInterpolation(z, log.(profile), extrapolation_bc=Line())
    itpprofile = itp(newz)

    return exp.(itpprofile)
end
extrapolate(profile, newz) = extrapolate(ALTITUDE[begin]:1000:ALTITUDE[end], profile, newz)    

function extrapolate(z::AbstractRange, profile::AbstractMatrix, newz::AbstractRange)
    expprofile = Vector{eltype(profile)}()
    nrows = 0
    for i in axes(profile, 2)
        p = extrapolate(z, view(profile,:,i), newz)
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