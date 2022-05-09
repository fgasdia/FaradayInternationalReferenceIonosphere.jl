# FaradayInternationalReferenceIonosphere.jl

[![DOI](https://zenodo.org/badge/332354802.svg)](https://zenodo.org/badge/latestdoi/332354802)

**Julia tools for working with FIRI ionosphere profiles.**

The Faraday-International Reference Ionosphere (FIRI) is a semiempirical model of the nonauroral ionosphere that enhances and extends the IRI model down to 60 km and densities above 10⁶ e-/m³. FIRI blends results from a simple ion-chemical model with sounding rocket profiles made by Langmuir probe and HF Faraday rotation measurements. FIRI-2018 is the most recent version of the model by M. Friedrich, C. Pock, and K. Torkar, but it was originally published in 2001 by Friedrich and Torkar. See the [Citations](#citations) section below.

## Usage

Only one function is exported from FaradayInternationalReferenceIonosphere: `firi`. `firi` returns the average profile after filtering the model by the keyword arguments to the function.

```julia
profile = firi(; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))
```

Each of the keyword arguments can either be a `Tuple` that inclusively brackets the range of values to be included in the average profile or can be a single value. Note that no interpolation occurs! If the value specified for `f10_7` is `100` (which is not one of the values in the model output) then a warning will be printed and an empty profile will be 
returned.

> :star: **Tip:**
>
> For convenience, load the package into your environment like
> ```julia
> using FaradayInternationalReferenceIonosphere
> import FaradayInternationalReferenceIonosphere as FIRI
> ```
> The first line brings `firi` into scope. If that's all you need, then you can stop there. However, the second line lets you reference the package using the shorthand `FIRI`. I assume the package has been loaded this way in the examples below.

The raw model data can be accessed as `FIRI.ALTITUDE` for a `Vector` of the altitude in meters, `FIRI.DATA` for a `Matrix` of each electron density profile in electrons per cubic meter, and `FIRI.HEADER` for a `Table` with each combination of FIRI model parameters. To find the acceptable values of the `f10_7` field, for example, you can use

```julia
unique(FIRI.HEADER.f10_7)
```

Although day of year ("DOY") appears as an independent field in the original published FIRI model file, this is entirely redundant with the month number. Each day of year simply corresponds to the day number in the middle of each month. Therefore, we use only the month field.

The keyword defaults bracket the full range of each field so that no filtering occurs.

### Interpolation across latitude and solar zenith angle

The latitude and solar zenith angle (`chi`) fields can be interpolated with the function form

```julia
profile = firi(chi, lat; f10_7=(75, 200), month=(1, 12))
```

This performs a simple linear interpolation at each altitude of the profile for the specified `chi` and `lat`. The latitude can also be linearly extrapolated beyond 60°. If a `chi` or `lat` less than 0° is provided, a warning is thrown. According to Friedrich et al., 2018, although FIRI is only valid for the Northern Hemisphere, it can be used to predict values for the Southern Hemisphere by adding 6 months. Any `chi` greater than 130° (night) returns the profiles for exactly 130°.

### Quantiles

The average profile returned by `firi` is skewed towards higher electron densities. If preferred, profile quantiles can be returned with the function `FIRI.quantile`. The quantile function is not exported to avoid collisions with `Statistics.quantile`. The usage of `FIRI.quantile` is similar to `firi`, except the quantile(s) are specified by `p` on the interval [0, 1].

```julia
profile = FIRI.quantile(p; chi=(0, 130), lat=(0, 60), f10_7=(75, 200), month=(1, 12))
```

There is also an interpolating form

```julia
profile = FIRI.quantile(chi, lat, p; f10_7=(75, 200), month=(1, 12))
```

### Extrapolation by altitude

To extrapolate the bottom of a profile to lower altitudes, use `FIRI.extrapolate`. This function will also interpolate to new (finer or coarser) altitudes in `newz`.

The function call is

```julia
newprofile = FIRI.extrapolate([z,] profile, newz)
```

## Citations

Friedrich, M., and Torkar, K. M. (2001). FIRI: A semiempirical model of the lower ionosphere. Journal of Geophysical Research: Space Physics, 106, 21409-21418. [doi: 10.1029/2001JA900070](https://doi.org/10.1029/2001JA900070)

Friedrich, M., Pock, C., & Torkar, K. (2018). FIRI-2018, an updated empirical model of the lower ionosphere. Journal of Geophysical Research: Space Physics, 123, 6737-6751. [doi: 10.1029/2018JA025437](https://doi.org/10.1029/2018JA025437)
