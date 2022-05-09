using FaradayInternationalReferenceIonosphere
import FaradayInternationalReferenceIonosphere as FIRI
using Plots

const alt = FIRI.ALTITUDE

# Compare mean and median profiles
day_mean = firi(chi=(0, 95))
day_med = FIRI.quantile(0.5, chi=(0, 95))
day_30 = FIRI.quantile(0.3, chi=(0, 95))
day_70 = FIRI.quantile(0.7, chi=(0, 95))

plot(day_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Day", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(day_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(day_30, alt/1000, label="30%", xscale=:log10)
plot!(day_70, alt/1000, label="70%", xscale=:log10)
savefig("day_mean.png")

night_mean = firi(chi=(100, 180))
night_med = FIRI.quantile(0.5, chi=(100, 180))
night_30 = FIRI.quantile(0.3, chi=(100, 180))
night_70 = FIRI.quantile(0.7, chi=(100, 180))

plot(night_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Night", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(night_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(night_30, alt/1000, label="30%", xscale=:log10)
plot!(night_70, alt/1000, label="70%", xscale=:log10)
savefig("night_mean.png")

# Compare mean and median around interpolation points
day_mean = firi(chi=(0, 95), lat=(30, 45))
day_med = FIRI.quantile(0.5, chi=(0, 95), lat=(30, 45))
day_30 = FIRI.quantile(0.3, chi=(0, 95), lat=(30, 45))
day_70 = FIRI.quantile(0.7, chi=(0, 95), lat=(30, 45))

plot(day_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Day, limited latitude", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(day_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(day_30, alt/1000, label="30%", xscale=:log10)
plot!(day_70, alt/1000, label="70%", xscale=:log10)
savefig("day_mean_latzoom.png")

day_mean = firi(chi=(30, 45))
day_med = FIRI.quantile(0.5, chi=(30, 45))
day_30 = FIRI.quantile(0.3, chi=(30, 45))
day_70 = FIRI.quantile(0.7, chi=(30, 45))

plot(day_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Day, limited SZA", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(day_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(day_30, alt/1000, label="30%", xscale=:log10)
plot!(day_70, alt/1000, label="70%", xscale=:log10)
savefig("day_mean_chizoom.png")

# Check to see if interpolation is acting as expected

real_lats = unique(FIRI.HEADER.lat)
lats = 0:5:90
interp_lat = setdiff(lats, real_lats)
N = length(lats)
cmap = palette(:rainbow, N, rev=true)

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10, legend=:topleft,
         ylims=(55, 110), legendtitle="Latitude")
for la in real_lats
    prof = firi(45, la)
    plot!(p, prof, alt/1000, label=la, color=cmap[findfirst(isequal(la), lats)])
end
for la in interp_lat
    prof = firi(45, la)
    plot!(p, prof, alt/1000, label=nothing, color=cmap[findfirst(isequal(la), lats)],
          linestyle=:dash)
end
display(p)
savefig("lat.png")

real_sza = unique(FIRI.HEADER.chi)
chis = 0:5:130
interp_sza = setdiff(chis, real_sza)
# interp_sza = [105, 110, 115, 120]
N = length(chis)
cmap = palette(:rainbow, N, rev=true)

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10,
         ylims=(55, 110), legend=:topleft, legendtitle="SZA")
for sza in real_sza
    prof = firi(sza, 45)
    plot!(p, prof, alt/1000, label=sza, color=cmap[findfirst(isequal(sza), chis)])
end

for sza in interp_sza
    prof = firi(sza, 45)
    plot!(p, prof, alt/1000, label=nothing, color=cmap[findfirst(isequal(sza), chis)],
          linestyle=:dash)
end
display(p)
savefig(p, "sza.png")


real_sza = [75, 80, 85, 90, 95, 100]
chis = 75:100
interp_sza = setdiff(chis, real_sza)
N = length(chis)
cmap = palette(:rainbow, N, rev=true)

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10,
         ylims=(55, 110), legend=:topleft, legendtitle="SZA")
for sza in real_sza
    prof = firi(sza, 45)
    plot!(p, prof, alt/1000, label=sza, color=cmap[findfirst(isequal(sza), chis)])
end

for sza in interp_sza
    prof = firi(sza, 45)
    plot!(p, prof, alt/1000, label=nothing, color=cmap[findfirst(isequal(sza), chis)],
          linestyle=:dash)
end
display(p)
savefig("sza_zoom.png")


###
# With extrapolation to ground

real_lats = unique(FIRI.HEADER.lat)
lats = 0:5:90
interp_lat = setdiff(lats, real_lats)
N = length(lats)
cmap = palette(:rainbow, N, rev=true)
newalt = 0:110

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10, legend=:topleft,
         ylims=(0, 110), yticks=0:20:110, legendtitle="Latitude", xticks=exp10.([0, 3, 6, 9, 12]))
for la in real_lats
    prof = FIRI.extrapolate(firi(45, la), newalt*1e3)
    plot!(p, prof, newalt, label=la, color=cmap[findfirst(isequal(la), lats)])
end
for la in interp_lat
    prof = FIRI.extrapolate(firi(45, la), newalt*1e3)
    plot!(p, prof, newalt, label=nothing, color=cmap[findfirst(isequal(la), lats)],
          linestyle=:dash)
end
display(p)
savefig("lat_extrap.png")

real_sza = unique(FIRI.HEADER.chi)
chis = 0:1:130
interp_sza = setdiff(chis, real_sza)
# interp_sza = [105, 110, 115, 120]
N = length(chis)
cmap = palette(:rainbow, N, rev=true)

p = plot(ylabel="Altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10, xlims=(1e3, 5e11),
         ylims=(40, 110), yticks=0:10:110, legend=:topleft, legendtitle="SZA", xticks=exp10.([0, 4, 6, 8, 10]))
for sza in real_sza
    prof = FIRI.extrapolate(firi(sza, 45), newalt*1e3)
    plot!(p, prof, newalt, label=sza, color=cmap[findfirst(isequal(sza), chis)], linewidth=2)
end

for sza in interp_sza
    prof = FIRI.extrapolate(firi(sza, 45), newalt*1e3)
    plot!(p, prof, newalt, label=nothing, color=cmap[findfirst(isequal(sza), chis)],
          linestyle=:dash)
end
annotate!([(5e6, 48,
    ("FIRI profiles at 45°N,\ninterpolated over solar zenith angle,\naveraged over all months and F10.7 indices,\nand extrapolated below 60km.", 9, :black, :left))])
display(p)
savefig(p, "sza_extrap.png")

# Using median
real_sza = unique(FIRI.HEADER.chi)
chis = 0:5:130
interp_sza = setdiff(chis, real_sza)
# interp_sza = [105, 110, 115, 120]
N = length(chis)
cmap = palette(:rainbow, N, rev=true)

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10,
         ylims=(0, 110), yticks=0:20:110, legend=:topleft, legendtitle="SZA", xticks=exp10.([0, 3, 6, 9, 12]))
for sza in real_sza
    prof = FIRI.extrapolate(FIRI.quantile(sza, 45, 0.5), newalt*1e3)
    plot!(p, prof, newalt, label=sza, color=cmap[findfirst(isequal(sza), chis)])
end

for sza in interp_sza
    prof = FIRI.extrapolate(FIRI.quantile(sza, 45, 0.5), newalt*1e3)
    plot!(p, prof, newalt, label=nothing, color=cmap[findfirst(isequal(sza), chis)],
          linestyle=:dash)
end
display(p)
savefig(p, "sza_median_extrap.png")
