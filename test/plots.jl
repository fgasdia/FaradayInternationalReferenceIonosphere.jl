using FIRITools
using Plots

const alt = FIRITools.ALTITUDE

# Compare mean and median profiles
day_mean = firi(chi=(0, 95))
day_med = FIRITools.quantile(0.5, chi=(0, 95))
day_30 = FIRITools.quantile(0.3, chi=(0, 95))
day_70 = FIRITools.quantile(0.7, chi=(0, 95))

plot(day_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Day", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(day_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(day_30, alt/1000, label="30%", xscale=:log10)
plot!(day_70, alt/1000, label="70%", xscale=:log10)

night_mean = firi(chi=(100, 180))
night_med = FIRITools.quantile(0.5, chi=(100, 180))
night_30 = FIRITools.quantile(0.3, chi=(100, 180))
night_70 = FIRITools.quantile(0.7, chi=(100, 180))

plot(night_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Night", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(night_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(night_30, alt/1000, label="30%", xscale=:log10)
plot!(night_70, alt/1000, label="70%", xscale=:log10)

# Compare mean and median around interpolation points
day_mean = firi(chi=(0, 95), lat=(30, 45))
day_med = FIRITools.quantile(0.5, chi=(0, 95), lat=(30, 45))
day_30 = FIRITools.quantile(0.3, chi=(0, 95), lat=(30, 45))
day_70 = FIRITools.quantile(0.7, chi=(0, 95), lat=(30, 45))

plot(day_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Day, limited latitude", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(day_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(day_30, alt/1000, label="30%", xscale=:log10)
plot!(day_70, alt/1000, label="70%", xscale=:log10)

day_mean = firi(chi=(30, 45))
day_med = FIRITools.quantile(0.5, chi=(30, 45))
day_30 = FIRITools.quantile(0.3, chi=(30, 45))
day_70 = FIRITools.quantile(0.7, chi=(30, 45))

plot(day_mean, alt/1000, label="mean", xscale=:log10, linewidth=2, ylims=(55, 110),
     legend=:topleft, title="Day, limited SZA", ylabel="altitude (km)", xlabel="Ne (m⁻³)")
plot!(day_med, alt/1000, label="median", xscale=:log10, linewidth=2)
plot!(day_30, alt/1000, label="30%", xscale=:log10)
plot!(day_70, alt/1000, label="70%", xscale=:log10)

# Check to see if interpolation is acting as expected

real_lats = unique(FIRITools.HEADER[!,"Lat, deg"])

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10, legend=:topleft,
         title="Latitude", ylims=(55, 110))
for la in real_lats
    prof = firi(45, la)
    plot!(p, prof, alt/1000, label=la)
end
display(p)


real_sza = unique(FIRITools.HEADER[!,"Chi, deg"])
interp_sza = setdiff(0:5:180, real_sza)
N = length(0:5:180)
cmap = palette(:rainbow, N)

p = plot(ylabel="altitude (km)", xlabel="Ne (m⁻³)", xscale=:log10, legend=false,
         title="SZA", ylims=(55, 110))
for sza in real_sza
    prof = firi(sza, 45)
    plot!(p, prof, alt/1000, label=sza, color=cmap[findfirst(isequal(sza), 0:5:180)])
end

for sza in interp_sza
    prof = firi(sza, 45)
    plot!(p, prof, alt/1000, label=sza, color=cmap[findfirst(isequal(sza), 0:5:180)],
          linestyle=:dash)
end

display(p)
