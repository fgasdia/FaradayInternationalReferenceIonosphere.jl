using FaradayInternationalReferenceIonosphere
import FaradayInternationalReferenceIonosphere as FIRI
using Test, Statistics, DelimitedFiles
using TypedTables

@testset "FaradayInternationalReferenceIonosphere.jl" begin
    ALTITUDE_LENGTH = 89
    @test size(FIRI.DATA) == (ALTITUDE_LENGTH, 1980)  # for firi2018
    @test length(FIRI.HEADER) == 1980
    @test length(columnnames(FIRI.HEADER)) == 5

    @test FIRI.values(:f10_7) == FIRI.values("f10_7")
    @test FIRI.values(:f10_7) == unique(FIRI.HEADER.f10_7)
    @test FIRI.values(:chi) == unique(FIRI.HEADER.chi)
    @test FIRI.values(:lat) == unique(FIRI.HEADER.lat)
    @test FIRI.values(:month) == unique(FIRI.HEADER.month)

    @test FIRI.twoclosest(FIRI.HEADER.chi, 90) == [90, 90]
    @test FIRI.twoclosest(FIRI.HEADER.chi, 89) == [85, 90]
    @test FIRI.twoclosest(FIRI.HEADER.chi, 135) == [100, 130]

    @test size(FIRI.selectprofiles()) == (ALTITUDE_LENGTH, 1980)
    @test size(FIRI.selectprofiles(month=(1, 6))) == (ALTITUDE_LENGTH, 1980÷2)
    @test_logs (:warn, "15 is not in model parameters. Looking for interpolating form of `selectprofiles`?") FIRI.selectprofiles(chi=15) 
    
    # Try interpolating chi and lat. Should always have same size
    @test size(FIRI.selectprofiles(15, 10)) == (ALTITUDE_LENGTH, 36)  # interpolate both
    @test size(FIRI.selectprofiles(15.1, 82.3)) == (ALTITUDE_LENGTH, 36)  # extrapolate
    @test size(FIRI.selectprofiles(30, 30)) == (ALTITUDE_LENGTH, 36)  # no interpolation
    @test size(FIRI.selectprofiles(30, 10)) == (ALTITUDE_LENGTH, 36) # interpolate lat
    @test size(FIRI.selectprofiles(29, 30)) == (ALTITUDE_LENGTH, 36) # interpolate chi

    profs85 = FIRI.selectprofiles(85, 31)
    profs89 = FIRI.selectprofiles(89, 31)
    profs90 = FIRI.selectprofiles(90, 31)
    @test all(profs85 .> profs89 .> profs90)
    # This test doesn't work for lat because the profiles are largely overlapping.
    # Even using `firi` instead of `selectprofiles` doesn't work because the mean of the
    # profiles causes there to be some crossing density regions

    @test size(firi()) == (ALTITUDE_LENGTH,)
    @test size(firi(month=(1, 6))) == (ALTITUDE_LENGTH,)
    @test size(firi(15, 10)) == (ALTITUDE_LENGTH,)

    @test firi(160, 45) == firi(130, 45)  # floor chi to 130°

    @test_throws ArgumentError firi(-10, 45)
    @test_throws ArgumentError firi(30, -15)

    @test FIRI.quantile(1, chi=(0, 90)) ==
        Statistics.quantile.(eachrow(FIRI.selectprofiles(chi=(0, 90))), 1)

    q1 = FIRI.quantile(chi=(0, 90), (0.3, 0.5))
    @test q1[:,1] == FIRI.quantile(chi=(0, 90), 0.3)
    @test q1[:,2] == FIRI.quantile(chi=(0, 90), 0.5)

    function quantilewoextrap(ref, p)
        # Compare to reference profiles without extrapolation
        v, h = ref  # values, header
        refmask = [a in FIRI.ALTITUDE[:,1] for a in v[:,1]]
        fmask = [a in v[:,1] for a in FIRI.ALTITUDE[:,1]] 
        isapprox(v[refmask,2],
                 FIRI.quantile(chi=(0, 95), p)[fmask], rtol=0.01)
    end

    @test quantilewoextrap(readdlm("FIRIday_30.csv", ','; header=true), 0.3)
    @test quantilewoextrap(readdlm("FIRIday_90.csv", ','; header=true), 0.9)

    # plot(ref30.ne, ref30.alt, label="Wei 30", xlims=(10^4, 4e10), ylims=(55, 110), xscale=:log10)
    # plot!(FIRI.quantile(chi=(0, 95), 0.3), FIRI.ALTITUDE, label=30, xscale=:log10)

    firialt = minimum(FIRI.ALTITUDE):1000:maximum(FIRI.ALTITUDE)
    p1 = FIRI.extrapolate(firialt, firi(), 30e3:1e3:120e3)
    p3 = FIRI.extrapolate(firi(), 30e3:1e3:120e3)
    p5 = FIRI.extrapolate(firi(), 70e3:1e3:120e3)

    @test !any(isinf, p1)
    @test !any(isnan, p1)

    @test length(p1) == length(30:120)
    @test length(p5) == length(70:120)

    @test p1 == p3
    @test p1[65 .< 30:120] ≈ firi()[65e3 .< FIRI.ALTITUDE .<= 120e3]
    @test p1[70 .<= 30:120] == p5

    # plot(firi(), FIRI.ALTITUDE/1000, xscale=:log10, legend=false)
    # plot!(FIRI.extrapolate(firi(), (30:148)*1e3; N=6), 30:148)
    # plot!(p5, 70:120)

    # ep = FIRI.extrapolate(FIRI.quantile(chi=(0, 95), 0.5), 40:110; max_altitude=60)
    # plot(ep, 40:110, xscale=:log10)

    p11 = FIRI.extrapolate(FIRI.DATA[:,[1, 2]], 30e3:1e3:120e3)
    @test p11[:,1] == FIRI.extrapolate(FIRI.DATA[:,1], 30e3:1e3:120e3)
    @test p11[:,2] == FIRI.extrapolate(FIRI.DATA[:,2], 30e3:1e3:120e3)

    newz = 30e3:250:110e3
    @test length(FIRI.extrapolate(firi(), newz)) == length(newz)
end
