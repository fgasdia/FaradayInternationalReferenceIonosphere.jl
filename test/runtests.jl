using FIRITools, Test, Statistics
using CSV, DataFrames


@testset "FIRITools.jl" begin
    @test size(FIRITools.DATA) == (94, 1980)  # for firi2018
    @test size(FIRITools.HEADER) == (1980, 5)

    @test FIRITools.twoclosest(FIRITools.HEADER[!,"Chi, deg"], 90) == [90, 90]
    @test FIRITools.twoclosest(FIRITools.HEADER[!,"Chi, deg"], 89) == [85, 90]
    @test FIRITools.twoclosest(FIRITools.HEADER[!,"Chi, deg"], 135) == [100, 130]

    @test size(FIRITools.selectprofiles()) == (94, 1980)
    @test size(FIRITools.selectprofiles(month=(1, 6))) == (94, 1980÷2)
    @test_logs (:warn, "15 is not in model parameters. Looking for interpolating form of `selectprofiles`?") FIRITools.selectprofiles(chi=15) 
    
    # Try interpolating chi and lat. Should always have same size
    @test size(FIRITools.selectprofiles(15, 10)) == (94, 36)  # interpolate both
    @test size(FIRITools.selectprofiles(15.1, 82.3)) == (94, 36)  # extrapolate
    @test size(FIRITools.selectprofiles(30, 30)) == (94, 36)  # no interpolation
    @test size(FIRITools.selectprofiles(30, 10)) == (94, 36) # interpolate lat
    @test size(FIRITools.selectprofiles(29, 30)) == (94, 36) # interpolate chi

    profs85 = FIRITools.selectprofiles(85, 31)
    profs89 = FIRITools.selectprofiles(89, 31)
    profs90 = FIRITools.selectprofiles(90, 31)
    @test all(profs85 .> profs89 .> profs90)
    # This test doesn't work for lat because the profiles are largely overlapping.
    # Even using `firi` instead of `selectprofiles` doesn't work because the mean of the
    # profiles causes there to be some crossing density regions

    @test size(firi()) == (94,)
    @test size(firi(month=(1, 6))) == (94,)
    @test size(firi(15, 10)) == (94,)

    @test_throws ArgumentError firi(-10, 45)
    @test_throws ArgumentError firi(30, -15)
    @test_logs (:warn, "`chi` greater than 130° uses `chi = 130°`") firi(160, 45)

    @test FIRITools.quantile(1, chi=(0, 90)) ==
        Statistics.quantile.(eachrow(FIRITools.selectprofiles(chi=(0, 90))), 1)

    q1 = FIRITools.quantile(chi=(0, 90), (0.3, 0.5))
    @test q1[:,1] == FIRITools.quantile(chi=(0, 90), 0.3)
    @test q1[:,2] == FIRITools.quantile(chi=(0, 90), 0.5)

    function quantilewoextrap(ref, p)
        # Compare to reference profiles without extrapolation
        refmask = [a in FIRITools.ALTITUDE[:,1] for a in ref.alt]
        fmask = [a in ref.alt for a in FIRITools.ALTITUDE[:,1]] 
        isapprox(ref.ne[refmask],
                 FIRITools.quantile(chi=(0, 95), p)[fmask], rtol=0.01)
    end

    @test quantilewoextrap(CSV.read("FIRIday_30.csv", DataFrame), 0.3)
    @test quantilewoextrap(CSV.read("FIRIday_90.csv", DataFrame), 0.9)

    # plot(ref30.ne, ref30.alt, label="Wei 30", xlims=(10^4, 4e10), ylims=(55, 110), xscale=:log10)
    # plot!(FIRITools.quantile(chi=(0, 95), 0.3), FIRITools.ALTITUDE, label=30, xscale=:log10)

    @test_throws ArgumentError FIRITools.extrapolate(firi(), 30e3:1e3:120e3)

    firialt = minimum(FIRITools.ALTITUDE):1000:maximum(FIRITools.ALTITUDE)
    p1 = FIRITools.extrapolate(firialt, firi(), 30e3:1e3:120e3; max_altitude=60e3)
    p2 = FIRITools.extrapolate(firialt, firi(), 30e3:1e3:120e3; N=6)
    p3 = FIRITools.extrapolate(firi(), 30e3:1e3:120e3; max_altitude=60e3)
    p4 = FIRITools.extrapolate(firi(), 30e3:1e3:120e3; N=6)

    p5 = FIRITools.extrapolate(firi(), 70e3:1e3:120e3; max_altitude=60e3)

    @test !any(isinf, p1)
    @test !any(isnan, p1)

    @test length(p1) == length(30:120)
    @test length(p5) == length(70:120)

    @test p1 == p2
    @test p1 == p3
    @test p1 == p4
    @test p1[60 .< 30:120] == firi()[60e3 .< FIRITools.ALTITUDE .<= 120e3]
    @test p1[70 .<= 30:120] == p5

    # plot(firi(), FIRITools.ALTITUDE, xscale=:log10)
    # plot!(FIRITools.extrapolate(firi(), 30:148, N=6), 30:148, xscale=:log10)
    # plot!(p5, 70:120, xscale=:log10)

    # ep = FIRITools.extrapolate(FIRITools.quantile(chi=(0, 95), 0.5), 40:110; max_altitude=60)
    # plot(ep, 40:110, xscale=:log10)

    p11 = FIRITools.extrapolate(FIRITools.DATA[:,[1, 2]], 30e3:1e3:120e3, N=6)
    @test p11[:,1] == FIRITools.extrapolate(FIRITools.DATA[:,1], 30e3:1e3:120e3, N=6)
    @test p11[:,2] == FIRITools.extrapolate(FIRITools.DATA[:,2], 30e3:1e3:120e3; N=6)

    newz = 30e3:250:110e3
    @test length(FIRITools.extrapolate(firi(), newz, N=6)) == length(newz)
end
