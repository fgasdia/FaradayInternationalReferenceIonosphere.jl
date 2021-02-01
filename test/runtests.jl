using FIRITools, Test, Statistics
using CSV, DataFrames

@testset "FIRITools.jl" begin
    @test size(profile()) == (94, 1980)  # for firi2018
    @test size(header()) == (1980, 7)

    @test all(buildmask())
    @test size(buildmask()) == (1980,)
    @test size(buildmask(doy=(1, 180))) == (1980,)

    @test profile()[!, buildmask(chi=(0, 90))] == profile(buildmask(chi=(0, 90)))

    @test FIRITools.quantile(buildmask(chi=(0, 90)), 1) ==
        Statistics.quantile.(eachrow(profile(buildmask(chi=(0, 90)))), 1)

    q1 = FIRITools.quantile(buildmask(chi=(0, 90)), (0.3, 0.5))
    @test q1[:,1] == FIRITools.quantile(buildmask(chi=(0, 90)), 0.3)
    @test q1[:,2] == FIRITools.quantile(buildmask(chi=(0, 90)), 0.5)

    function quantilewoextrap(ref, p)
        # Compare to reference profiles without extrapolation
        refmask = [a in FIRITools.ALTITUDE[:,1] for a in ref.alt]
        fmask = [a in ref.alt for a in FIRITools.ALTITUDE[:,1]] 
        isapprox(ref.ne[refmask],
                FIRITools.quantile(buildmask(chi=(0, 95)), p)[fmask], rtol=0.01)
    end

    @test quantilewoextrap(CSV.read("FIRIday_30.csv", DataFrame), 0.3)
    @test quantilewoextrap(CSV.read("FIRIday_90.csv", DataFrame), 0.9)

    # plot(ref30.ne, ref30.alt, label="Wei 30", xlims=(10^4, 4e10), ylims=(55, 110), xscale=:log10)
    # plot!(FIRITools.quantile(buildmask(chi=(0, 95)), 0.3), FIRITools.ALTITUDE[:,1], label=30, xscale=:log10)

    p1 = FIRITools.extrapolate(FIRITools.altitude(), profile()[!,"1"], 30:120; max_altitude=60)
    p2 = FIRITools.extrapolate(FIRITools.altitude(), profile()[!,"1"], 30:120; N=6)
    p3 = FIRITools.extrapolate(profile()[!,"1"], 30:120; max_altitude=60)
    p4 = FIRITools.extrapolate(profile()[!,"1"], 30:120; N=6)

    p5 = FIRITools.extrapolate(profile()[!,"1"], 70:120; max_altitude=60)

    @test length(p1) == length(30:120)
    @test length(p5) == length(70:120)

    @test p1 == p2
    @test p1 == p3
    @test p1 == p4
    @test p1[60 .< 30:120] == profile()[60 .< FIRITools.altitude() .<= 120,"1"]
    @test p1[70 .<= 30:120] == p5

    # plot(profile()[!,"1"], FIRITools.altitude(), xscale=:log10)
    # plot!(FIRITools.extrapolate(profile()[!,"1"], 30:148, N=6), 30:148, xscale=:log10)
    # plot!(p5, 70:120, xscale=:log10)

    # ep = FIRITools.extrapolate(FIRITools.quantile(buildmask(chi=(0, 95)), 0.5),
    #                            40:110; max_altitude=60)
    # plot(ep, 40:110, xscale=:log10)

    p11 = FIRITools.extrapolate(Matrix(profile()[!,["1", "2"]]), 30:120, N=6)
    @test p11[:,1] == p4
    @test p11[:,2] == FIRITools.extrapolate(profile()[!,"2"], 30:120; N=6)
end