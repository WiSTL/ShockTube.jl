using ShockTube
using Test

@testset "ShockTube.jl" begin
    st = shockcalc(Species("N2"), Species("Ar"), 1.8)
    @test isapprox(ustrip(pressure(st.driver)), 1.426894182009389e6, rtol=2e-3)
    @test isapprox(ustrip(st.u2), 300.0825207636464, rtol=2e-3)
    @test isapprox(ustrip(density(st.shocked)), 3.3867027987446874, rtol=2e-3)
    @test isapprox(ustrip(soundspeed(st.reflected)), 544.5317089952606, rtol=5e-3)
end
