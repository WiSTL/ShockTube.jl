module ShockTube

using PyCall 
using Printf
using Unitful
using Dates
using CSV
using DataFrames: DataFrame, mapcols!
using DSP: digitalfilter, Lowpass, Butterworth, filtfilt
using ImageFiltering: mapwindow
using Statistics: median

import Base: convert, ==, isequal, hash, getindex, setindex!, haskey, keys, show

export Species, Mixture
export γ, temperature, pressure, density, molar_density, R_specific, soundspeed
export shockjump!, shockjump, driverpressure!, driverpressure, shockcalc!, shockcalc
export PressureTrace, xt_data

include("init.jl")
# PyCopy = pyimport("copy")
# chem = pyimport_conda("thermo.chemical", "thermo", "conda-forge")
# fluids = pyimport_conda("fluids", "fluids", "conda-forge")

abstract type Fluid end

struct Species <: Fluid
    o::PyObject
end
Species(chemname::String; kwargs...) = Species(PyChem.Chemical(chemname; kwargs...))
Base.show(io::IO, species::Species) = @printf(io, "Species(%s, %0.1f K, %0.3e Pa)", species.ID, species.T, species.P)

struct Mixture <: Fluid
    o::PyObject
end
Mixture(chemnames::Vector{String}; kwargs...) = Mixture(PyChem.Mixture(chemnames; kwargs...))

PyObject(f::Fluid) = getfield(f, :o)
convert(::Type{T}, o::PyObject) where {T <: Fluid} = T(o)
==(f1::Fluid, f2::Fluid) = PyObject(f1) == PyObject(f2)
hash(f::Fluid) = hash(PyObject(f))
Base.copy(f::T) where {T <: Fluid} = T(PyCopy.copy(PyObject(f)))
# Base.Docs.doc(f::Fluid) = Base.Docs.doc(PyObject(f))

Base.getproperty(f::Fluid, s::Symbol) = getproperty(PyObject(f), s)
Base.setproperty!(f::Fluid, s::Symbol, x) = setproperty!(PyObject(f), s, x)
Base.hasproperty(f::Fluid, s::Symbol) = hasproperty(PyObject(f), s)
Base.propertynames(f::Fluid) = propertynames(PyObject(f))
haskey(f::Fluid, x) = haskey(PyObject(f), x)

const MOL_PER_SECOND = typeof(1.0u"mol/s")
mutable struct Flow{F <: Fluid, DU <: Unitful.FreeUnits}
    fluid::F
    rate::MOL_PER_SECOND
    display_unit::DU
end

# Fluid property accessors
temperature(f::Fluid)  = f.T * u"K"
pressure(f::Fluid)  = f.P * u"Pa"
density(f::Fluid)  = f.rho * u"kg/m^3"
molar_density(f::Fluid) = f.rhom * u"mol/m^3"
γ(f::Fluid)  = f.isentropic_exponent
R_specific(f::Fluid)  = f.R_specific * u"J/(kg*K)"
soundspeed(f::Fluid) = sqrt(γ(f) * R_specific(f) * temperature(f)) |> u"m/s"  

# Shock jump conditions
function shockjump!(gas, Mach)
    α1 = (γ(gas) + 1) / (γ(gas) - 1)
    PR = (Mach^2 * (1 + α1) - 1) / α1
    TR = PR * (PR + α1) / (1 + α1 * PR)
    u2 = (soundspeed(gas) * (α1 - 1) * (PR - 1) /
          #-----------------------------
           √((1 + α1) * (1 + α1 * PR)))
    gas.P = gas.P * PR
    gas.T = gas.T * TR
    return gas, uconvert(u"m/s", u2)
end

shockjump(gas, Mach) = shockjump!(copy(gas), Mach)

function driverpressure!(driver, driven, Ms)
    γ1, γ4 = γ(driven), γ(driver)
    a1, a4 = soundspeed(driven), soundspeed(driver)
    driver.P = driven.P * (1 + 2γ1 / (γ1 + 1) * (Ms^2 - 1)) *
                    (1 + a1 / a4 * (γ4 - 1) / (γ1 + 1) * (1/Ms - Ms)) ^ (2γ4 / (1 - γ4))
    return driver
end
driverpressure(driver, driven, Ms) = driverpressure!(copy(driver), driven, Ms)

function shockcalc!(driver, driven, Ms)
    shocked, u2 = shockjump(driven, Ms)
    driverpressure!(driver, driven, Ms)
    return (driver = driver, driven = driven, shocked=shocked, u2=u2)
end
shockcalc(driver, driven, Ms) = shockcalc!(copy(driver), driven, Ms)

struct PressureTrace
    t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    timestamp::DateTime
    data::DataFrame
end

function PressureTrace(filepath, filter=(median_size = 15, lowpass_mul=100, filt=Butterworth(2)))
    open(filepath) do f
        header = Dict(Pair(split(readline(f), '\t', limit=3)[1:2]...) for i = 1:23)

        t = range(parse(Float64, header["X0"]), 
        length = parse(Int, header["Samples"]), 
          step = parse(Float64, header["Delta_X"]))

        date = Date(header["Date"], "yyyy/mm/dd")
        time = Time(first(split(header["Time"], '.')))
        timestamp = DateTime(date, time)

        seekstart(f)
        data = CSV.File(f, datarow=24, drop=[1], delim='\t', 
            header="PT" .* string.(0:12)) |> DataFrame
        
        
        if !isnothing(filter)
            freq = 1/step(t)
            lowpass = digitalfilter(Lowpass(freq/filter.lowpass_mul, fs=freq), filter.filt)
            mapcols!(data) do s
                s̄ = mapwindow(median, s, filter.median_size)
                ŝ = filtfilt(lowpass, s)
                [median((s[i], s̄[i], ŝ[i])) for i in eachindex(s)]
            end
        end
        PressureTrace(t, timestamp, data)
    end
end

function xt_data(ptrace, driver, driven, ptloc_filepath, trigger::Pair, ref_PT::Symbol)
    ptlocs = CSV.File(ptloc_filepath) |> DataFrame
    trigger_PT, trigger_thresh = trigger

    # Get PT indices (PT1 => 1, PT9 => 9, etc) within ptlocs
    trigger_idx = parse(Int, string(trigger_PT)[3:end])
    ref_idx = parse(Int, string(ref_PT)[3:end])

    # find trigger time index
    trigger_sensor = >(ustrip(trigger_thresh |> u"psi"))
    trigger_ptrace_idx = findfirst(trigger_sensor, ptrace.data[!, trigger_PT])

    ref_ptrace_idx = findfirst(trigger_sensor, ptrace.data[!, ref_PT])
    Δt = ptrace.t[ref_ptrace_idx] - ptrace.t[trigger_ptrace_idx]
    Δx = ptlocs[ref_idx, :x_m] - ptlocs[trigger_idx, :x_m]

    # calculate Mach number between trigger and reference PTs
    W₀ = (Δx / Δt)*u"m/s"
    M₀ = W₀/soundspeed(driven)

    # determine gas states
    states = shockcalc(driver, driven, M₀)
    shocked_psi = pressure(states.shocked) |> u"psi" |> ustrip
    t_shock = [isnothing(i) ? NaN : ptrace.t[i] for i in findfirst.(>(0.5*shocked_psi), eachcol(ptrace.data))]

    return states, DataFrame(:x => ptlocs[!, :x_m], :t_shock => t_shock)
end

end # module