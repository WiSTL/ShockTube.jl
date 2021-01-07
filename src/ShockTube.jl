module ShockTube

using PyCall 
using Printf
using Reexport
@reexport using Unitful
using Dates
using CSV
using DataFrames: DataFrame, mapcols!
using DSP: digitalfilter, Lowpass, Butterworth, filtfilt
using ImageFiltering: mapwindow
using Statistics: median

# using Reexport
# @reexport using Unitful: @u_str

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

"""
    xt_data(ptrace, ptlocs, driver, driven, trigger, ref_PT)

Calculate the time of shock arrival at each pressure transducer.
Uses the driven gas properties and the Δx/Δt between the trigger
and reference pressure transducers to calculate an approximate Mach
number, which is used to determine driver conditions and shocked 
gas pressure for determination of shock arrival time at each transducer.

    Arguments:
    - `ptrace::PressureTrace`
    - `ptlocs::DataFrame`
    - `driver::Fluid`, e.g. `Species("Helium")`
    - `driven::Fluid`, e.g. `Mixture(["Helium", "Acetone"], zs=[0.95, 0.05])`
    - `trigger::Pair{Symbol, Unitful.Quantity}`, e.g. `:PT3 => 10u"psi"`
    - `ref_PT::Symbol`, e.g. `:PT6`

# Examples
```julia-repl
julia> ptlocs = CSV.File("PT_locations.csv") |> DataFrame
12×3 DataFrame
 Row │ name    x_m       σ_m
     │ String  Float64   Float64
─────┼──────────────────────────────
   1 │ PT1     0.506412  0.00254
   2 │ PT2     1.36207   0.0035921
   3 │ PT3     1.8415    0.0035921
   4 │ PT4     2.05422   0.0035921
   5 │ PT5     2.33997   0.00299529
   6 │ PT6     3.84492   0.00392726
   7 │ PT7     4.49897   0.00392726
   8 │ PT8     5.10222   0.00392726
   9 │ PT9     5.42607   0.00392726
  10 │ PT10    5.98329   0.00392726
  11 │ PT11    6.50716   0.00392726
  12 │ PT12    6.75481   0.00392726

julia> xt_data(PressureTrace("run3/ptrace.lvm"), ptlocs, 
        Species("N2"), Species("Ar"), :PT3 => 10u"psi", :PT6)
((driver = Species(N2, 298.1 K, 4.275e+06 Pa), 
  driven = Species(Ar, 298.1 K, 1.013e+05 Pa), 
 shocked = Species(Ar, 712.5 K, 6.039e+05 Pa), 
      u2 = 429.3754150350808 m s^-1), 
      12×2 DataFrame
 Row │ x         t_shock
     │ Float64   Float64
─────┼────────────────────
   1 │ 0.506412  0.008039
   2 │ 1.36207   0.009323
   3 │ 1.8415    0.010003
   4 │ 2.05422   0.010307
   5 │ 2.33997   0.010695
   6 │ 3.84492   0.012797
   7 │ 4.49897   0.013711
   8 │ 5.10222   0.014637
   9 │ 5.42607   0.015019
  10 │ 5.98329   0.015796
  11 │ 6.50716   0.016535
  12 │ 6.75481   0.016881)

```
"""
function xt_data(ptrace::PressureTrace, ptlocs::DataFrame, driver::Fluid, driven::Fluid, trigger::Pair, ref_PT::Symbol)
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

    return (;states, shock=(;W₀, M₀), xt=DataFrame(:x => ptlocs[!, :x_m], :t_shock => t_shock))
end

"""
An alternate `ptrace` syntax allows for Base Julia-only arguments, namely filepaths, strings, and symbols.

# Examples
```julia-repl
julia> xt_data("run1/ptrace.lvm", "PTlocations.csv", "N2", "Ar", :PT3 => 10, :PT6)
((driver = Species(N2, 298.1 K, 1.633e+06 Pa), 
  driven = Species(Ar, 298.1 K, 1.013e+05 Pa), 
 shocked = Species(Ar, 563.7 K, 4.085e+05 Pa), 
      u2 = 316.0718112278267 m s^-1), 
  12×2 DataFrame
   Row │ x         t_shock
       │ Float64   Float64
  ─────┼────────────────────
     1 │ 0.506412  0.007734
     2 │ 1.36207   0.009218
     3 │ 1.8415    0.010002
     4 │ 2.05422   0.010372
     5 │ 2.33997   0.010854
     6 │ 3.84492   0.013367
     7 │ 4.49897   0.014464
     8 │ 5.10222   0.023391
     9 │ 5.42607   0.016024
    10 │ 5.98329   0.016962
    11 │ 6.50716   0.017844
    12 │ 6.75481   0.018261)
```
"""
function xt_data(ptrace_path, ptloc_path, drivergas, drivengas, trigger, ref_PT; 
    ptrace_filter = (median_size = 15, lowpass_mul=100, filt=Butterworth(2)))

    ptrace = PressureTrace(ptrace_path, ptrace_filter)
    ptlocs = CSV.File(ptloc_path) |> DataFrame
    driver = Species(drivergas)
    driven = Species(drivengas)
    trigger_psi = first(trigger) => last(trigger)*u"psi"
    xt_data(ptrace, ptlocs, driver, driven, trigger_psi, ref_PT)
end
end # module