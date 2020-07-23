module ShockTube

using PyCall 
using Printf
using Unitful
import Base: convert, ==, isequal, hash, getindex, setindex!, haskey, keys, show
export Species, Mixture
export γ, T, p, ρ, ρm, R_specific, soundspeed
export shockjump!, shockjump

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
T(f::Fluid)  = f.T * u"K"
p(f::Fluid)  = f.P * u"Pa"
ρ(f::Fluid)  = f.rho * u"kg/m^3"
ρm(f::Fluid) = f.rhom * u"mol/m^3"
γ(f::Fluid)  = f.isentropic_exponent
R_specific(f::Fluid)  = f.R_specific * u"J/(kg*K)"
soundspeed(f::Fluid) = sqrt(γ(f) * R_specific(f) * T(f))

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
    return gas, u2
end

shockjump(gas, Mach) = shockjump!(copy(gas), Mach)

function driverstate(driver, driven, Ms)
    γ1, γ4 = γ(driven), γ(driver)
    a1, a4 = soundspeed(driven), soundspeed(driver)
    P4 = driven.P * (1 + 2γ1 / (γ1 + 1) * (Ms^2 - 1)) *
                    (1 + a1 / a4 * (γ4 - 1) / (γ1 + 1) * (1/Ms - Ms)) ^ (2γ4 / (1 - γ4))

end

end
