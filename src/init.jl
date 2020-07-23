# global PyObject constants that get initialized at runtime.  We
# initialize them here (rather than via "global foo = ..." in __init__)
# so that their type is known at compile-time.

const PyThermo = PyNULL()
const PyChem = PyNULL()
const PyFluids = PyNULL()
const PyCopy = PyNULL()

function __init__()
    copy!(PyThermo, pyimport_conda("thermo", "thermo", "conda-forge"))
    copy!(PyChem, pyimport("thermo.chemical"))
    copy!(PyFluids, pyimport_conda("fluids", "fluids", "conda-forge"))
    copy!(PyCopy, pyimport("copy"))
end