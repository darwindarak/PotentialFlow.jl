# PotentialFlow

*a scaffolding for building 2D inviscid models*

| Documentation | Build Status |
|:---:|:---:|
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://darwindarak.github.com/PotentialFlow.jl) | [![Build Status](https://img.shields.io/travis/darwindarak/PotentialFlow.jl/master.svg?label=linux)](https://travis-ci.org/darwindarak/PotentialFlow.jl) [![Build status](https://img.shields.io/appveyor/ci/darwindarak/PotentialFlow-jl/master.svg?label=windows)](https://ci.appveyor.com/project/darwindarak/potentialflow-jl/branch/master) [![codecov](https://codecov.io/gh/darwindarak/PotentialFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/darwindarak/PotentialFlow.jl) |

## Installation

This package requires Julia `0.6-` and above.
It is not a registered package, so it should be installed with:
```julia
julia> Pkg.clone("git@github.com:darwindarak/PotentialFlow.jl.git")
```
Since it is still under heavy development, you should run
```julia
julia> Pkg.test("PotentialFlow") # might take some time
```
to make sure things are working as intended and
```julia
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

Examples can be found in the [documentation](https://darwindarak.github.io/PotentialFlow.jl) and the [Jupyter notebooks](https://github.com/darwindarak/PotentialFlow.jl/tree/master/examples).
