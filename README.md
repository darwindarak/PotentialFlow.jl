# VortexModel

*a scaffolding for building vortex models*

| Documentation | Build Status |
|:---:|:---:|
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://darwindarak.github.com/VortexModel.jl) | [![Build Status](https://travis-ci.org/darwindarak/VortexModel.jl.svg?branch=master)](https://travis-ci.org/darwindarak/VortexModel.jl) [![Build status](https://ci.appveyor.com/api/projects/status/7fdnt699n6l0g92q/branch/master?svg=true)](https://ci.appveyor.com/project/darwindarak/vortexmodel-jl/branch/master) [![codecov](https://codecov.io/gh/darwindarak/VortexModel.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/darwindarak/VortexModel.jl) |

## Installation

This package requires Julia `0.6-` and above.
It is not a registered package, so it should be installed with:
```julia
julia> Pkg.clone("git@github.com:darwindarak/VortexModel.jl.git")
```
Since it is still under heavy development, you should run
```julia
julia> Pkg.test("VortexModel") # might take some time
```
to make sure things are working as intended and
```julia
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

Examples can be found in the [documentation](https://darwindarak.github.io/VortexModel.jl) and the [Jupyter notebooks](https://github.com/darwindarak/VortexModel.jl/tree/master/examples).
