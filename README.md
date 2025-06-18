[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/MolSimToolkit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/MolSimToolkit.jl/dev)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/MolSimToolkit.jl/actions)
[![codecov](https://codecov.io/gh/m3g/MolSimToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/m3g/MolSimToolkit.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# MyDevMolSimToolkit

[MolSimToolkit.jl](https://github.com/m3g/MolSimToolkit.jl) provides a set of tools to 
analyse molecular dynamics simulations, and a framework for the development of custom
analysis tools. This repository is a experimental version of the package for Work in Progress functionalities

## Installation

```julia
import Pkg; Pkg.add("MolSimToolkit")
```

## Using

```julia
using MolSimToolkit
```

## Recommended

Also install `PDBTools`:

```julia
import Pkg; Pkg.add("PDBTools")
using PDBTools
```







