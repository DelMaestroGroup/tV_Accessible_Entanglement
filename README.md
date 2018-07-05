# tV_Operational_Entanglement

A Julia exact diagonalization code for the [tV model](https://en.wikipedia.org/wiki/Bose%E2%80%93Hubbard_model) with a focus on operational entanglement entropy under a spatial bipartition [Wiseman & Vaccaro, 2003](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.91.097902) & [Barghathi et al., 2018](https://arxiv.org/pdf/1804.01114.pdf). 

The basis is enumerated using the method of [Szabados et al., 2011](http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf).
See also [Zhang et al., 2011](http://arxiv.org/pdf/1102.4006v1.pdf).


## Requirements

* [ArgParse](https://github.com/carlobaldassi/ArgParse.jl) (`Pkg.add("ArgParse")`)
* [JeszenszkiBasis](https://github.com/0/JeszenszkiBasis.jl) (`Pkg.clone("https://github.com/0/JeszenszkiBasis.jl.git")`)

### Optional

* [LsqFit](https://github.com/JuliaOpt/LsqFit.jl) (`Pkg.add("LsqFit")`)
* [Qutilities](https://github.com/0/Qutilities.jl) (`Pkg.clone("https://github.com/0/Qutilities.jl.git")`)


## Examples

* `julia tV_main.jl --help`
* `julia tV_main.jl --out output.dat --ee 1 --site-max 1 4 2`
* `julia tV_main.jl --out output.dat --pbc --u-log --u-min 2 --u-max 2 --u-num 1 --ee 1 --site-max 1 --u-num 1 --probs 4 2`
