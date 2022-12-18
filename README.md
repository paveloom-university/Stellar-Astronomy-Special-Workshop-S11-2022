### Notices

#### Mirrors

Repository:
- [Codeberg](https://codeberg.org/paveloom-university/Stellar-Astronomy-Special-Workshop-S11-2022)
- [GitHub](https://github.com/paveloom-university/Stellar-Astronomy-Special-Workshop-S11-2022)
- [GitLab](https://gitlab.com/paveloom-g/university/s11-2022/stellar-astronomy-special-workshop)

#### Prerequisites

Make sure you have installed

- [Julia](https://julialang.org)
- [Tectonic](https://tectonic-typesetting.github.io)
- TeX Live:
    - Packages (Fedora Linux):
        - `poppler-utils`
        - `texlive-pgfplots`
        - `texlive-standalone`
        - `texlive-xetex`

#### Run

First, instantiate the project with

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Then, run one of the following

```bash
# Run without a daemon
julia --project=. scripts/compute.jl
julia --project=. scripts/compute.jl --postfix "'Custom postfix'"

# Run with a daemon
./julia.bash scripts/compute.jl
./julia.bash scripts/compute.jl --postfix "'Custom postfix'"
```

Use the `--help` flag to get help information.

To kill the daemon, run `./julia.bash kill`.

#### Report

To compile the report, run

```bash
tectonic -X compile report/report.tex
```
