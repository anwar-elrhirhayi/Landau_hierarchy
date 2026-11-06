# Landau minimal fluctuation hierarchy

This code simulates the relaxation to equilibrium of an anisotropic, long-range interacting system of N particles, using a deterministic *minimal hierarchy* of kinetic equations. On a 2D velocity grid, it solves:

- the **Landau equation** for the mean distribution function $F$, which gives the ensemble-averaged evolution;
- the **Lyapunov equation** for the second-order cumulant $S_2$, which gives the finite-$N$ variance of observables;
- the **Lyapunov equation** for the third-order cumulant $S_3$ (with a closure), which gives the third moment.

Together, these predict the finite-$N$ fluctuations of observables such as the temperatures $T_{xx}$ and $T_{yy}$, without running N-body simulations.

## Installation

The code is written in [Julia](https://julialang.org/downloads/). The required packages can be installed from the Julia REPL with:

```julia
using Pkg
Pkg.add(["ArgParse", "BenchmarkTools", "HDF5", "FFTW", "MuladdMacro",
         "FilePathsBase", "BlockArrays", "LoopVectorization", "Plots"])
```

`Random` and `LinearAlgebra` are part of Julia's standard library and need no installation.

## Usage

From a terminal, move to the code directory:

```bash
cd code_Minimal_hierarchy
```

Each run writes its data as an HDF5 file in the `data/` subdirectory, and saves the corresponding plot in the same folder. The level of the hierarchy is selected with `--hierarchy` (`F`, `S2` or `S3`).

> **Note.** The commands below use a coarse velocity grid (`--Nv`) and a large timestep (`--dt`) so that they run quickly. All three runs integrate up to the same final time (`dt × Nsteps = 1000`). The resulting plots may differ slightly from those in the paper.

### 1. Mean evolution: Landau equation for $F$ (~1 s)

```bash
julia --check-bounds=no run/Run_Hierarchy.jl --Mtot 1.0 --Nv 20 --Npart 100 --kmin 10 --kmax 100 --vmax 6.0 --dt 1 --Nsteps 1000 --hierarchy F --IC MAXW --method FFT
```

Output: `plot_TXX_AVG_HIERARCHY_F_NPART_100_NV_20_IC_MAXW_METHOD_FFT.png`

<img width="600" height="400" alt="plot_TXX_AVG_HIERARCHY_F_NPART_100_NV_20_IC_MAXW_METHOD_FFT" src="https://github.com/user-attachments/assets/2a0149ec-5dfd-4d38-8140-6a6860d4f104" />

### 2. Variance: Lyapunov equation for $S_2$ (~60 s)

```bash
julia --check-bounds=no run/Run_Hierarchy.jl --Mtot 1.0 --Nv 20 --Npart 100 --kmin 10 --kmax 100 --vmax 6.0 --dt 5 --Nsteps 200 --hierarchy S2 --IC MAXW --method FFT
```

Output: `plot_TXX_VAR_HIERARCHY_S2_NPART_100_NV_20_IC_MAXW_METHOD_FFT.png`

<img width="600" height="400" alt="plot_TXX_VAR_HIERARCHY_S2_NPART_100_NV_20_IC_MAXW_METHOD_FFT" src="https://github.com/user-attachments/assets/5dca030e-5bca-4ebe-9a06-35833a2bbf1b" />

### 3. Skewness: Lyapunov equation for $S_3$ (~75 s)

```bash
julia --check-bounds=no run/Run_Hierarchy.jl --Mtot 1.0 --Nv 14 --Npart 100 --kmin 10 --kmax 100 --vmax 6.0 --dt 20 --Nsteps 50 --hierarchy S3 --IC MAXW --method FFT
```

Output: `plot_TXX_M3_HIERARCHY_S3_NPART_100_NV_14_IC_MAXW_METHOD_FFT.png`

<img width="600" height="400" alt="plot_TXX_M3_HIERARCHY_S3_NPART_100_NV_14_IC_MAXW_METHOD_FFT" src="https://github.com/user-attachments/assets/930e0c93-8da9-4bde-829b-5cdb06dc98ec" />

## Contact

Should you encounter any issues, please contact us at anwar.el-rhirhayi@univ-orleans.fr.
