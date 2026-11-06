# Landau minimal fluctuation hierarchy

Numerical simulations of physical anisotropic N-particle's systems relaxing to equilibrium using the minimal hierarchy.

Should you encounter any issues, please do not hesitate to contact us at anwar.el-rhirhayi@univ-orleans.fr

## Packages
This code is written in Julia. To install Julia: https://julialang.org/downloads/platform/

List of Julia packages required to run the code [to be installed via
"using ..." in the REPL] :

- ArgParse                            
- BenchmarkTools                     
- Random                             
- LinearAlgebra                       
- HDF5
- FFTW
- MuladdMacro                                 
- FilePathsBase: isdir, mkdir, parent
- BlockArrays                       
- LoopVectorization
- Plots

## Run Landau simulations

I. Open the minimal hierarchy folder in the terminal.

```julia
cd code_Minimal_hierarchy
```

II. Generate the data and the figure.

The generated data will be saved as an HDF5 file in the data/ subdirectory, while the corresponding plots will be output in the same data/ folder.
To simulate the average behaviour throughout the Landau equation, run the code with the following command (~1s to execute)
```julia
julia --check-bounds=no run/Run_Hierarchy.jl --Mtot 1.0 --Nv 20 --Npart 100  --kmin 10 --kmax 100 --vmax 6.0 --dt 1 --Nsteps 1000 --hierarchy F --IC MAXW --method FFT
```
This generates the plot named ``plot_TXX_AVG_HIERARCHY_F_NPART_100_NV_20_IC_MAXW_METHOD_FFT.png``:

<img width="600" height="400" alt="plot_TXX_AVG_HIERARCHY_F_NPART_100_NV_20_IC_MAXW_METHOD_FFT" src="https://github.com/user-attachments/assets/2a0149ec-5dfd-4d38-8140-6a6860d4f104" />

**Note**: To accelerate the code here, we purposely choose a small value of NV and a big value of the timestep DT.

## Run the Lyapunov equation for $S_2$

To simulate the Lyapunov equation for S_2, run the following command (~60s to execute)

```julia
julia --check-bounds=no run/Run_Hierarchy.jl --Mtot 1.0 --Nv 20 --Npart 100  --kmin 10 --kmax 100 --vmax 6.0 --dt 5 --Nsteps 200 --hierarchy F --IC MAXW --method FFT
```

This generates the plot named ``plot_TXX_M3_HIERARCHY_S3_NPART_100_NV_14_IC_MAXW_METHOD_FFT.png``:

<img width="600" height="400" alt="plot_TXX_VAR_HIERARCHY_S2_NPART_100_NV_20_IC_MAXW_METHOD_FFT" src="https://github.com/user-attachments/assets/5dca030e-5bca-4ebe-9a06-35833a2bbf1b" />


**Note**: To accelerate the code here, we purposely choose a small value of NV and a big value of the timestep DT.

## Run the Lyapunov equation for $S_3$

To simulate the kinetic equation for S_3, run the following command (~75s to execute)

```julia
julia --check-bounds=no run/Run_Hierarchy.jl --Mtot 1.0 --Nv 14 --Npart 100  --kmin 10 --kmax 100 --vmax 6.0 --dt 20 --Nsteps 50 --hierarchy F --IC MAXW --method FFT
```

This generates the plot named ``plot_TXX_VAR_HIERARCHY_S2_NPART_100_NV_20_IC_MAXW_METHOD_FFT.png``:

<img width="600" height="400" alt="plot_TXX_M3_HIERARCHY_S3_NPART_100_NV_14_IC_MAXW_METHOD_FFT" src="https://github.com/user-attachments/assets/930e0c93-8da9-4bde-829b-5cdb06dc98ec" />

**Note**: To accelerate the code here, we purposely choose a small value of NV and a big value of the timestep DT. This may create a plot that silghtly differs from the paper's one.
