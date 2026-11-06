using ArgParse                               # Parse command lines
using BenchmarkTools                         # To benchmark precisely
using Random                                 # Draw random numbers
using HDF5                                   # To use HDF5 files
using FilePathsBase: isdir, mkdir, parent    # To create folders
using MuladdMacro                            # FMA 
using Base.Threads                           # Parallelisation
using LinearAlgebra
using FFTW                                   # FFT library
using SpecialFunctions                       # Used to compute Gamma functions in DOOR

