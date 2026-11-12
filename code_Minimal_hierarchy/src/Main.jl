##################################################
include("Packages.jl")       # Loading the packages
include("Args.jl")           # Parsing the command-line
include("Arrays.jl")         # All arrays for the grid for Landau eq
include("Utils.jl")          # integration and computation tools
include("Constants.jl")      # Global physical constants
include("Observables.jl")    # To compute observables
include("Gradients.jl")      # Various utility functions
include("Interaction.jl")    # f(u) = interaction term
include("FFT.jl")            # FFT part
include("Initialisation_F.jl") # Initilaisation of F
include("Integration.jl")    # Integration of F
include("Kernel.jl")         # Computing the interaction kernel once and for all
# include("GC.jl")             # Test of the Gallavotti-Cohen fluctuation theorem


