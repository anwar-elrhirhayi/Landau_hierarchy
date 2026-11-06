##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--G"
    help = "Gravitational constant"
    arg_type = Float64
    default = 1.0
    "--Mtot"
    help = "Total mass"
    arg_type = Float64
    default = 1.0
    "--Nv"
    help = "Number of grid discretization"
    arg_type = Int64
    default = 25
    "--Npart"
    help = "Number of particles"
    arg_type = Int64
    default = 100
    "--kmin"
    help = "Minimum resonance number"
    arg_type = Int64
    default = 10
    "--kmax"
    help = "Maximum resonance number"
    arg_type = Int64
    default = 100
    "--vmax"
    help = "Maximum velocity cutoff in the grid: VMAX ~ 5 sigma_0"
    arg_type = Float64
    default = 5.0
    "--dt"
    help = "Integration timestep"
    arg_type = Float64
    default = 5e-2
    "--Nsteps"
    help = "Number of timesteps"
    arg_type = Int64
    default = 20
    "--hierarchy"
    help = "Solve Lypunov hierarchy for F, S2 or S3"
    arg_type = String
    default = "S2"
    "--IC"
    help = "Initial condition for F: Maxwellian or Door"
    arg_type = String
    default = "MAXW"
    "--method"
    help = "Method used to compute the L_v operator: FD (finite difference) or FFT"
    arg_type = String
    default = "FFT"
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################
const G           = parsed_args["G"]             # Gravitational constant
const MASS_TOT    = parsed_args["Mtot"]          # Total mass
const NV          = parsed_args["Nv"]            # Number of points discretization
const NPART       = parsed_args["Npart"]         # Number of particles
const KMIN        = parsed_args["kmin"]          # Minimum resonance number
const KMAX        = parsed_args["kmax"]          # Maximum resonance number
const VMAX        = parsed_args["vmax"]          # Maximum velocity cutoff in the grid
const DT          = parsed_args["dt"]            # Integration timestep
const NSTEPS      = parsed_args["Nsteps"]        # Number of timesteps
const HIERARCHY   = parsed_args["hierarchy"]     # Compute the hierarchy 
const IC          = parsed_args["IC"]            # Initial Condition
const METHOD      = parsed_args["method"]        # Method used to compute the L_v operator
#####
if !(0 < KMIN <= KMAX)
    error("Interaction should satisfy: 0 < KMIN <= KMAX")
end
#####
if  (HIERARCHY != "F") &&
    (HIERARCHY != "S2") &&
    (HIERARCHY != "S3")  
    #####
    error("Unsupported ANSWER")
    #####
end
#####