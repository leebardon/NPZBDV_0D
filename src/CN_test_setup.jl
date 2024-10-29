
include("CN_test_model.jl")
include("test_helper_funcs.jl")
# using Plots


#------------------------------------------------------------------------------------------------------------#
#                                           CORE PARAMS 
#------------------------------------------------------------------------------------------------------------#
nn = 1
nc = 1
ndon = 1
ndoc = 1
nb = 2

years = 1
days = years * 366
nrec = days * 20
dt = 0.01
nt = Int(days / dt)

# -----------------------------------------------------------------------------------------------------------#
#                                    HETEROTROPHIC BACTERIA PARAMS
#------------------------------------------------------------------------------------------------------------#
# CM = [1.0;;]
CM = [1 1 
      0 0 ] 

y_i = ones(ndon) * 0.3
y_ij = broadcast(*, y_i, CM)

# Fg_b = [1.0]
Fg_b = [1.0, 1.0]
Fa_b = 1.0 .- Fg_b

# umax_i = [1.0]      # Fg_b and umax_i are dummy vals - already provided from Emily's previous work (trade-off applied)
Vmax_i = [1.0, 1.0]  
Km_i = Vmax_i ./ 10
# vmax_ij = [32.0;;]
Vmax_ij =  [32.0  4.0
            0.0  0.0 ] 
# Km_ij = [1.54;;]
Km_ij = [1.54  0.124
         0.0  0.0 ]


# -----------------------------------------------------------------------------------------------------------#
#                                     MORTALITY RATES (mmol/m3/day)
#------------------------------------------------------------------------------------------------------------#
m_lb = ones(nb) * 0.01
m_qb = ones(nb) * 1.0
om_dist_mort = [1.0]


# -----------------------------------------------------------------------------------------------------------#
#   INITIAL CONDITIONS
#------------------------------------------------------------------------------------------------------------#
nIC = ones(nn) * 10.0
cIC = ones(nc) * 10.0
bIC = ones(nb) * 0.1
donIC = ones(ndon) * 1.0
docIC = ones(ndoc) * 5.0
c2nIC = docIC ./ donIC


# -----------------------------------------------------------------------------------------------------------#
#   INITIATE PARAMS
#------------------------------------------------------------------------------------------------------------#
DON_source = 0.06
DOC_source = 1.0
rsource = 0.0
OM_sink = 1.0
CNr = 5

struct test_params
    CNr::Int64                     # environmental carbon to nitrogen ratio
    years::Int64                   # num years (run time)
    days::Int64                    # num days (run time)
    nrec::Int64                    # num timepoints to record
    dt::Float64                    # length of one time-step
    nt::Int64                      # num days / num timesteps
    nn::Int64                      # number of inorganic N pools
    nc::Int64                      # number of inorganic C pools
    nb::Int64                      # number of bacteria
    ndon::Int64                    # number of DON pools
    ndoc::Int64                    # number of DOC pools
    nIC::Array{Float64,1}          # initial condition for n at time 0
    cIC::Array{Float64,1}          # initial condition for c at time 0
    bIC::Array{Float64,1}          # initial condition for b at time 0
    donIC::Array{Float64,1}        # initial condition for don at time 0
    docIC::Array{Float64,1}        # initial condition for doc at time 0
    c2nIC::Array{Float64,1}        # initial ratio of DOC to DON
    Vmax_i::Array{Float64,1}       # max uptake rate overall for d_i
    Vmax_ij::Array{Float64,2}      # max uptake rate of bacteria j on d_i
    Km_i::Array{Float64,1}         # half saturation rate overall for d_i 
    Km_ij::Array{Float64,2}        # half saturation rate of bacteria j on d_i 
    y_ij::Array{Float64,2}         # yield rate of bacteria j on d_i
    m_lb::Array{Float64,1}         # linear mort of b
    m_qb::Array{Float64,1}         # quadratic mort of b
    CM::Array{Bool,2}              # consumption matrix: nd X nb
    Fg_b::Array{Float64,1}         # fraction of proteome devoted to growth for each B
    om_dist_mort::Array{Float64,1} # Distribution of OM from mortality
    rsource::Float64               # General OM supply rate (mmol N/m3/day)
    DON_source::Float64            # DON supply rate (mmol N/m3/day)
    DOC_source::Float64            # DOC supply rate (mmol C/m3/day)
    OM_sink::Float64               # OM sinking rate (1/day)
end


params = test_params(
    CNr, years, days, nrec, dt, nt,
    nn, nc, nb, ndon, ndoc, nIC, cIC, bIC, donIC, docIC, c2nIC,
    Vmax_i, Vmax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, CM, Fg_b,
    om_dist_mort, rsource, DON_source, DOC_source, OM_sink
)

run_CN_test(params)