
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
np = 6
nb = 1

years = 1
days = years * 366
nrec = days * 20
dt = 0.01
nt = Int(days / dt)

# -----------------------------------------------------------------------------------------------------------#
#                                       PHYTOPLANKTON PARAMS 
#------------------------------------------------------------------------------------------------------------
CMp = [1 1 1 1 1 1]

Fg_p = [0.1, 0.25, 0.5, 0.68, 0.79, 0.91]      # fraction of proteome optimized to growth
Fa_p = 1.0 .- Fg_p                              # fraction optimized to substrate affintiy

vmax_i = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0]        # max growth rates (per day)
Kp_i = vmax_i ./ 10                              # half saturation of P_i

vmax_ij = set_vmax_ij(nn, np, vmax_i, Fg_p)    # growth rate of P_i on N
Kp_ij = set_Kp_ij(nn, np, Fa_p, CMp, vmax_ij)  # half saturation of P_i on N


# -----------------------------------------------------------------------------------------------------------#
#                                    HETEROTROPHIC BACTERIA PARAMS
#------------------------------------------------------------------------------------------------------------#
CM = [1.0;;]

y_i = ones(ndon) * 0.3
y_ij = broadcast(*, y_i, CM)

Fg_b = [1.0]
Fa_b = 1.0 .- Fg_b

umax_i = [1.0]      # Fg_b and umax_i are dummy vals - already provided from Emily's previous work (trade-off applied)
Km_i = umax_i ./ 10

umax_ij = [32.0;;]

Km_ij = [1.54;;]


# -----------------------------------------------------------------------------------------------------------#
#                                     MORTALITY RATES (mmol/m3/day)
#------------------------------------------------------------------------------------------------------------#
m_lp = ones(np) * 0.1
m_qp = ones(np) * 1.0  # (.1 if explicit grazers, if not, 1)

m_lb = ones(nb) * 0.01
m_qb = ones(nb) * 1.0

om_dist_mort = [1.0]


# -----------------------------------------------------------------------------------------------------------#
#   INITIAL CONDITIONS
#------------------------------------------------------------------------------------------------------------#
nIC = ones(nn) * 5.0
pIC = ones(np) * 0.1
bIC = ones(nb) * 0.01
cIC = ones(nc) * 20.0
donIC = ones(ndon) * 0.1
docIC = ones(ndoc) * 0.1

# Uncomment if only using one B
# bIC[2] = 0.0
# donIC[2] = 0.0
# docIC[2] = 0.0

# -----------------------------------------------------------------------------------------------------------#
#   INITIATE PARAMS
#------------------------------------------------------------------------------------------------------------#
rsource = 0
rsink = 0
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
    np::Int64                      # number of phytoplankton
    nb::Int64                      # number of bacteria
    ndon::Int64                    # number of DON pools
    ndoc::Int64                    # number of DOC pools
    nIC::Array{Float64,1}          # initial condition for n at time 0
    cIC::Array{Float64,1}          # initial condition for c at time 0
    pIC::Array{Float64,1}          # initial condition for p at time 0
    bIC::Array{Float64,1}          # initial condition for b at time 0
    donIC::Array{Float64,1}        # initial condition for don at time 0
    docIC::Array{Float64,1}        # initial condition for doc at time 0
    vmax_i::Array{Float64,1}       # max uptake rate overall for p_i
    vmax_ij::Array{Float64,2}      # max uptake rate of phyto j on n
    Kp_i::Array{Float64,1}         # half saturation rate for p_i 
    Kp_ij::Array{Float64,2}        # half saturation rate of phyto j on n 
    m_lp::Array{Float64,1}         # linear mort of p (1/d) 
    m_qp::Array{Float64,1}         # quadratic mort of p (m3/mmol/d)
    CMp::Array{Bool,2}             # consumption matrix: nn X np
    Fg_p::Array{Float64,1}        # fraction of proteome devoted to growth for each P
    umax_i::Array{Float64,1}       # max uptake rate overall for d_i
    umax_ij::Array{Float64,2}      # max uptake rate of bacteria j on d_i
    Km_i::Array{Float64,1}         # half saturation rate overall for d_i 
    Km_ij::Array{Float64,2}        # half saturation rate of bacteria j on d_i 
    y_ij::Array{Float64,2}         # yield rate of bacteria j on d_i
    m_lb::Array{Float64,1}         # linear mort of b
    m_qb::Array{Float64,1}         # quadratic mort of b
    CM::Array{Bool,2}              # consumption matrix: nd X nb
    Fg_b::Array{Float64,1}        # fraction of proteome devoted to growth for each B
    om_dist_mort::Array{Float64,1} # Distribution of OM from mortality
    rsource::Float64               # OM supply rate (mmol N/m3/day)
    rsink::Float64                 # OM sinking rate (1/day)
end


params = test_params(
    CNr, years, days, nrec, dt, nt,
    nn, nc, np, nb, ndon, ndoc, nIC, cIC, pIC, bIC, donIC, docIC,
    vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, CMp, Fg_p,
    umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, CM, Fg_b,
    om_dist_mort, rsource, rsink
)

run_CN_test(params)