# @autodoc

struct Prms
    years::Int64
    days::Int64                    # num days (run time)
    nrec::Int64                    # num timepoints to record
    dt::Float64                    # length of one time-step
    nt::Int64                      # num days / num timesteps
    nn::Int64                      # number of inorganic N pools
    nc::Int64                      # number of inorganic C pools
    np::Int64                      # number of phytoplankton
    nz::Int64                      # number of zooplankton
    nb::Int64                      # number of bacteria
    ndn::Int64                     # number of DOM pools
    ndc::Int64                     # number of DOC pools
    nv::Int64                      # number of viral pools 
    CNr::Int64                     # carbon to nitrogen ratio
    nIC::Array{Float64,1}          # initial condition for n at time 0
    cIC::Array{Float64,1}          # initial condition for c at time 0
    pIC::Array{Float64,1}          # initial condition for p at time 0
    zIC::Array{Float64,1}          # initial condition for z at time 0
    bIC::Array{Float64,1}          # initial condition for b at time 0
    dnIC::Array{Float64,1}         # initial condition for dom at time 0
    dcIC::Array{Float64,1}         # initial condition for doc at time 0
    vIC::Array{Float64,1}          # initial condition for v at time 0
    vmax_i::Array{Float64,1}       # max uptake rate overall for p_i
    vmax_ij::Array{Float64,2}      # max uptake rate of phyto j on n
    Kp_i::Array{Float64,1}         # half saturation rate for p_i 
    Kp_ij::Array{Float64,2}        # half saturation rate of phyto j on n 
    m_lp::Array{Float64,1}         # linear mort of p (1/d) 
    m_qp::Array{Float64,1}         # quadratic mort of p (m3/mmol/d)
    CMp::Array{Bool,2}             # consumption matrix: nn X np
    Fg_p::Array{Float64, 1}        # fraction of proteome devoted to growth for each P
    umax_i::Array{Float64,1}       # max uptake rate overall for d_i
    umax_ij::Array{Float64,2}      # max uptake rate of bacteria j on d_i
    Km_i::Array{Float64,1}         # half saturation rate overall for d_i 
    Km_ij::Array{Float64,2}        # half saturation rate of bacteria j on d_i 
    y_ij::Array{Float64,2}         # yield rate of bacteria j on d_i
    m_lb::Array{Float64,1}         # linear mort of b
    m_qb::Array{Float64,1}         # quadratic mort of b
    CM::Array{Bool,2}              # consumption matrix: nd X nb
    Fg_b::Array{Float64, 1}        # fraction of proteome devoted to growth for each B
    g_max::Array{Float64,1}        # max grazing rate of z
    K_g::Array{Float64,1}          # half saturation rate of z
    γ::Array{Float64,1}            # fraction of assimilation (assimilation efficiency) for z 
    m_lz::Array{Float64,1}         # linear mort of z
    m_qz::Array{Float64,1}         # quadratic mort of z
    GrM::Array{Bool,2}             # grazing matrix: nz x (np + nb)
    vly::Float64                   # viral lysis rate - L/(virus * day)
    vbs::Int64                     # viral burst size
    vde::Float64                   # viral decay rate - 1/day
    VM::Array{Bool,2}              # lysis matrix: nv * (np + nb)
    fsaven::String                 # save file name
    om_dist_mort::Array{Float64,1} # Distribution of OM from mortality
    om_dist_lys::Array{Float64,1}  # Dist of OM from viral lysis 
    om_dist_vde::Array{Float64,1}  # Dist of OM from viral decay
    rsource::Float64               # OM supply rate (mmol N/m3/day)
    rsink::Float64                 # OM sinking rate (1/day)
    pulse::Int64                   # 0 = no nutrient pulsing, 1 = periodic pulse, 2 = semi-stochastc pulse
    graze::Int64                   # 1 = explicit grazing, 2 = implicit grazing
    lysis::Int64                   # 1 = explicit viral lysis, 2 = implicit viral lysis (incorporated into quadratic mort term)
    carbon::Int64                  # 1 = carbon pools included, 2 = not included
end

