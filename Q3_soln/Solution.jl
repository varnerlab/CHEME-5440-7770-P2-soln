include("Include.jl")

# load the parameters -
function load_parameter_dictionary(path_to_dictionary::String)::Dict{AbstractString,Any}
    return TOML.parsefile(path_to_dictionary)
end

function calculate_translation_gain(parameter_dictionary::Dict{AbstractString,Any})::Float64

    # constants -
    AVN = 6.02e23

    # ribsosome_copy_number=26000.0   # units: copies/cell source: BIND:101441
    # ribosome_elongation_rate=14.0   # units: aa/s source: BIND:108487
    # characteristic_read_length=333  # units: aa source: problem
    # protein_read_length=300         # units: aa source: problem
    # protein_half_life=24            # units: hr source: problem
    # translation_sat_constant=200.0  # units: muM source:  problem
    # mass_ecoli_cell=4.3e-13         # units: g/cell source: problem
    # volume_ecoli_cell=1e-15         # units: L/cell source: problem
    # ecoli_doubling_time=40          # units: min source: problem   
    # ecoli_fraction_water=0.70       # units: dimensionless source: problem

    # get stuff from the parameter dictionary -
    V = parameter_dictionary["volume_ecoli_cell"]                   # L/cell
    M = parameter_dictionary["mass_ecoli_cell"]                     # g/cell
    WF = parameter_dictionary["ecoli_fraction_water"]               # dimensionless
    KL = parameter_dictionary["translation_sat_constant"]           # muM
    RL = parameter_dictionary["ribsosome_copy_number"]              # copies/cell
    L = parameter_dictionary["protein_read_length"]                 # aa
    Lo = parameter_dictionary["characteristic_read_length"]         # aa
    eL = parameter_dictionary["ribosome_elongation_rate"]           # aa/s
    tau_L = parameter_dictionary["tau_L"]                           # dimensionless
    dt = parameter_dictionary["ecoli_doubling_time"]                # min
    protein_half_life = parameter_dictionary["protein_half_life"]   # hr

    # convert KL -
    KL = KL*(V)*(1/M)*(1/(1-WF))            # mumol/gDW

    # compute Vmax_L -
    RL = RL*(1/AVN)*(1/M)*(1/(1-WF))*1e6    # mumol/gDW`
    kE = eL*(1/Lo)*3600                     # hr^-1
    LF = (Lo/L)                             # dimensionless
    Vmax_L = kE*LF*RL                       # mumol/gDW-hr

    # degradation term -
    mu = (log(2)/dt)*60                     # hr^-1
    kd_L = -1*(log(0.5)/protein_half_life)  # hr^-1
    degradation_term = (mu+kd_L)

    # compute the gain -
    translation_gain = Vmax_L/(degradation_term*tau_L*KL)

    # return -
    return translation_gain
end