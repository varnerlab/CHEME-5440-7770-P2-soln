using TOML
using PyPlot
using DelimitedFiles

function build_data_dictionary(path_to_parameter_file)
    # load parameter file -
    return TOML.parsefile(path_to_parameter_file)
end

function binding_function(actor,K,n)
    return (((actor/K)^n)/(1+(actor/K)^n))
end

function calculate_kinetic_limit(data_dictionary)

    # get some parameters from the data dictionary -
    kcat = data_dictionary["turnover_number"]*(3600)                # convert to hr^-1
    enzyme_concentration = data_dictionary["enzyme_concentration"]  # units: μmol/L
    F6P = data_dictionary["concentration_F6P"]                      # units: mmol/L
    ATP = data_dictionary["concentration_ATP"]                      # units: mmol/L
    K_ATP = data_dictionary["saturation_constant_ATP"]              # units: mmol/L
    K_F6P = data_dictionary["saturation_constant_F6P"]              # units: mmol/L

    # compute the kinetic limit -
    kinetic_limit = (kcat*enzyme_concentration)*(F6P/(K_F6P+F6P))*(ATP/(K_ATP+ATP))

    return kinetic_limit
end

function calculate_allosteric_correction(ADP,AMP,data_dictionary)

    # compute the correction -

    # f_ADP -
    n = data_dictionary["control_parameters"]["ADP"]["order_parameter"]
    K = data_dictionary["control_parameters"]["ADP"]["binding_parameter"]
    f_ADP = binding_function(ADP,K,n)

    # f_AMP -
    n = data_dictionary["control_parameters"]["AMP"]["order_parameter"]
    K = data_dictionary["control_parameters"]["AMP"]["binding_parameter"]
    f_AMP = binding_function(AMP,K,n)

    # get K's -
    K1 = data_dictionary["control_parameters"]["K1"]
    K2 = data_dictionary["control_parameters"]["K2"]
    K3 = data_dictionary["control_parameters"]["K3"]

    # compute -
    vv = (K1+K2*f_AMP)/(1+K1+K2*f_AMP+K3*f_ADP)

    # return -
    return vv
end

function calculate_rate_function(ADP,AMP,data_dictionary)

    # compute the kinetic limit -
    kinetic_limit = calculate_kinetic_limit(data_dictionary)

    # compute the allosteric correction -
    vv = calculate_allosteric_correction(ADP,AMP,data_dictionary)

    # compute the overall rate -
    rate = kinetic_limit*vv

    # return -
    return rate
end

function compute_AMP_b_array(data_dictionary)

    # load the AMP data -
    AMP_data = readdlm(data_dictionary["data"]["path_AMP_data"])

    # get K1 -
    K1 = data_dictionary["control_parameters"]["K1"]

    # calculate the kinetic_limit -
    kinetic_limit = calculate_kinetic_limit(data_dictionary)

    # init -
    b_vector = Float64[]

    # size of data -
    (nr,nc) = size(AMP_data)
    for sample_index = 1:nr
        AMP_value = AMP_data[sample_index,1]
        rate_value = AMP_data[sample_index,2]

        # compute alpha -
        alpha_value = (rate_value/kinetic_limit)

        # compute b-value -
        b_value = (K1-(1+K1)*alpha_value)/(alpha_value - 1)

        # cache -
        push!(b_vector, b_value)
    end

    return (b_vector, AMP_data)
end

# load data_dictionary -
path_to_parameter_file = "./Parameters.toml"
dd = build_data_dictionary(path_to_parameter_file)

# calculate the kinetic limit -
kinetic_limit = calculate_kinetic_limit(dd)

# compute_AMP_b_array -
(b_array, AMP_data) = compute_AMP_b_array(dd)

# run the model for AMP -
AMP_array = collect(range(0.0,1.0,length=100))
ADP = 0.0
rate_AMP = Float64[]
for AMP_value in AMP_array
    r_value = calculate_rate_function(ADP,AMP_value,dd)
    push!(rate_AMP,r_value)
end

# simulation -
plot(AMP_array,rate_AMP)
errorbar(AMP_data[:,1],AMP_data[:,2],yerr=AMP_data[:,3],fmt="ro")

# labels -
xlabel("3-5-AMP [mM]",fontsize=16)
ylabel(L"Rate [$\mu$M/h]",fontsize=16)
