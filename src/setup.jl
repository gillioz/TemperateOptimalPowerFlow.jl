export import_model, split_nondispatchable!, get_ordered_gen_ids, setup, scale_down


Network = Dict{String, Dict{String, Any}}


function import_model(file::String) :: Dict{String, Dict{String, Any}}
    network = JSON.parsefile(file)

    # TODO: check the presence of critical elements in the network

    return Dict(key => network[key] for key in ["bus", "gen", "load", "branch"])
end


function split_nondispatchable!(network::Network, nondisp_ids::Vector{String})

    nondisp_ids = intersect(nondisp_ids, keys(network["gen"]))
    if "gen_nd" ∉ keys(network)
        network["gen_nd"] = Dict{String, Any}()
    end
    for id ∈ nondisp_ids
        network["gen_nd"][id] = network["gen"][id]
        delete!(network["gen"], id)
    end
    nothing
end


function get_ordered_ids(network::Network, key::String)
    return string.(sort(parse.(Int, keys(network[key]))))
end


get_ordered_gen_ids(network::Network) = get_ordered_ids(network, "gen")


function setup(output_directory::String, network::Network, loads::AbstractArray{<:Real},
        gen_costs::AbstractArray{<:Real}, nondispatch_ids::Vector{String} = String[],
        nondispatch_series::AbstractArray{<:Real} = Real[];
        overwrite::Bool = false)

    if isdir(output_directory)
        if overwrite
            @info "Overwriting directory '$(output_directory)'" _group = ""
            rm(output_directory, recursive = true)
            mkdir(output_directory)
        else
            throw(ErrorException("Cannot overwrite existing directory '$(output_directory)'"))
        end
    else
        @info "Creating directory '$(output_directory)'" _group = ""
        mkdir(output_directory)
    end

    @info "Loading network" _group = ""
    N_buses = length(network["bus"])
    N_loads = length(network["load"])
    N_gens = length(network["gen"])
    N_lines = length(network["branch"])
    N_nondispatch = length(network["gen_nd"])

    T = size(loads, 2)

    @assert size(loads) == (N_loads, T)
    @assert size(gen_costs) == (N_gens, T)
    @assert length(nondispatch_ids) == N_nondispatch
    @assert size(nondispatch_series) == (N_nondispatch, T)

    bus_ids = get_ordered_ids(network, "bus")
    gen_ids = get_ordered_ids(network, "gen")
    load_ids = get_ordered_ids(network, "load")
    line_ids = get_ordered_ids(network, "branch")
    nondispatch_ids = get_ordered_ids(network, "gen_nd")

    bus_ids_map = Dict(id => i for (i, id) in enumerate(bus_ids))
    gen_ids_map = Dict(id => i for (i, id) in enumerate(gen_ids))

    @info "Computing susceptance matrix" _group = ""
    # incidence matrix
    lines_from_bus = [bus_ids_map[string(network["branch"][id]["f_bus"])] for id ∈ line_ids]
    lines_to_bus = [bus_ids_map[string(network["branch"][id]["t_bus"])] for id ∈ line_ids]
    M = sparse([lines_from_bus; lines_to_bus], [1:N_lines; 1:N_lines],
        [-ones(N_lines); ones(N_lines)], N_buses, N_lines)
    b = [network["branch"][id]["br_b"] for id ∈ line_ids];
    B = M * Diagonal(b)
    DataDrop.store_matrix("$(output_directory)/susceptance.h5", b)

    @info "Computing PTDF matrix (this may take some time)" _group = ""
    # pseudo-inverse
    Bi = pinv(Matrix(B))
    # PTDF matrix
    L = Diagonal(b) * Bi
    DataDrop.store_matrix("$(output_directory)/PTDF_matrix.h5", L)

    @info "Distributing generator on buses" _group = "" _group = ""
    gen_buses = [bus_ids_map[string(network["gen"][id]["gen_bus"])] for id ∈ gen_ids]
    A_gen = sparse(gen_buses, 1:N_gens, ones(N_gens), N_buses, N_gens)
    DataDrop.store_matrix("$(output_directory)/A_gen.h5", A_gen)

    @info "Distributing loads on buses" _group = "" _group = ""
    load_buses = [bus_ids_map[string(network["load"][id]["load_bus"])] for id ∈ load_ids]
    A_load = sparse(load_buses, 1:N_loads, ones(N_loads), N_buses, N_loads)
    DataDrop.store_matrix("$(output_directory)/A_load.h5", A_load)

    @info "Distributing non-dispatchable generators on buses" _group = "" _group = ""
    A_nondispatch = sparse(
        indexin([string(network["gen_nd"][id]["gen_bus"]) for id in nondispatch_ids], bus_ids),
        1:N_nondispatch, ones(N_nondispatch), length(bus_ids), N_nondispatch)
    DataDrop.store_matrix("$(output_directory)/A_nondispatch.h5", A_nondispatch)
    DataDrop.store_matrix("$(output_directory)/P_nondispatch.h5", nondispatch_series)

    @info "Computing max capacity for each generator" _group = ""
    gen_pmax = [convert(Float64, network["gen"][id]["pmax"]) for id ∈ gen_ids]
    DataDrop.store_matrix("$(output_directory)/P_max.h5", gen_pmax)

    @info "Computing expected production for each generator" _group = ""
    gen_pexp = [convert(Float64, network["gen"][id]["pexp"]) for id ∈ gen_ids]
    DataDrop.store_matrix("$(output_directory)/P_exp.h5", gen_pexp)

    @info "Computing thermal limits" _group = ""
    thermal_limits = [convert(Float64, network["branch"][id]["rate_a"]) for id ∈ line_ids]
    DataDrop.store_matrix("$(output_directory)/thermal_limits.h5", thermal_limits)

    @info "Computing quadratic line costs" _group = ""
    line_costs = 1. ./ thermal_limits
    LA = L * A_gen
    quadratic_cost = LA' * (line_costs .* LA)
    DataDrop.store_matrix("$(output_directory)/quadratic_cost.h5", quadratic_cost)

    @info "Computing linear line costs" _group = ""
    injections = A_nondispatch * nondispatch_series - A_load * loads
    linear_cost = 2 * LA' * (line_costs .* L * injections)
    DataDrop.store_matrix("$(output_directory)/linear_line_cost.h5", linear_cost)

    @info "Computing linear generation costs" _group = ""
    DataDrop.store_matrix("$(output_directory)/linear_gen_cost.h5", gen_costs)

    @info "Computing total load constraints" _group = ""
    DataDrop.store_matrix("$(output_directory)/P_load.h5", loads)
    total_load = sum(loads, dims=1)[1, :]
    total_nondispatch = sum(nondispatch_series, dims=1)[1, :]
    DataDrop.store_matrix("$(output_directory)/P_total.h5", total_load - total_nondispatch)

    @info "Computing ramp constraints" _group = ""
    ramp_gen_ids = [id for id in gen_ids if network["gen"][id]["max_ramp_rate"] > 0]
    N_ramp = length(ramp_gen_ids)
    A_ramp = sparse(1:N_ramp, [gen_ids_map[id] for id ∈ ramp_gen_ids], ones(N_ramp), N_ramp, N_gens)
    DataDrop.store_matrix("$(output_directory)/A_ramp.h5", A_ramp)
    ramp_max = [network["gen"][id]["max_ramp_rate"] for id in ramp_gen_ids]
    DataDrop.store_matrix("$(output_directory)/ramp_max.h5", ramp_max)
end


function scale_down_file(file::String, new_timesteps::Int)
    data = DataDrop.retrieve_matrix(file)
    data_size = size(data)
    dims = length(data_size)
    old_timesteps = data_size[end]

    @info "Scaling down file '$file' from $(old_timesteps) to $(new_timesteps) time steps" _group = ""

    if old_timesteps % new_timesteps != 0
        throw(ArgumentError("The number of time steps ($(old_timesteps)) " *
            "cannot be divided into $(new_timesteps) parts"))
    end
    n = old_timesteps ÷ new_timesteps

    if dims == 1
        data = dropdims(sum(reshape(data, (n, new_timesteps)), dims=1), dims=1) ./ n
    elseif dims == 2
        data = dropdims(sum(reshape(data, (:, n, new_timesteps)), dims=2), dims=2) ./ n
    else
        throw(ArgumentError("Only arrays with up to two dimensions are supported"))
    end

    rm(file)
    DataDrop.store_matrix(file, data)

    nothing
end


function scale_down(data_directory::String, new_timesteps::Int)
    scale_down_file("$(data_directory)/P_total.h5", new_timesteps)
    scale_down_file("$(data_directory)/linear_line_cost.h5", new_timesteps)

    if isfile("$(data_directory)/linear_gen_cost.h5")
        scale_down_file("$(data_directory)/linear_gen_cost.h5", new_timesteps)
    end

    nothing
end


# TODO: recombination of dispatchable and non-dispatchable series in one file
# TODO: export results as dataframe in csv format (both gens and lines)
