
export retrieve_gen_results, retrieve_injections, retrieve_loads
export retrieve_line_flows, retrieve_line_rates, retrieve_line_angles
export update_network!


function retrieve_gen_results(data_directory::String, result_file::String = "P_result",
        include_nondispatch::Bool = true)
    P_gen = DataDrop.retrieve_matrix("$(data_directory)/$(result_file).h5")
    if !include_nondispatch
        return P_gen
    end
    P_nondispatch = DataDrop.retrieve_matrix("$(data_directory)/P_nondispatch.h5")
    P_all = vcat(P_gen, P_nondispatch)

    # compute the permutation matrix that re-orders the generators
    gen_ids = DataDrop.retrieve_matrix("$(data_directory)/gen_ids.h5")
    sorted_gen_ids = sort_strings(gen_ids)
    n = length(gen_ids)
    A = sparse([findfirst(==(id), sorted_gen_ids) for id in gen_ids], 1:n, ones(n), n, n)

    return A * P_all
end


function retrieve_injections(data_directory::String, result_file::String = "P_result")
    P_gen = DataDrop.retrieve_matrix("$(data_directory)/$(result_file).h5")
    A_gen = DataDrop.retrieve_matrix("$(data_directory)/A_gen.h5")
    A_nondispatch = DataDrop.retrieve_matrix("$(data_directory)/A_nondispatch.h5")
    P_nondispatch = DataDrop.retrieve_matrix("$(data_directory)/P_nondispatch.h5")
    A_load = DataDrop.retrieve_matrix("$(data_directory)/A_load.h5")
    P_load = DataDrop.retrieve_matrix("$(data_directory)/P_load.h5")
    return A_gen * P_gen + A_nondispatch * P_nondispatch - A_load * P_load
end


function retrieve_line_flows(data_directory::String, result_file::String = "P_result")
    P = retrieve_injections(data_directory, result_file)
    L = DataDrop.retrieve_matrix("$(data_directory)/PTDF_matrix.h5")
    return L * P
end


function retrieve_line_rates(data_directory::String, result_file::String = "P_result")
    flows = retrieve_line_flows(data_directory, result_file)
    thermal_limits = DataDrop.retrieve_matrix("$(data_directory)/thermal_limits.h5")
    return abs.(flows) ./ thermal_limits
end


function retrieve_line_angles(data_directory::String, result_file::String = "P_result")
    flows = retrieve_line_flows(data_directory, result_file)
    susceptance = DataDrop.retrieve_matrix("$(data_directory)/susceptance.h5")
    return abs.(flows) ./ susceptance
end


function retrieve_loads(data_directory::String)
    return DataDrop.retrieve_matrix("$(data_directory)/P_load.h5")
end


function update_network!(network::Dict{String, Any}, loads::AbstractVector{<:Real},
        gens::AbstractVector{<:Real}, lines::AbstractVector{<:Real})
    # update loads
    for (i, id) = enumerate(get_ordered_ids(network, "load"))
        network["load"][id]["pd"] = loads[i]
        network["load"][id]["qd"] = 0
    end
    # update gens
    for (i, id) = enumerate(get_ordered_ids(network, "gen"))
        network["gen"][id]["pg"] = gens[i]
        network["gen"][id]["qg"] = 0
    end
    # update lines
    for (i, id) = enumerate(get_ordered_ids(network, "branch"))
        line = network["branch"][id]
        line["pt"] = lines[i]
        line["pf"] = -lines[i]
        line["qt"] = 0
        line["qf"] = 0
    end
end