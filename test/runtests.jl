using Test
using PowerModels
using TemperateOptimalPowerFlow
using DataFrames
using Ipopt
using MathOptInterface


@testset "TemperateOptimalPowerFlow.jl" begin

    @testset "test_network_import_and_preparation" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # verify presence of branches
        @test "branch" in keys(network)
        branches = network["branch"]
        @test isa(branches, Dict)
        @test length(branches) > 0
        # verify that branches do not have a "cost" property
        for (key, branch) in network["branch"]
            @test "cost" ∉ keys(branches)
        end
        # add a branch cost
        cost = 123.4
        add_line_costs!(network, cost)
        # verify that all branches now have a "cost" property
        for (key, branch) in network["branch"]
            @test "cost" in keys(branch)
            @test branch["cost"] == cost
        end
    end

    @testset "test_create_list_of_loads" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # create a list of loads
        @test issetequal(create_list_of_loads(network), ["1", "2", "3"])
    end

    @testset "test_create_list_of_gens" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # create a list of loads
        @test issetequal(create_list_of_gens(network), ["1", "2", "3"])
    end

    @testset "test_get_loads_info" begin
        # import PanTaGruEl
        network = parse_file("pantagruel.json")
        @test isa(network, Dict)
        # create a list of loads
        list_of_loads = create_list_of_loads(network)
        @test length(list_of_loads) > 0
        # get info about the load's country and population
        df = get_loads_info(network, list_of_loads, ["country", "load_prop"])
        @test isa(df, DataFrame)
        # for each country, check that the 'load_prop' adds up to one
        countries = unique(df[:, "country"])
        for country in countries
            @test abs(1.0 - sum(df[df[:, "country"] .== country, "load_prop"])) < 1.0e-6
        end
    end

    @testset "test_get_gens_info" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # create a list of gens
        list_of_gens = create_list_of_gens(network)
        @test length(list_of_gens) == 3
        # get info about the gen's max capacity
        df = get_gens_info(network, list_of_gens, ["pmax"])
        @test isa(df, DataFrame)
        @test issetequal(df[:, "pmax"], [0, 15.0, 20.0])
    end

    @testset "test_assign_loads" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # verify that the model uses per units
        @test network["per_unit"]
        @test network["baseMVA"] == 100.0
        # create a list of loads
        list_of_loads = create_list_of_loads(network)
        @test length(list_of_loads) == 3
        # create loads in MW
        loads = [100.0, 200.0, 300.0]
        assign_loads!(network, list_of_loads, loads)
        for i = 1:3
            @test network["load"][list_of_loads[i]]["pd"] == i
        end
    end

    @testset "test_assign_ramp_max" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # assign types to the 3 generators in the network
        network["gen"]["1"]["type"] = "hydro"
        network["gen"]["2"]["type"] = "nuclear"
        network["gen"]["3"]["type"] = "gas"
        # assign ramp max values for each types 
        assign_ramp_max!(network, 3.0, "nuclear")
        assign_ramp_max!(network, 1.0, ["hydro", "gas"])
        # verify that all generators have now a ramp_max attribute
        for gen in values(network["gen"])
            @test haskey(gen, "ramp_max")
        end
        # verify that the ramp maxima have been properly assigned
        @test network["gen"]["1"]["ramp_max"] == 1.0
        @test network["gen"]["2"]["ramp_max"] == 3.0
        @test network["gen"]["3"]["ramp_max"] == 1.0
    end

    @testset "test_assign_ramp_max_ratio" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # assign types to the 3 generators in the network
        for gen in values(network["gen"])
            gen["type"] = "nuclear"
        end
        # assign ramp max values for each types 
        assign_ramp_max_ratio!(network, 0.1, "nuclear")
        # verify that all generators have now the correct ramp max attribute
        for gen in values(network["gen"])
            @test haskey(gen, "ramp_max")
            @test abs(gen["ramp_max"] - 0.1 * gen["pmax"]) <= 1e-8
        end
    end

    @testset "test_get_optimzer" begin
        optimizer = get_optimizer()
        @test isa(optimizer(), Ipopt.Optimizer)
        silent_optimizer = get_silent_optimizer()
        @test isa(silent_optimizer, MathOptInterface.OptimizerWithAttributes)
    end

    @testset "test_compute_line_rates" begin
        # import a simple model with 3 buses (see PowerModels.jl)
        network = parse_file("case3.m")
        @test isa(network, Dict)
        # get optimiser
        optimizer = get_silent_optimizer()
        # compute line rates
        line_rates = compute_line_rates(network, optimizer)
        @test length(line_rates) == 3
        # verify that all rates are between zero and one
        @test all(values(line_rates) .>= -1e6)
        @test all(values(line_rates) .<= 1 + 1e6)
    end

end
