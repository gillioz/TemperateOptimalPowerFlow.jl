module TemperateOptimalPowerFlowGurobiExt

using TemperateOptimalPowerFlow
using Gurobi
import MathOptInterface as MOI


function __init__()
    global gurobi_env = Gurobi.Env()
end

function TemperateOptimalPowerFlow.get_optimizer()
    MOI.instantiate(() -> Gurobi.Optimizer(gurobi_env))
end

function TemperateOptimalPowerFlow.get_silent_optimizer()
    MOI.instantiate(MOI.OptimizerWithAttributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => 0))
end

end # module TemperateOptimalPowerFlowGurobiExt
