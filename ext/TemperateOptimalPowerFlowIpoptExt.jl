module TemperateOptimalPowerFlowIpoptExt

using TemperateOptimalPowerFlow
using Ipopt
import MathOptInterface as MOI


function TemperateOptimalPowerFlow.get_optimizer()
    Ipopt.Optimizer
end

function TemperateOptimalPowerFlow.get_silent_optimizer()
    MOI.instantiate(MOI.OptimizerWithAttributes(get_optimizer(), "print_level" => 0))
end

end # module TemperateOptimalPowerFlowIpoptExt