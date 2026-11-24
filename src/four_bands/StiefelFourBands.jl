module StiefelFourBands

using LinearAlgebra
using Printf
using TOML

export FilterModel, solve_polygon, synthesize, evaluate, wall_values, elliptic_constants,
       ThetaCurve, theta_coordinate, theta_controlled, riemann_constant_residual,
       θ, branch_characteristic, theta_digits,
       PaperFilter, paper_period_matrix, paper_branchpoints, paper_filter, paper_fourbands,
       ZolotarevFourBands, zolotarev_fourbands, FoldedFourBands, folded_fourbands,
       folded_period_matrix, folded_from_periods, PickFourBands, pick_fourbands, pick_branchpoints,
       pick_period_matrix, pick_from_periods, CubicFourBands, cubic_fourbands,
       filter_from_periods, period_parameters, read_period_filter, source_genus,
       save_filter, load_filter, plot_filter, oscillating_threebands, welded_value

include("geometry.jl")
include("inverse.jl")
include("polygon.jl")
include("theta.jl")
include("evaluate.jl")
include("digits.jl")
include("paper_filter.jl")
include("paper_inverse.jl")
include("period_input.jl")
include("other_bands.jl")
include("pick_periods.jl")
include("folded_periods.jl")
include("cubic_inverse.jl")
include("period_file.jl")
include("storage.jl")
include("plot.jl")
include("../three_bands/weighted_welding.jl")
include("../three_bands/oscillating_threebands.jl")

end
