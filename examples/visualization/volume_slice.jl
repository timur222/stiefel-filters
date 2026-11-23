include(joinpath(@__DIR__, "../../src/visualization/plot_two_surfaces.jl"))

plot_volume_slice_transformed(2; N=4, nx=50, ny=50, pixels=500, #??? Возвращена сетка, с которой график раньше строился быстро.
# !!! plot_volume_slice_transformed(1_000; N=4, nx=5, ny=5, pixels=500,
                              graph_name="slice_1_000_000") #??? Измените имя перед запуском.
