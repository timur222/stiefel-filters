#??? На Linux сечение сохраняется в PNG без зависимости от окна GLFW.
if Sys.islinux() && !haskey(ENV, "SHTIFFEL_PLOT_BACKEND") #???
    ENV["SHTIFFEL_PLOT_BACKEND"] = "cairo" #???
end #???
include(joinpath(@__DIR__, "../../src/visualization/gamma_slice.jl")) #???

N = 11 #???
M = 2 #??? По этому M выбираем γ, затем рисуем все M от 1 до floor(N/2).
# !!! x_reference = 4.0
# !!! y_reference = 0.026
# !!! gamma0 = gamma_value(x_reference, y_reference; N=N, M=M)
reference = try #??? Сначала берём ту же точку, что и в исходном примере.
    (; x=4.0, y=0.026, gamma=gamma_value(4.0, 0.026; N=N, M=M)) #???
catch err #???
    @warn "Опорная точка (4, 0.026) не подходит; ищем другую" причина=err #???
    find_gamma_reference(; N=N, M=M, #???
        x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0)) #???
end #???
gamma0 = reference.gamma #???
println("Опорная точка: x=", reference.x, ", y=", reference.y, ", γ=", gamma0) #???
reference_data = _gamma_sample(reference.x, reference.y; N=N, M=M, primes=true) #???
if reference_data !== nothing #???
    println("Исходные точки x₁…x₆: ", reference_data.points) #???
    println("Их образы T(x₁)…T(x₆): ", reference_data.normalization.images) #???
    t4prime = reference_data.normalization.transform(reference_data.x4prime) #???
    t5prime = reference_data.normalization.transform(reference_data.x5prime) #???
    println("Образы дополнительных границ: T(x₄′)=", t4prime, #???
            ", T(x₅′)=", t5prime) #???
    t4, t5 = reference_data.normalization.images[4:5] #???
    first_visible = isfinite(t4prime) && -1 < t5 < 1 && #???
        max(-1.0, min(t4, t4prime)) < min(1.0, t5, max(t4, t4prime)) #???
    second_visible = isfinite(t5prime) && -1 < t4 < 1 && #???
        max(-1.0, t4, min(t5, t5prime)) < min(1.0, max(t5, t5prime)) #???
    println("В опорной точке есть допустимый интервал: первая семья=", #???
            first_visible, ", вторая семья=", second_visible) #???
end #???

# !!! plot_gamma_slice(gamma0; N=N, M=M,
gamma0 = 100
plot_gamma_slice(gamma0; N=N, #??? M=nothing объединяет все M=1:floor(N/2).
    x_bounds=(1e-6, 110.0), y_bounds=(1e-6, 110.0), #???
    nx=1000, ny=1000, pixels=2000, nt=15, #???
    graph_name="gamma_11_1.0001") #???
