include("../../src/stiefel/func.jl")
using Plots
using QuadGK
# import res

#setprecision(512)
try
    τ = 0.6
    N = 10
    M = 6
    

    
    #=
    u_11 = 0+0*im
    u_12 = 0+0*im

    u_21 = (0.25-1*τ*im/4)

    u_22 = path_on_jacobian_new(u_11, u_12, u_21, p; tol=1e-11, step_num=100, frac_for_tang=100000)

    println("[u_11, u_12] = [", u_11, " ", u_12, "]")
    println("θ_shtiff([u_11, u_12]) = ", θ_shtiff([u_11, u_12], p; char_shtiff=[3, 5]))
    println()
    println("[u_21, u_22] = [", u_21, " ", u_22, "]")
    println("θ_shtiff([u_21, u_22]) = ", θ_shtiff([u_21, u_22], p; char_shtiff=[3, 5]))

    u_3 = project_to_image([u_21, u_22], p; shtiff_chars=[[3, 5]], act_vars=[2], tol=1e-10, maxiter=100)
    println()
    println("[u_31, u_32] = [", u_3)
    println("θ_shtiff([u_31, u_32]) = ", θ_shtiff([u_3[1], u_3[2]], p; char_shtiff=[3, 5]))
    =#
    #=
    f(t) = begin
        #println(typeof(p))
        p_copy = copy(p)
        p_copy[1, 1] = p_copy[1, 1] + im*t
        return θ_shtiff([u_1, u_2], p_copy, char_shtiff = [3, 5])
    end

    df(t) = begin
        #println(typeof(p))
        p_copy = copy(p)
        p_copy[1,1] = p_copy[1, 1] + im*t
        return im*θ_shtiff([u_1, u_2], p_copy, derivs_t=[[1, 1]], char_shtiff = [3, 5])
    end

    a = 1
    b = 2
    I, err = quadgk(df, a, b; rtol=1e-12, atol=0)
    println("Integral: ", I, " with error: ", err)
    println("True Integral: ", f(b)-f(a))

    println("[u_1, u_2] = [", u_1, " ", u_2, "]")
    u = project_to_image([u_1, u_2], p; shtiff_chars = [[3, 5]], act_vars = [2])
    println("[u_1', u_2'] = [", u[1], " ", u[2], "]")

    println(θ_shtiff(u, p; char_shtiff = [3, 5]))
    =#


    #=
    time_reim = time()
    for i=1:1000
        p_reim = RiemannMatrix(p, siegel=false)
    end
    time_reim = time() - time_reim
    println("time ", time_reim/1000)

    p_reim = RiemannMatrix(p, siegel=false)


    time_theta = time()
    for i=1:1000
        θ_faster([1, 3+im*4], p_reim)
    end
    time_theta = time() - time_theta
    println("time ", time_theta/1000)
    p_reim = RiemannMatrix(p, siegel=false)
    =#

    #println("τ ", τ)

    #=
    f(x, y) = norm(dθ_dΩ([x + im * y, 1], p_reim, 2, 2))
    xs = range(-1, 1; length=100)
    ys = range(-1, 1; length=100)
    Z = [f(x, y) for y in ys, x in xs]
    Plots.surface(xs, ys, Z;
        xlabel="x",
        ylabel="y",
        zlabel="f(x,y)",
        title="Graph of f(x,y)"
    )
    =#
catch e
    println("Error type: ", typeof(e))
    println("Error message:")
    showerror(stdout, e)
    println("\n\nStacktrace:")
    Base.show_backtrace(stdout, catch_backtrace())
    println()
end
