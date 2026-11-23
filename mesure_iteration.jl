const tol = 1e-6
const max_iter = 100

function w(x, y, z)
    return exp(-x^2-y^2-z^2)
end

function f(x)
    return x/(x-1)
end 

function f_der(x)
    return -1/(x-1)^2
end

function iter_arg_w(x, y, z; max_iter=100)
    for i in 1:max_iter
        x_new = f(x)/f(y)
        y_new = f(x)/f(z)
        z_new = f(x)
        
        x = x_new
        y = y_new
        z = z_new
    end
    return x, y, z
end

function iter_w(x, y, z, w; max_iter=100)
    for i in 1:max_iter
        x_new = f(x) / f(y)
        y_new = f(x) / f(z)
        z_new = f(x)

        old_w = w

        w_new = let old_w = old_w
            (x, y, z) -> old_w(
                f(x) / f(y),
                f(x) / f(z),
                f(x)
            ) * f_der(x) * f(x)^2
        end

        x = x_new
        y = y_new
        z = z_new
        w = w_new
    end

    return x, y, z, w
end

x = 1.23
y = 2.13
z = 3.41
#w = (x, y, z) -> 1/(1+x^2+y^2+z^2)
mx_iter_1 = 20
mx_iter_2 = 21
mx_iter_3 = 22

w_val_1 = iter_w(x, y, z, w; max_iter=mx_iter_1)[4](x, y, z)
w_val_2 = iter_w(x, y, z, w; max_iter=mx_iter_2)[4](x, y, z)
w_val_3 = iter_w(x, y, z, w; max_iter=mx_iter_3)[4](x, y, z)

println("x = ", x, ", y = ", y, ", z = ", z, ", max_iter1 = ", mx_iter_1, ": ", iter_arg_w(x, y, z; max_iter=100), " w_val_1 = ", w_val_1)
println("x = ", x, ", y = ", y, ", z = ", z, ", max_iter2 = ", mx_iter_2, ": ", iter_arg_w(x, y, z; max_iter=101), " w_val_2 = ", w_val_2)
println("x = ", x, ", y = ", y, ", z = ", z, ", max_iter3 = ", mx_iter_3, ": ", iter_arg_w(x, y, z; max_iter=102), " w_val_3 = ", w_val_3)
println()
