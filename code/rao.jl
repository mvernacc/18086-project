# Rao nozzle geometry parabolic approximation
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

using PyPlot
function find_parameters(r_c, r_t, θ_1, θ_2, r_e)
    # Suggested r_1 and r_2 values from Huzel and Huang
    r_1 = 1.5 * r_t
    r_2 = 0.4 * r_t
    # Convergent straight section parameters
    c_1 = -tan(θ_1)
    c_2 = r_c
    x_c = - 1 / c_1 * (r_c - r_t - r_1*(1-cos(θ_1)))
    # Throat location
    x_t = x_c + r_1 * sin(θ_1)
    # Supersonic circle-parabola joint location
    x_2 = x_t + r_2 * sin(θ_2)
    # Divergent section length
    # Equivalent to length for a 15-degree conical nozzle
    L_e = 1 / tan(deg2rad(15)) * (r_e - r_t)
    # Nozzle length
    x_e = x_t + L_e
    # Supersonic parabola
    A = [(x_e)^2 x_e 1;
         2*x_2 1 0;
         x_2^2 x_2 1]
    b = [r_e, tan(θ_2), r_t + (1 - cos(θ_2)) * r_2]
    d = A\b
    return (r_c, r_t, r_1, r_2, c_1, c_2, x_c, x_t, x_2, x_e, d[1], d[2], d[3])
end


function contour(x, parameters)
    r_c, r_t, r_1, r_2, c_1, c_2, x_c, x_t, x_2, x_e, d_1, d_2, d_3 = parameters
    if x < x_c
        # Convergent line section
        return c_1 * x + c_2
    elseif x < x_t
        # Convergent circle section
        return r_1 + r_t - (r_1^2 - (x - x_t)^2)^0.5
    elseif x < x_2
        # Divergent circle section
        return r_2 + r_t - (r_2^2 - (x - x_t)^2)^0.5
    elseif x < x_e
        # Divergent parabolic section
        return d_1 * x^2 + d_2 * x + d_3
    else
        return 0
    end
end

function contour_derivative(x, parameters)
    r_c, r_t, r_1, r_2, c_1, c_2, x_c, x_t, x_2, x_e, d_1, d_2, d_3 = parameters
    if x < x_c
        # Convergent line section
        return c_1
    elseif x < x_t
        # Convergent circle section
        return (x - x_t) / (r_1^2 - (x - x_t)^2)^0.5
    elseif x < x_2
        # Divergent circle section
        return (x - x_t) / (r_2^2 - (x - x_t)^2)^0.5
    elseif x < x_e
        # Divergent parabolic section
        return d_1 * x + d_2
    else
        return 0
    end
end

p = find_parameters(4, 1, deg2rad(30), deg2rad(25), 8)
dump(p)
println(p)
x = collect(0:1e-2:p[10])
y = [contour(x1, p) for x1 in x]
dydx = [contour_derivative(x1, p) for x1 in x]
dump(x)
dump(y)
subplot(2,1,1)
plot(x, y)
subplot(2,1,2)
plot(x, dydx, marker="x")
show()

