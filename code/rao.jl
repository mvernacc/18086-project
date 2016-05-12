# Rao nozzle geometry parabolic approximation
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

export nozzle_parameters, nozzle_contour, nozzle_contour_derivative

function nozzle_parameters(r_c, r_t, θ_1, θ_2, r_e)
    # Suggested r_1 and r_2 values from Huzel and Huang are
    # r_1 = 1.5 * r_t
    # r_2 = 0.4 * r_t
    # I increase these a bit; the sharp curvature gives the solver
    # problems.
    r_1 = 3.0 * r_t
    r_2 = 1.0 * r_t
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


function nozzle_contour(x, parameters)
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
    else
        # Divergent parabolic section
        return d_1 * x^2 + d_2 * x + d_3
    end
end

function nozzle_contour_derivative(x, parameters)
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
    else
        # Divergent parabolic section
        return d_1 * x + d_2
    end
end

