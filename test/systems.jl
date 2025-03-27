##
##      Systems to test the library.
##
##      1. Bernoulli Shifted Generalized (Beta-X)
function beta_x(beta::Float64, len::Int; transient = 10 * len)
    result = zeros(Float64, len)

    before = rand(Float64, 1)[1]
    for i = 1:(transient + len)
        after = before * beta       # beta * x_n
        while (after > 1.0)         # mod 1
            after -= 1.0
        end

        before = after
        if (i > transient)
            result[i - transient] = before
        end
    end

    return result
end
##
##      2. Logistic Map
function logistic_map(r::Float64, len::Int; transient = 10 * len)
    result = zeros(Float64, len)

    before = rand(Float64, 1)[1]
    for i = 1:(transient + len)
        before = r * before * (1 - before)
        if (i > transient)
            result[i - transient] = before
        end
    end

    return result
end
##
##      3. Lorenz System
function lorenz!(du, u, p, dt)
    x, y, z = u
    sigma, rho, beta = p

    du[1] = sigma * (y - x)
    du[2] = x * (rho - z) - y
    du[3] = x * y - beta * z
end
##
##      4. Rössler System
function rossler!(du, u, p, dt)
    x, y, z = u
    a, b, c = p

    du[1] = - y - z
    du[2] = x + a * y
    du[3] = b + z * (x - c)
end