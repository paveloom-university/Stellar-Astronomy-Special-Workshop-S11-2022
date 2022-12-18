# This script performs weighted orthogonal regression

"Parse the string, taking more arguments if it's quoted"
function parse_string(i)::String
    # Start from the first argument after the flag
    j = i
    # If the arguments starts with an apostrophe,
    s = if startswith(ARGS[i], "'")
        # Find the index of the argument
        # which ends with an apostrophe
        while !endswith(ARGS[j], "'")
            j += 1
        end
        # Join the arguments in one string
        # and remove the apostrophes
        chop(join(ARGS[i:j], ' '), head=1, tail=1)
    else
        # Return the next argument
        ARGS[i]
    end
    return s
end

# Define default values for optional arguments
POSTFIX = ""
FORCE = false

# Parse the options
for i in eachindex(ARGS)
    # A postfix for the names of output files
    if ARGS[i] == "--postfix"
        try
            global POSTFIX = " ($(parse_string(i+1)))"
        catch
            println("Couldn't parse the value of the `--postfix` argument.")
            exit(1)
        end
    end
end

# Prepare color codes
RESET = "\e[0m"
GREEN = "\e[32m"
YELLOW = "\e[33m"

# Check for required arguments
if "--help" in ARGS
    println("""
        $(YELLOW)USAGE:$(RESET)
        { julia --project=. | ./julia.bash } scripts/compute.jl [--postfix <POSTFIX>]

        $(YELLOW)OPTIONS:$(RESET)
            $(GREEN)--postfix <POSTFIX>$(RESET)    A postfix for the names of output files"""
    )
    exit(1)
end

"Padding in the output"
pad = " "^4

"Floating point type used across the script"
F = Float64

"Integer type used across the script"
I = Int64

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Optim
using Plots
using Roots

# Use the PGFPlotsX backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(
    fontfamily="Computer Modern",
    dpi=300,
    size=(400, 400),
    legend=nothing,
)

# Define the data set
X = [0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
ω_X = [1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1.0]
σ_X = @. 1 / √(ω_X)
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]
ω_Y = [1, 1.8, 4, 8, 20, 20, 70, 70, 100, 500]
σ_Y = @. 1 / √(ω_Y)

# Compute the size of the data set
N = length(X)

"Compute a weighted squared difference"
function diff(x::F, a::F, b::F, i::I)
    return ω_Y[i] * (a + b * x - Y[i])^2 + ω_X[i] * (x - X[i])^2
end

"Compute a sum of weighted squared differences"
function wsum(x::Vector{F}, a::F, b::F)
    return sum(i -> diff(x[i], a, b, i), 1:N)
end

"""
Find the best parameters (and their average errors) of
the linear model `y = a + bx` using the D. York's method.
Return also the arrays of deviations for projections

The arguments here are the approximate value of `b` and
the required precision.
"""
function york(; b₀::F, ε::F, linear::Bool=false)::Tuple{F,F,F,F,Vector{F},Vector{F}}
    # Start with the provided approximate value of `b`
    b = b₀
    b_prev = b₀ + 2ε
    # Keep several global variables for
    # getting the results out of the loop
    a = σ_a = σ_b = 0
    x = y = similar(X)
    # Until the required precision is reached
    while abs(b_prev - b) > ε
        # Save the current `b` as the previous one
        b_prev = b
        # Compute the weights of each pair
        W = @. ω_X * ω_Y / (ω_X + b^2 * ω_Y)
        # Compute the sum of the weights
        W_sum = sum(W)
        # Compute the gravity center
        G_X = sum(W .* X) / W_sum
        G_Y = sum(W .* Y) / W_sum
        # Compute the deviations from the gravity center
        U = X .- G_X
        V = Y .- G_Y
        # If the linear approximation is requested, do that.
        # Otherwise, use the cubic equation for the approximation
        if linear
            # Compute the next approximation of the `b` parameter
            b = sum(@. W^2 * V * (U / ω_Y + b * V / ω_X)) /
                sum(@. W^2 * U * (U / ω_Y + b * V / ω_X))
        else
            # Compute the coefficients of the cubic equation
            d = 3 * sum(@. W^2 * U^2 / ω_X)
            α = 2 * sum(@. W^2 * U * V / ω_X) / d
            β = (sum(@. W^2 * V^2 / ω_X) - sum(@. W * U^2)) / d
            γ = -sum(@. W * U * V) * 3 / d
            # Compute the next approximation of the `b` parameter
            φ = acos((α^3 - 1.5 * α * β + 0.5 * γ) / (α^2 - β)^(1.5))
            b = α + 2 * √(α^2 - β) * cos((φ + 4π) / 3)
        end
        # Compute the next approximation of the `a` parameter
        a = G_Y - b * G_X
        # Compute the average errors
        σ²_b = 1 / (N - 2) * sum(@. W * (b * U - V)^2) / sum(@. W * U^2)
        σ²_a = σ²_b * sum(@. W * X^2) / W_sum
        σ_a = √(σ²_a)
        σ_b = √(σ²_b)
        # Compute the coordinates of the projections
        x = @. X + -b * W / ω_X * (a + b * X - Y)
        y = @. Y + W / ω_Y * (a + b * X - Y)
    end
    return (a, b, σ_a, σ_b, x, y)
end

println(pad, "> Fitting the linear model using the D. York's method...")

# Fit the linear model
a_l, b_l, σ_a_l, σ_b_l, x_l, y_l = york(b₀=-0.5, ε=1e-8, linear=true)
a_c, b_c, σ_a_c, σ_b_c, x_c, y_c = york(b₀=-0.5, ε=1e-8)

# Compute the sums of the weighted squared differences
∑_l = wsum(x_l, a_l, b_l)
∑_c = wsum(x_c, a_c, b_c)

println('\n', pad, pad, "Linear:")
println(pad, pad, pad, "a: ", a_l, " ± ", σ_a_l)
println(pad, pad, pad, "b: ", b_l, " ± ", σ_b_l)
println(pad, pad, pad, "∑: ", ∑_l, '\n')
println(pad, pad, "Cubic:")
println(pad, pad, pad, "a: ", a_c, " ± ", σ_a_c)
println(pad, pad, pad, "b: ", b_c, " ± ", σ_b_c)
println(pad, pad, pad, "∑: ", ∑_c, '\n')

println(pad, "> Plotting the results of the D. York's method...")

# Plot the data
scatter(
    X,
    Y,
    xerror=σ_X,
    yerror=σ_Y,
    markershape=:xcross,
    aspect_ratio=:equal,
    linewidth=0.1,
    markersize=1,
    markerstrokewidth=0.1,
    xlabel=L"x",
    ylabel=L"y",
)

# Plot the regression line using the D. York's method
plot!(x -> a_c + b_c * x, linewidth=0.1)

# Plot the projections
scatter!(
    x_c,
    y_c,
    markershape=:xcross,
    markersize=1,
    markerstrokewidth=0.1,
)
for (X, Y, x, y) in zip(X, Y, x_c, y_c)
    plot!([X, x], [Y, y], linewidth=0.2)
end

# Save the figure
savefig(joinpath(PLOTS_DIR, "york.pdf"))

# Define the initial values for the numerical optimization
θ₀ = [6, -0.5]

"""
The target function for numerical optimization

It can be shown that this function can be treated as a
non-constant part of the log likelihood function
"""
function target(θ::Vector{F})::F
    # The input parameters define the line
    a, b = θ
    # We go through each point of the input data and
    # find a point on this line which is the closest
    return 0.5 * sum(
        i -> optimize(θ -> diff(θ[1], a, b, i), [X[i]]).minimum,
        1:N
    )
end

"Find the best parameters of the linear model `y = a + bx`
by numerically minimizing the target function"
function minimize()::Tuple{F,F,F}
    res = optimize(target, θ₀)
    return res.minimizer..., res.minimum
end

"Compute the minimum of the target function
while one of the parameters is frozen"
function minimize_frozen(
    idx::I,
    value::F,
    θ₀::Vector{F},
)::Float64
    # Exclude the frozen parameter from the active ones
    θ₀ = [θ₀[1:idx-1]; θ₀[idx+1:end]]
    # Recreate the target function, handling a frozen parameter
    function target_inner(θ::Vector{F})
        return target([θ[1:idx-1]; value; θ[idx:end]])
    end
    # Optimize the new target function
    res = optimize(target_inner, θ₀)
    # Return the minimum
    return res.minimum
end

"Compute the confidence intervals for the parameter"
function intervals(idx::I, value::F, L₀::F)::Tuple{F,F,F,F}
    # Create an alias function for the frozen parameter
    minimize_with(p) = minimize_frozen(idx, p, θ₀)
    # Find the first root
    r = find_zero((p) -> minimize_with(p) - L₀ - 0.5, value)
    # Find the second one nearby; compute the confidence intervals
    r_left = r_right = σ_left = σ_right = 0
    # If the root is to the left from the value
    if value - r > 0
        r_left = r
        σ_left = abs(value - r_left)
        r_right = find_zero((p) -> minimize_with(p) - L₀ - 0.5, value + σ_left)
        σ_right = abs(value - r_right)
    else
        r_right = r
        σ_right = abs(value - r_right)
        r_left = find_zero((p) -> minimize_with(p) - L₀ - 0.5, value - σ_right)
        σ_left = abs(value - r_left)
    end
    return r_left, r_right, σ_left, σ_right
end

"Plot the profile of the parameter"
function plot_profile(
    idx::I,
    value::F,
    name::Symbol,
    r_left::F,
    r_right::F,
    σ_left::F,
    σ_right::F,
    L₀::F,
)
    # Create an alias function for the frozen parameter
    minimize_with(p) = minimize_frozen(idx, p, θ₀)
    # Plot the profile of the parameter
    p = plot(
        minimize_with,
        r_left - σ_left,
        r_right + σ_right;
        xlabel=latexstring("$(name)"),
        ylabel=latexstring("L_$(name)($(name))")
    )
    # Add vertical lines to the plot
    plot!(p, [r_left, r_left], [L₀ - 0.2, L₀ + 0.5]; linestyle=:dash)
    plot!(p, [value, value], [L₀ - 0.2, L₀]; linestyle=:dash)
    plot!(p, [r_right, r_right], [L₀ - 0.2, L₀ + 0.5]; linestyle=:dash)
    # Add points to the plot
    scatter!(p, [r_left], [L₀ + 0.5,])
    scatter!(p, [value], [L₀,])
    scatter!(p, [r_right], [L₀ + 0.5,])
    # Add horizontal lines to the plot
    hline!(p, [L₀ + 0.5,]; linestyle=:dash)
    hline!(p, [L₀,]; linestyle=:dash)
    # Add annotations to the plot
    annotate!(
        p,
        [
            (v + 0.2 * (value - r_left), L₀ - 0.1, text("$(round(v; digits=2))", 8))
            for v in [r_left, value, r_right]
        ]
    )
    annotate!(
        p,
        [
            (value - 1.85 * (value - r_left), L₀ + 0.06, text(L"L_0", 8)),
            (value - 1.65 * (value - r_left), L₀ + 0.5 + 0.06, text(L"L_0 + 1/2", 8)),
        ]
    )

    # Save the figure
    savefig(p, "$(PLOTS_DIR)/$(name).pdf")
end

println(pad, "> Fitting the linear model by numerical optimization...")

# Fit the linear model
a_n, b_n, ∑_n = minimize()

println('\n', pad, pad, "a: ", a_n)
println(pad, pad, "b: ", b_n)
println(pad, pad, "∑: ", ∑_n, '\n')

println(pad, "> Computing the confidence intervals...")

# Compute the confidence intervals
a_left, a_right, a_σ_left, a_σ_right = intervals(1, a_n, ∑_n)
b_left, b_right, b_σ_left, b_σ_right = intervals(2, b_n, ∑_n)

println('\n', pad, pad, "a_left: ", a_left)
println(pad, pad, "a_right: ", a_right)
println(pad, pad, "σ_left: ", a_σ_left)
println(pad, pad, "σ_right: ", a_σ_right, '\n')
println(pad, pad, "b_left: ", b_left)
println(pad, pad, "b_right: ", b_right)
println(pad, pad, "σ_left: ", b_σ_left)
println(pad, pad, "σ_right: ", b_σ_right, '\n')

println(pad, "> Plotting the profiles...")

# Plot the profiles
plot_profile(1, a_n, :a, a_left, a_right, a_σ_left, a_σ_right, ∑_n)
plot_profile(2, b_n, :b, b_left, b_right, b_σ_left, b_σ_right, ∑_n)

println()
