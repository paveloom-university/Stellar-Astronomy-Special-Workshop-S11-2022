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
I = UInt64

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots

# Use the PGFPlotsX backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(
    fontfamily="Computer Modern",
    dpi=300,
    size=(400, 400),
    legend=nothing,
    linewidth=0.1,
    markersize=2,
    markerstrokewidth=0.1,
    xlabel=L"x",
    ylabel=L"y",
    draw_arrow=true,
    aspect_ratio=:equal,
)

# Define the data set
X = [0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4]
ω_X = [1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1.0]
σ_X = @. 1 / √(ω_X)
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]
ω_Y = [1, 1.8, 4, 8, 20, 20, 70, 70, 100, 500]
σ_Y = @. 1 / √(ω_Y)

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
    a = 0
    σ_a = 0
    σ_b = 0
    ΔX = similar(X)
    ΔY = similar(Y)
    # Compute the size of the data set
    N = length(X)
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
        # Compute the deviations
        ΔX = @. -b * W / ω_X * (a + b * X - Y)
        ΔY = @. W / ω_Y * (a + b * X - Y)
    end
    return (a, b, σ_a, σ_b, ΔX, ΔY)
end

println(pad, "> Fitting the linear model using the D. York's method...")

# Fit the linear model
a_l, b_l, σ_a_l, σ_b_l, ΔX_l, ΔY_l = york(b₀=-0.5, ε=1e-8, linear=true)
a_c, b_c, σ_a_c, σ_b_c, ΔX_c, ΔY_c = york(b₀=-0.5, ε=1e-8)

println('\n', pad, pad, "Linear:")
println(pad, pad, pad, "a: ", a_l)
println(pad, pad, pad, "b: ", b_l, '\n')
println(pad, pad, "Cubic:")
println(pad, pad, pad, "a: ", a_c)
println(pad, pad, pad, "b: ", b_c, '\n')

println(pad, "> Plotting the results of the D. York's method...")

# Plot the data
plot(X, Y, xerror=σ_X, yerror=σ_Y)

# Plot the regression line using the D. York's method
plot!(x -> a_c + b_c * x)

# Plot the projections
for (x, y, δx, δy) in zip(X, Y, ΔX_c, ΔY_c)
    plot!([x, x + δx], [y, y + δy], arrow=true)
end

# Save the figure
savefig(joinpath(PLOTS_DIR, "york.pdf"))

println()
