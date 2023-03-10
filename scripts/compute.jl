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
REPORT_DIR = joinpath(ROOT_DIR, "report")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Optim
using Plots
using Printf
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
??_X = [1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1.0]
??_X = @. 1 / ???(??_X)
Y = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5]
??_Y = [1, 1.8, 4, 8, 20, 20, 70, 70, 100, 500]
??_Y = @. 1 / ???(??_Y)

# Compute the size of the data set
N = length(X)

"Compute a weighted squared difference"
function diff(x::F, a::F, b::F, i::I)
    return ??_Y[i] * (a + b * x - Y[i])^2 + ??_X[i] * (x - X[i])^2
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
function york(; b???::F, ??::F, linear::Bool=false)::Tuple{F,F,F,F,Vector{F},Vector{F}}
    # Number of digits to round the output to
    digits = 9
    # Print the header
    if linear
        println('\n', pad, pad, "Linear:")
    else
        println('\n', pad, pad, "Cubic:")
    end
    println('\n', pad, pad, pad, "b???: ", b???, '\n')
    # Start with the provided approximate value of `b`
    b = b???
    b_prev = b??? + 2??
    # Keep several global variables for
    # getting the results out of the loop
    a = ??_a = ??_b = 0
    x = y = similar(X)
    # Store the iterations in a TeX file
    io = open(
        joinpath(
            REPORT_DIR,
            "iterations_" *
            (linear ? "linear" : "cubic") *
            ".tex"
        ),
        "w"
    )
    println(
        io,
        raw"""
        \begin{table}[h]
          \centering
          \begin{tabular}{cccccc}
            \toprule
            $ i $ &
            $ a $ &
            $ ??_a $ &
            $ b $ &
            $ ??_b $ &
            $ \sum $ \\
            \midrule"""
    )
    # Until the required precision is reached
    i = 0
    while abs(b_prev - b) > ??
        i += 1
        # Save the current `b` as the previous one
        b_prev = b
        # Compute the weights of each pair
        W = @. ??_X * ??_Y / (??_X + b^2 * ??_Y)
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
            b = sum(@. W^2 * V * (U / ??_Y + b * V / ??_X)) /
                sum(@. W^2 * U * (U / ??_Y + b * V / ??_X))
        else
            # Compute the coefficients of the cubic equation
            d = 3 * sum(@. W^2 * U^2 / ??_X)
            ?? = 2 * sum(@. W^2 * U * V / ??_X) / d
            ?? = (sum(@. W^2 * V^2 / ??_X) - sum(@. W * U^2)) / d
            ?? = -sum(@. W * U * V) * 3 / d
            # Compute the next approximation of the `b` parameter
            ?? = acos((??^3 - 1.5 * ?? * ?? + 0.5 * ??) / (??^2 - ??)^(1.5))
            b = ?? + 2 * ???(??^2 - ??) * cos((?? + 4??) / 3)
        end
        # Compute the next approximation of the `a` parameter
        a = G_Y - b * G_X
        # Compute the average errors
        ????_b = 1 / (N - 2) * sum(@. W * (b * U - V)^2) / sum(@. W * U^2)
        ????_a = ????_b * sum(@. W * X^2) / W_sum
        ??_a = ???(????_a)
        ??_b = ???(????_b)
        # Compute the coordinates of the projections
        x = @. X + -b * W / ??_X * (a + b * X - Y)
        y = @. Y + W / ??_Y * (a + b * X - Y)
        # Compute the sum of weighted squared differences
        ??? = wsum(x, a, b)
        # Print the results in this iteration
        println(pad, pad, pad, "i: ", i)
        println(pad, pad, pad, pad, "a: ", a, " ?? ", ??_a)
        println(pad, pad, pad, pad, "b: ", b, " ?? ", ??_b)
        println(pad, pad, pad, pad, "???: ", ???)
        Printf.format(
            io,
            Printf.Format(
                "    %d" *
                repeat(" & %.$(digits)f", 5) *
                " \\\\\n",
            ),
            i, a, ??_a, b, ??_b, ???,
        )
    end
    # Close the file stream
    println(
        io,
        raw"""
            \bottomrule
          \end{tabular}
        """,
        if linear
            "  \\caption{???????????????? ??????????????????}\n"
        else
            "  \\caption{???????????????????? ??????????????????}\n"
        end,
        raw"\end{table}",
    )
    close(io)
    return (a, b, ??_a, ??_b, x, y)
end

println(pad, "> Fitting the linear model using the D. York's method...")

# Fit the linear model
a_l, b_l, ??_a_l, ??_b_l, x_l, y_l = york(b???=-0.5, ??=1e-8, linear=true)
a_c, b_c, ??_a_c, ??_b_c, x_c, y_c = york(b???=-0.5, ??=1e-8)

println('\n', pad, "> Plotting the results of the D. York's method...")

# Plot the data
scatter(
    X,
    Y,
    xerror=??_X,
    yerror=??_Y,
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
    plot!(
        [X, x],
        [Y, y],
        linewidth=0.2,
        color=palette(:default).colors[1],
    )
end

# Save the figure
savefig(joinpath(PLOTS_DIR, "york.pdf"))

# Define the initial values for the numerical optimization
????? = [6, -0.5]

"""
The target function for numerical optimization

It can be shown that this function can be treated as a
non-constant part of the log likelihood function
"""
function target(??::Vector{F})::F
    # The input parameters define the line
    a, b = ??
    # We go through each point of the input data and
    # find a point on this line which is the closest
    return 0.5 * sum(
        i -> optimize(?? -> diff(??[1], a, b, i), [X[i]]).minimum,
        1:N
    )
end

"Find the best parameters of the linear model `y = a + bx`
by numerically minimizing the target function"
function minimize()::Tuple{F,F,F}
    res = optimize(target, ?????)
    return res.minimizer..., res.minimum
end

"Compute the minimum of the target function
while one of the parameters is frozen"
function minimize_frozen(
    idx::I,
    value::F,
    ?????::Vector{F},
)::Float64
    # Exclude the frozen parameter from the active ones
    ????? = [?????[1:idx-1]; ?????[idx+1:end]]
    # Recreate the target function, handling a frozen parameter
    function target_inner(??::Vector{F})
        return target([??[1:idx-1]; value; ??[idx:end]])
    end
    # Optimize the new target function
    res = optimize(target_inner, ?????)
    # Return the minimum
    return res.minimum
end

"Compute the confidence intervals for the parameter"
function intervals(idx::I, value::F, L???::F)::Tuple{F,F,F,F}
    # Create an alias function for the frozen parameter
    minimize_with(p) = minimize_frozen(idx, p, ?????)
    # Find the first root
    r = find_zero((p) -> minimize_with(p) - L??? - 0.5, value)
    # Find the second one nearby; compute the confidence intervals
    r_left = r_right = ??_left = ??_right = 0
    # If the root is to the left from the value
    if value - r > 0
        r_left = r
        ??_left = abs(value - r_left)
        r_right = find_zero((p) -> minimize_with(p) - L??? - 0.5, value + ??_left)
        ??_right = abs(value - r_right)
    else
        r_right = r
        ??_right = abs(value - r_right)
        r_left = find_zero((p) -> minimize_with(p) - L??? - 0.5, value - ??_right)
        ??_left = abs(value - r_left)
    end
    return r_left, r_right, ??_left, ??_right
end

"Plot the profile of the parameter"
function plot_profile(
    idx::I,
    value::F,
    name::Symbol,
    r_left::F,
    r_right::F,
    ??_left::F,
    ??_right::F,
    L???::F,
)
    # Create an alias function for the frozen parameter
    minimize_with(p) = minimize_frozen(idx, p, ?????)
    # Plot the profile of the parameter
    p = plot(
        minimize_with,
        r_left - ??_left,
        r_right + ??_right;
        xlabel=latexstring("$(name)"),
        ylabel=latexstring("L_$(name)($(name))")
    )
    # Add vertical lines to the plot
    plot!(p, [r_left, r_left], [L??? - 0.2, L??? + 0.5]; linestyle=:dash)
    plot!(p, [value, value], [L??? - 0.2, L???]; linestyle=:dash)
    plot!(p, [r_right, r_right], [L??? - 0.2, L??? + 0.5]; linestyle=:dash)
    # Add points to the plot
    scatter!(p, [r_left], [L??? + 0.5,])
    scatter!(p, [value], [L???,])
    scatter!(p, [r_right], [L??? + 0.5,])
    # Add horizontal lines to the plot
    hline!(p, [L??? + 0.5,]; linestyle=:dash)
    hline!(p, [L???,]; linestyle=:dash)
    # Add annotations to the plot
    annotate!(
        p,
        [
            (v + 0.2 * (value - r_left), L??? - 0.1, text("$(round(v; digits=2))", 8))
            for v in [r_left, value, r_right]
        ]
    )
    annotate!(
        p,
        [
            (value - 1.85 * (value - r_left), L??? + 0.06, text(L"L_0", 8)),
            (value - 1.65 * (value - r_left), L??? + 0.5 + 0.06, text(L"L_0 + 1/2", 8)),
        ]
    )

    # Save the figure
    savefig(p, "$(PLOTS_DIR)/$(name).pdf")
end

println(pad, "> Fitting the linear model by numerical optimization...")

# Fit the linear model
a_n, b_n, ???_n = minimize()

println('\n', pad, pad, "a: ", a_n)
println(pad, pad, "b: ", b_n)
println(pad, pad, "???: ", ???_n, '\n')

println(pad, "> Computing the confidence intervals...")

# Compute the confidence intervals
a_left, a_right, a_??_left, a_??_right = intervals(1, a_n, ???_n)
b_left, b_right, b_??_left, b_??_right = intervals(2, b_n, ???_n)

println('\n', pad, pad, "a_left: ", a_left)
println(pad, pad, "a_right: ", a_right)
println(pad, pad, "??_left: ", a_??_left)
println(pad, pad, "??_right: ", a_??_right, '\n')
println(pad, pad, "b_left: ", b_left)
println(pad, pad, "b_right: ", b_right)
println(pad, pad, "??_left: ", b_??_left)
println(pad, pad, "??_right: ", b_??_right, '\n')

# Print the results to a TeX file
open(joinpath(REPORT_DIR, "numerical.tex"), "w") do io
    # Number of digits to round the output to
    digits = 6
    println(
        io,
        raw"""
        \begin{table}[h]
          \centering
          \begin{tabular}{ccccccc}
            \toprule
            $ a $ &
            $ ??_a^+ $ &
            $ ??_a^- $ &
            $ b $ &
            $ ??_b^+ $ &
            $ ??_b^- $ &
            $ \sum $ \\
            \midrule""",
    )
    Printf.format(
        io,
        Printf.Format("    %.$(digits)f" * repeat(" & %.$(digits)f", 6) * " \\\\\n"),
        a_n, a_??_right, a_??_left, b_n, b_??_right, b_??_left, ???_n
    )
    println(
        io,
        raw"""
            \bottomrule
          \end{tabular}
        \end{table}
        """
    )
end

println(pad, "> Plotting the profiles...")

# Plot the profiles
plot_profile(1, a_n, :a, a_left, a_right, a_??_left, a_??_right, ???_n)
plot_profile(2, b_n, :b, b_left, b_right, b_??_left, b_??_right, ???_n)

println()
