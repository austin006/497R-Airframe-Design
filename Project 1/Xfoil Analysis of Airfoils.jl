using Xfoil, Plots, Printf
pyplot()

function plotCl(n, a)
    # set operating conditions
    alpha = -9:1:14 # range of angle of attacks, in degrees
    re = n # Reynolds number

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)
    
    # read airfoil coordinates from a file
    x, y = open(a, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end

    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)
    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane()

    # determine airfoil coefficients across a range of angle of attacks
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
    end

    plot!(alpha, c_l, label="", xlabel="Angle of Attack (degrees)", ylabel="Lift Coefficient", show=true)

end
function plotCd(n, a)
    # set operating conditions
    alpha = -9:1:14 # range of angle of attacks, in degrees
    re = n # Reynolds number

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)
    
    # read airfoil coordinates from a file
    x, y = open(a, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end

    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)
    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane()

    # determine airfoil coefficients across a range of angle of attacks
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
    end

    plot!(alpha, c_d, label="", xlabel="Angle of Attack (degrees)", ylabel="Drag Coefficient",
    overwrite_figure=true, show=true)

end
function plotCm(n, a)
    # set operating conditions
    alpha = -9:1:14 # range of angle of attacks, in degrees
    re = n # Reynolds number

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)
    
    # read airfoil coordinates from a file
    x, y = open(a, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end

    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)
    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane()

    # determine airfoil coefficients across a range of angle of attacks
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
    end

    plot!(alpha, c_m, label="", xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient",
    overwrite_figure=false, show=true)

end
function plotClCdRatio(n, a)
    # set operating conditions
    alpha = -9:1:14 # range of angle of attacks, in degrees
    re = n # Reynolds number

    # initialize outputs
    n_a = length(alpha)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)
    
    # read airfoil coordinates from a file
    x, y = open(a, "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        end
        x, y
    end

    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)
    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane()

    # determine airfoil coefficients across a range of angle of attacks
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; iter=100, reinit=true)
    end

    plot!(c_d, c_l, label="", xlabel="Drag Coefficient", ylabel="Lift Coefficient", show=true)

end

ReN = 1e5
plot()
naca0008 = ".\\Airfoils.dat\\naca0008.dat"
naca0012 = ".\\Airfoils.dat\\naca0012.dat"
naca0021 = ".\\Airfoils.dat\\naca0021.dat"
naca2412 = ".\\Airfoils.dat\\naca2412.dat"
naca4412 = ".\\Airfoils.dat\\naca4412.dat"

plotCm(ReN, naca0008)
plotCm(ReN, naca0012)
plotCm(ReN, naca0012)
plotCm(ReN, naca2412)
plotCm(ReN, naca4412)






