using VortexLattice
using SNOW
using Plots
using FLOWMath: akima, Akima, derivative, gradient

# Optimization for twist on a rectangular wing

#fitting wing geometry to a spline
function splinegeometry(x, p, semispan)
    # Fit a spline to chord values
    spliney = x
    splinex = collect(range(0, semispan, length(x)))

    tspline = Akima(splinex, spliney)
    spanpt = range(0, semispan, p)

    yle = collect(spanpt)
    theta = tspline.(spanpt) #(2.0*pi/180)*ones(points)

    return yle, theta
end

# Aerodynamics with VortexLattice.jl
function aerodynamics(x, spanvalue, chordvalue)

    points = 50
    yle, theta = splinegeometry(x, points, spanvalue/2)

    chord = chordvalue*ones(points)
    xle = zeros(points)
    zle = zeros(points)
    phi = zeros(points)
    fc = fill((xc) -> 0, points) # camberline function for each section

    # discretization parameters
    ns = points - 1
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()

    # reference parameters
    cref = (sum(chord))/(length(chord))
    bref = spanvalue
    Sref = cref * bref
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # construct geometry with mirror image
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        fc=fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

    # symmetry is not used in the analysis
    symmetric = false

    # create vector containing all surfaces
    surfaces = [surface]

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    # perform far-field analysis
    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    #Stability Derivatives
    dCF, dCM = stability_derivatives(system)

    CDa, CYa, CLa = dCF.alpha
    Cla, Cma, Cna = dCM.alpha
    CDb, CYb, CLb = dCF.beta
    Clb, Cmb, Cnb = dCM.beta
    CDp, CYp, CLp = dCF.p
    Clp, Cmp, Cnp = dCM.p
    CDq, CYq, CLq = dCF.q
    Clq, Cmq, Cnq = dCM.q
    CDr, CYr, CLr = dCF.r
    Clr, Cmr, Cnr = dCM.r

    properties = get_surface_properties(system)

    Qinf = (1/2) * (1.204) * (Vinf)^2

    Di = CD * Qinf * Sref
    L = CL * Qinf * Sref

    return Di, L, CL, bref, Sref, CD, cref
end

# Objective function
function func!(g, x) # x = [twist1, tw2, tw3]
    spanvalue = 8
    chordvalue = 1
    Di, L, CL, bref, Sref, CD, cref = aerodynamics(x, spanvalue, chordvalue)

    # Constraint: Lift is constant
    g[1] = L / 0.9924523474492697 #1.18

    return Di
end

x0 = 0.5*ones(5) # starting point
lx = (0.0*pi/180)*ones(length(x0)) # lower bounds on x
ux = (2.0*pi/180)*ones(length(x0))  # upper bounds on x
ng = 1 # number of constraints
lg = [1.0] # lower bounds on g
ug = [1.0] # upper bounds on g

ip_options = Dict("max_iter" => 100, "tol" => 1e-8)
solver = IPOPT(ip_options)
options = Options(;solver)
# options = Options(solver=IPOPT())  # choosing IPOPT solver
# options = Options(derivatives=ForwardAD()) # forward-mode algorithmic differentiation

xopt, fopt, info = minimize(func!, x0, ng, lx, ux, lg, ug, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)