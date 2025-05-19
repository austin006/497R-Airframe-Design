using VortexLattice
using SNOW
using Plots
using FLOWMath: akima, Akima, derivative, gradient

# Optimization for an elliptic wing

#fitting wing geometry to a spline
function splinegeometry(x, p, semispan)
    # Fit a spline to chord values
    spliney = x
    splinex = collect(range(0, semispan, length(x)))

    cspline = Akima(splinex, spliney)
    spanpt = range(0, semispan, p)

    chords = cspline.(spanpt)
    xle = (cspline(0.0) .- cspline.(spanpt))./2
    yle = collect(spanpt)

    return xle, yle, chords
end

function getGrid(chords, nc = 5; semispan = 4.0)                                
    ns = length(chords)-1                                                       
    xyz = zeros(3, nc+1, ns+1)                                                  
    yvec = LinRange(0.0, semispan, ns+1)                                        
                                                                                
    for i in 1:nc+1                                                             
        xyz[2, i, :] .= yvec                                                    
    end                                                                         
    for i in 1:ns+1                                                             
        xyz[1, :, i] .= collect(LinRange(0.0, chords[i], nc+1))                 
    end                                                                         
                                                                                
    # Make the quarter-chord line zero sweep                                    
    for i in 1:ns+1                                                             
        xyz[1, :, i] .-= chords[i]*0.25                                         
    end                                                                         
                                                                                
    return xyz                                                                  
                                                                                
end

# Aerodynamics with VortexLattice.jl
function aerodynamics(chords, spanvalue)

    points = 50
    xle, yle, chords = splinegeometry(chords, points, spanvalue/2)
    zle = zeros(points)
    theta = zeros(points)
    phi = zeros(points)
    fc = fill((xc) -> 0, points) # camberline function for each section

    # discretization parameters
    ns = points - 1
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()

    # reference parameters
    cref = (sum(chords))/(length(chords))
    bref = spanvalue
    Sref = cref * bref
    rref = [0.0, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # construct geometry with mirror image
     grid, surface = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
         fc=fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

    # xyz = getGrid(chords, length(chords); semispan = spanvalue/2)
    # grid, surface = grid_to_surface_panels(xyz)

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
function func!(g, x) # x = [chord1, c2, c3, c4]
    spanvalue = 8
    Di, L, CL, bref, Sref, CD, cref = aerodynamics(x, spanvalue)

    # Constraint: Lift is constant
    g[1] = L - 1.7
    g[2:end] = diff(x)

    return Di
end

ns = 3 # number of chord distribution values
x0 = ones(ns) # starting point
lx = 0.0*ones(ns) # lower bounds on x
ux = 10.0*ones(ns)  # upper bounds on x

# constraints for limiting chord values
ng = 1+(ns-1) # number of constraints
lg = -Inf * ones(ng) # lower bounds on g
lg[1] = 0.0
ug = zeros(ng) # upper bounds on g

# constraints for not limiting chord values
# ng = 1 # number of constraints
# lg = [0.0] # lower bounds on g
# ug = [0.0] # upper bounds on g

ip_options = Dict("max_iter" => 300, "tol" => 1e-8)
solver = IPOPT(ip_options)
options = Options(;solver)
# options = Options(solver=IPOPT())  # choosing IPOPT solver
# options = Options(derivatives=ForwardAD()) # forward-mode algorithmic differentiation

xopt, fopt, info = minimize(func!, x0, ng, lx, ux, lg, ug, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)