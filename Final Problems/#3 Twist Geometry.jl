using VortexLattice
using Plots
using FLOWMath: akima, Akima, derivative, gradient

#[t1, t2, t3, etc]
x = [0.08689306092649526, 0.08517784073392665, 0.0780651844530549, 0.06617491944571376, 0.04172586118048137]
#x = [0.08613842882004043, 0.07843036402538514, 0.04525050917931647]

spanvalue = 8
chordvalue = 1
points = 50

spliney = x
splinex = collect(range(0, spanvalue/2, length(x)))

tspline = Akima(splinex, spliney)
spanpt = range(0, spanvalue/2, points)

yle = collect(spanpt)
theta = tspline.(spanpt)

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
cref = chordvalue
bref = spanvalue
Sref = cref * bref
rref = [0.0, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 0.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# construct geometry with mirror image
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc=fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=false)

# symmetry is not used in the analysis
symmetric = true

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


Qinf = (1/2) * (1.204) * (Vinf)^2

Di = CD * Qinf * Sref
L = CL * Qinf * Sref

@show Di, L
#@show CL, CD, (CL/CD), Qinf, Di, L

#lifting_line_coefficients
r, c = lifting_line_geometry([grid])
Cf, Cm = lifting_line_coefficients(system, r, c; frame=Wind())

#Plot span vs lift coefficient
Cf2 = Cf[1]
cl = Cf2[3,:]
plot(LinRange(0, spanvalue/2, length(cl)), cl, label="resulting curve", title="Elliptical Comparison", xlabel="Span (m)", ylabel="Coefficient of Lift", linewidth=2, show=true)

#Elliptical Lift distribution
Xell = LinRange(0, spanvalue/2, length(cl))
Yell =  cl[1] .* sqrt.( 1 .- ( (Xell.^2) ./ (spanvalue/2)^2 ) )
plot!(Xell, Yell, label="elliptical curve", linewidth=2)

#We can also generate files to visualize the results in Paraview using the function write_vtk.
# properties = get_surface_properties(system)
# write_vtk("mirrored-planar-wing", surfaces, properties; symmetric)
# println("Your wing has been plotted :)")