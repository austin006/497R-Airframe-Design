using VortexLattice
using Plots
using FLOWMath: akima, Akima, derivative, gradient

spanvalue = 8
#[semi-span, c1, c2, c3, etc]
#x = [spanvalue/2, 0.8803291996790129, 1.214471372743616, 0.06573723038099441]
#x = [spanvalue/2, 0.590942576756044, 1.4853407447066955, 0.8446925099785731, 1.078459164092597, 0.019775658941325845]
#the answer:
x = [spanvalue/2, 1.0, 0.968, 0.866, 0.661, 0]
# cspline = rootchord * sqrt( 1 - ( (x^2)/(b/2)^2 ) )

x = [spanvalue/2, 0.869334121079466, 0.9807533574315391, 0.7239840724839008, 0.7942531028174915, 0.03878982273662209]

points = 50

y = x[2:end]  # span, c1, c2, etc
splinex = collect(range(0, x[1], length(x)-1))

cspline = Akima(splinex, y)

spanpt = range(0, x[1], points)
chord = cspline.(spanpt)
xle = (cspline(0.0) .- cspline.(spanpt))/2
yle = collect(spanpt)

# scatter(splinex, y)
# # scatter(yle, chord)
# plot!(spanpt, cspline.(spanpt), aspect_ratio=:equal)

# geometry (right half of the wing)
zle = zeros(points)
theta = (2.0*pi/180)*ones(points)
phi = zeros(points)
fc = fill((xc) -> 0, points) # camberline function for each section

# discretization parameters
ns = points - 1
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
cref = (sum(x) - x[1])/(length(x)-1)
bref = x[1]*2
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


Qinf = (1/2) * (1.204) * (Vinf)^2

Di = CD * Qinf * Sref
L = CL * Qinf * Sref

# println("b = $bref, sref = $Sref, AR = $(bref*bref/Sref)")
# println("CL is: $CL and CD is: $CD")
# println("CL/CD is: $(CL/CD)")
# println("The Qinf is: $Qinf")
println("Di is: $Di and the L is: $L")

#@show CL, CD, (CL/CD), Qinf, Di, L

#We can also generate files to visualize the results in Paraview using the function write_vtk.
properties = get_surface_properties(system)
write_vtk("mirrored-planar-wing", surfaces, properties; symmetric)