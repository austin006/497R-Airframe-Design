using VortexLattice

spanvalue = 8
chordvalue = 1

# geometry (right half of the wing)
xle = [0.0, 0.0]
yle = [0.0, spanvalue/2]
zle = [0.0, 0.0]
chord = [chordvalue, chordvalue]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section

# discretization parameters
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
cref = chordvalue
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


Qinf = (1/2) * (1.204) * (Vinf)^2

Di = CD * Qinf * Sref
L = CL * Qinf * Sref

@show CL, CD, (CL/CD), Qinf, Di, L

#We can also generate files to visualize the results in Paraview using the function write_vtk.
properties = get_surface_properties(system)
write_vtk("mirrored-planar-wing", surfaces, properties; symmetric)