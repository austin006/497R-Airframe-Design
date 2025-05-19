using VortexLattice
using Plots
using FLOWMath: akima, Akima, derivative, gradient

#[c1, c2, c3, etc] the answer:
#x = [1.0, 0.968, 0.866, 0.661, 0]
# cspline = rootchord * sqrt( 1 - ( (x^2)/(b/2)^2 ) )

# Optimal solution:
#x = [0.990932018620797, 0.9901227276962783, 0.9876950391365055, 0.9836395846845459, 0.9779377415600808, 0.9705537573419005, 0.9614514416082949, 0.9505850247326799, 0.9378873506742835, 0.9232818775588912, 0.9066758650159614, 0.8879547755259892, 0.8669869631395322, 0.8435921206514606, 0.8175710615376739, 0.788653699646549, 0.7564954281158464, 0.720650363448121, 0.6805130642433932, 0.6352470697987135, 0.5837433109668767, 0.5247002199709955, 0.45702640872075895, 0.3773794195347345, 0.22969457677817823]

#x = [0.8641570361008406, 0.9560541053539827, 0.7339327253055038, 0.782333193543182, 0.033834191229674]
#x = [0.8282738550196137, 1.9998772385304404, 0.32353985307922817, 0.1020917917514297, 1.8666537472903408]
#x = [5.475105236202179, 2.4469774393077097, 7.294924245455563, 3.0065849582140887, 3.175114163671531, 7.563274546954158, 9.930158737567604, 0.8136768399450572, 3.217444515549613, 7.009232473083338, 9.999932190693734, 4.0323355086060175, 1.8748442130716778, 4.0962608090983075, 6.681614930209702, 5.312907071380102, 5.596538686642137, 5.320641471729242, 4.803053507819258, 0.19585596354371707, 8.092678045681712, 9.979556703111452, 6.195601454749648, 0.40073410624461037, 9.778526455266633]
#x = [0.9307524436528909, 0.9307523295670074, 0.7743711598753611, 0.7743710192306222, 0.3037182771899257]
#x = [0.9259052543175672, 0.9259052471943864, 0.15179738749121324]
x = [0.81853899, 0.965859761, 0.113387467967]

spanvalue = 8
points = 50

spliney = x
splinex = collect(range(0, spanvalue/2, length(x)))

cspline = Akima(splinex, spliney)
spanpt = range(0, spanvalue/2, points)

chord = cspline.(spanpt)
xle = (cspline(0.0) .- cspline.(spanpt))/2
yle = collect(spanpt)

# geometry (right half of the wing)
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
cref = (sum(x))/(length(x))
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

#@show Di, L
#@show CL, CD, (CL/CD), Qinf, Di, L

#lifting_line_coefficients
r, c = lifting_line_geometry([grid])
Cf, Cm = lifting_line_coefficients(system, r, c; frame=Wind())

#Plot span vs lift coefficient
Cf2 = Cf[1]
cl = Cf2[3,:]
plot(spanpt, cspline.(spanpt), label="resulting curve", title="Elliptical Comparison", xlabel="Span (m)", ylabel="Chord Length (m)", linewidth=2, show=true)

#Elliptical Lift distribution
Xell = LinRange(0, spanvalue/2, length(cl))
Yell =  x[1] .* sqrt.( 1 .- ( (Xell.^2) ./ (spanvalue/2)^2 ) )
plot!(Xell, Yell, label="elliptical curve", linewidth=2)

#We can also generate files to visualize the results in Paraview using the function write_vtk.
# properties = get_surface_properties(system)
# write_vtk("mirrored-planar-wing", surfaces, properties; symmetric)
# println("Your wing has been plotted :)")