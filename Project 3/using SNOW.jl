using SNOW
#Messing around with Cibin
#=
function func!(g, x)
    #Obj function
    f = (x[1]-1.5)^2 + 2

    #Constraint
    g[1] = 2*x[1]

    return f
end

x0 = [100.0]
ng = 1
g = zeros(ng)
xopt, fopt, info = minimize(func!, x0, g)
=#

#Example on the SNOW.jl quick start
#=
function simple!(g, x)
    # objective
    f = 4*x[1]^2 - x[1] - x[2] - 2.5

    # constraints
    g[1] = -x[2]^2 + 1.5*x[1]^2 - 2*x[1] + 1
    g[2] = x[2]^2 + 2*x[1]^2 - 2*x[1] - 4.25
    g[3] = x[1]^2 + x[1]*x[2] - 10.0

    return f
end

x0 = [1.0; 2.0]  # starting point
lx = [-5.0, -5]  # lower bounds on x
ux = [5.0, 5]  # upper bounds on x
ng = 3  # number of constraints
lg = -Inf*ones(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g
options = Options(solver=IPOPT())  # choosing IPOPT solver

xopt, fopt, info = minimize(simple!, x0, ng, lx, ux, lg, ug, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)
=#

# Practice problem 
#=
function func!(g, x)
    f = -x[1]*x[2]
    g[1] = x[1] + 4*x[2] - 240
    return f
end

x0 = [1.0; 2.0]  # starting point
lx = [0.0, 0]  # lower bounds on x
ux = [Inf, Inf]  # upper bounds on x
ng = 1  # number of constraints
lg = -Inf*ones(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g
# options = Options(solver=IPOPT())  # choosing IPOPT solver
options = Options(derivatives=ForwardAD()) # forward-mode algorithmic differentiation

xopt, fopt, info = minimize(func!, x0, ng, lx, ux, lg, ug, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)
=#