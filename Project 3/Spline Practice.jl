using Plots
using FLOWMath: akima, Akima, derivative, gradient

# x = [0.0, 1.0, 2.0, 3.0, 4.0] #0:pi/4:2*pi
# y = x.*sin.(x)
# xpt = 0:.1:4 #LinRange(0, 1, 20)
# ypt = Akima(x, y)

# #figure()
# scatter(x, y)
# plot!(xpt, ypt.(xpt))
# plot!(xpt, xpt.*sin.(xpt))
# #savefig("interp.svg"); nothing # hide

x0 = [10, 2.0 ,1.6, 1.0, .5, 0]
y = [x0[2], x0[3], x0[4], x0[5], x0[6]] # span, c1, c2, etc
x = [0.0, x0[1]/3, x0[1]*(2/3), x0[1]]
xpt = 0:.1:x0[1] #LinRange(0, 1, 20)
ypt = Akima(x, y)

#figure()
scatter(x, y)
plot!(xpt, ypt.(xpt), aspect_ratio=:equal)
#plot!(xpt, xpt.*sin.(xpt))