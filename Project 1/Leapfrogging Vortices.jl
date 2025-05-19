# This is the Leapfrogging Vortices project for the FLOW Lab onboarding project.
# FLOW Lab LTRAD Program
# Author: Austin McGlashan
# Date: October 13-17, 2023

#=
Planning the program structure...
    repeat for the number of timeSteps
        repeat for each point (four times)
            repeat for each other point (three times)
                calculate distance between points
                find Net velocity
                Put the Net velocity from each point into an array of x and y components
            end

            add all three velocities to find net velocity
            find final position = net velocity * Δt + initial position
            add final position to end of position array containing list of positions
        end
    end
    plot all the x and y values for each timestep
=#

using LinearAlgebra
using Plots

#Declaration of initial variables
gammaP14 = [0, 0, 1]
gammaP23 = [0, 0, -1]
d = 1.0
Δt = 0.01
timeSteps = 4000

positionP1 = [0.0, -0.5, 0]
positionP2 = [0.0, 0.5, 0]
positionP3 = [1.0, 0.5, 0]
positionP4 = [1.0, -0.5, 0]

"""
    calculateNetVelocity()

    This function takes in one primary point in addition to three points and their according gamma values.
    The output is the net net induced velocity of the other three points on the first point.
"""
function calculateNetVelocity(a, b, c, d, gammaB, gammaC, gammaD)
    velocity1 = calculateInducedVelocity(a, b, gammaB)
    velocity2 = calculateInducedVelocity(a, c, gammaC)
    velocity3 = calculateInducedVelocity(a, d, gammaD)

    velocity1 .+ velocity2 .+ velocity3
end

"""
    calculateInducedVelocity()

    This function takes in point a, point b, and the gamma value of point b
    It outputs the induced velocity of point b on point a
"""
function calculateInducedVelocity(a, b, gammaB)
    xDistance = b[1, end] - a[1,end]
    yDistance = b[2, end] - a[2,end]

    r = sqrt((xDistance)^2 + (yDistance)^2)
    vectorR = [xDistance, yDistance, 0]

    cross(gammaB, vectorR) ./ (2*pi*(r^2))
end

#For each timestep find the new position of each point and add it to the end of the list of positions
for i = 1:timeSteps 
    change1 = Δt .* calculateNetVelocity(positionP1, positionP2, positionP3, positionP4, gammaP23, gammaP23, gammaP14)
    change2 = Δt .* calculateNetVelocity(positionP2, positionP1, positionP3, positionP4, gammaP14, gammaP23, gammaP14)
    change3 = Δt .* calculateNetVelocity(positionP3, positionP2, positionP1, positionP4, gammaP23, gammaP14, gammaP14)
    change4 = Δt .* calculateNetVelocity(positionP4, positionP2, positionP3, positionP1, gammaP23, gammaP23, gammaP14)
    
    global positionP1 = hcat(positionP1, positionP1[:, end] .+ change1)
    global positionP2 = hcat(positionP2, positionP2[:, end] .+ change2)
    global positionP3 = hcat(positionP3, positionP3[:, end] .+ change3)
    global positionP4 = hcat(positionP4, positionP4[:, end] .+ change4)
end

#Plot all of the x and y values of the points at each timestep
plot(xlims = (0, 12), ylims = (-2.1, 2.1), zlim = (0, 5), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], yticks = [0,-2, 2])

plot!(positionP1[1, :], positionP1[2, :])
plot!(positionP2[1, :], positionP2[2, :])
plot!(positionP3[1, :], positionP3[2, :])
plot!(positionP4[1, :], positionP4[2, :]) 

#=
#Plot the points as an animation
@gif for t = 1:timeSteps 
    change1 = Δt .* calculateNetVelocity(positionP1, positionP2, positionP3, positionP4, gammaP23, gammaP23, gammaP14)
    change2 = Δt .* calculateNetVelocity(positionP2, positionP1, positionP3, positionP4, gammaP14, gammaP23, gammaP14)
    change3 = Δt .* calculateNetVelocity(positionP3, positionP2, positionP1, positionP4, gammaP23, gammaP14, gammaP14)
    change4 = Δt .* calculateNetVelocity(positionP4, positionP2, positionP3, positionP1, gammaP23, gammaP23, gammaP14)
    
    global positionP1 = hcat(positionP1, positionP1[:, end] .+ change1)
    global positionP2 = hcat(positionP2, positionP2[:, end] .+ change2)
    global positionP3 = hcat(positionP3, positionP3[:, end] .+ change3)
    global positionP4 = hcat(positionP4, positionP4[:, end] .+ change4)

    plot(positionP1[1, t], positionP1[2, t])
    plot!(positionP2[1, t], positionP2[2, t], aspect_ratio=:equal)
    plot!(positionP3[1, t], positionP3[2, t], aspect_ratio=:equal)
    plot!(positionP4[1, t], positionP4[2, t], aspect_ratio=:equal) 
end every 10
=#