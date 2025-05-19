#=
This is the beginner Julia Activities for the FLOW lab onboarding project.
FLOW Lab LTRAD Program
Author: Austin McGlashan
Date: October 13, 2023
=#

using Plots

""" 
    fourSeriesThicknessFormula(x, m)

    This function takes in two inputs: x, and m; and outpouts the thickness t.

"""
function fourSeriesThicknessFormula(x, m)
    5 * m * (0.2969 * sqrt(x) - 0.1260 * x - 0.3537 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

#=
Testing the fourSeriesThicknessFormula() function with for and while loops:

    m = 0.10

    for i = 0:10
        println(fourSeriesThicknessFormula((i/10), m))
        #println("x: " + (i/10) + " and t: " + fourSeriesThicknessFormula((i/10), m))
    end

    let i = 0.0
        while i < 1.0
            println(fourSeriesThicknessFormula(i, m))
            i += 0.1
        end
    end
=#

"""
    fourSeriesCamberFormula(x, p, c)

    This function takes in three inputs: x, p, and c; and outputs the camber vector z.

"""
function fourSeriesCamberFormula(x, p, c)
    if x <= p
        (c * (2*p*x - x^2)) / (p^2)
    else
        (c * (1 - 2*p + 2*p*x - x^2)) / (1-p)^2
    end
end

#=
Testing the fourSeriesCamberFormula() function with for and while loops:

    for i = 0:10
        println(fourSeriesCamberFormula((i/10), p/10, c/100))
    end

    let i = 0.0
        while i <= 1.0
            println(fourSeriesCamberFormula(i, p/10, c/100))
            i += 0.1
        end
    end
=#

"""
    airfoilCoordinate(c, p, t)

    This function takes in whole number values for inputs: c, p, t, and an array of x spanning (0,1); and outputs the upper and lower airfoil z coordinates as vectors. 

    max camber = c
    max camber position = p
    max thickness = t

"""
function fourSeriesAirfoilCoordinate(c, p, t, x)
    Zupper = fourSeriesCamberFormula.(x, p/10, c/100) + fourSeriesThicknessFormula.(x, t/100)
    Zlower = fourSeriesCamberFormula.(x, p/10, c/100) - fourSeriesThicknessFormula.(x, t/100)
    
    hcat(Zupper, Zlower)
end

#=
Testing the airfoilCoordinate() function
    println(fourSeriesAirfoilCoordinate(c, p, t, x))
=#

"""
    naca(n)

    This function takes in a NACA 4-series coordinate (n) and plots the airfoil.

"""
function naca(n)
    c = div(n, 1000)
    p = div((n % 1000), 100)
    t = n % 100
    x = collect(range(0.0, stop=1.0, step=0.01))
    #println(fourSeriesAirfoilCoordinate(c, p, t, x))
    plot!(x, fourSeriesAirfoilCoordinate(c, p, t, x), aspect_ratio=:equal)
end

naca(1224)
naca(008)