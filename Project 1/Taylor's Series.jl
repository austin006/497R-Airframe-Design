using Plots
function main()
    #=
    x = LinRange(-2*pi, 2*pi, 100)
    plot(x, cos.(x))
    y0(x) = 1
    plot!(x, y0.(x))
    y1(x) = 1 + 0*x
    plot!(x, y1.(x))
    y2(x) = 1 + 0*x - (1/2)x^2
    plot!(x, y2.(x))
    y3(x) = 1 + 0 - (1/2)x^2 + 0*x^3
    plot!(x, y3.(x))
    y4(x) = 1 + 0 - (1/2)x^2 + 0*x^3 - ()*x^4
    plot!(x, y4.(x))
    =#
    
    x = LinRange(-2, 2, 100)
    term(n, x) = (x.^n)./factorial(n)
    plot(x, ℯ.^x)
    plot!(x, term(0, x))
    plot!(x, term(0,x) + term(1, x))
    plot!(x, term(0,x) + term(1,x) + term(2, x))
    plot!(x, term(0,x) + term(1,x) + term(2, x) + term(3,x)) 
    #=
    x = LinRange(0, 2, 100)
    term(n, x) = (x.^n)./factorial(n)
    plot(x, sqrt.(x)) 
    =#
end
main()