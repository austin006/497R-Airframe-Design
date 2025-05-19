function main()
    f(x) = ℯ^(-x) - x
    errorTolarance = 10^-6

    function g(x) 
        return ℯ^(-x)
    end

    println("What is your guess?")
    x = 20
    xLast = 2*x
    
    counter = 1
    while abs(x - xLast) >= errorTolarance
        counter += 1
        xLast = x
        x = g(x)
        print(counter)
        print("  ")
        print(x)
        print("  ")
        println(x - xLast)
    end
end
main()