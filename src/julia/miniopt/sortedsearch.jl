using BenchmarkTools
a = rand(Int, 10000)
b = rand(Int)
_a = sort(a)


function search(a, b)

    if b in a 
        return 0
    end

end

function sortedsearch(a, b)
    if b in a
        return 0
    end

end

@benchmark search(a, b)
