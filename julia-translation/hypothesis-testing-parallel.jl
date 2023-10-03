include("quivers.jl");

Quivers = [GeneralizedKroneckerQuiver(i) for i in range(3,20)]
d = [5,11]
theta = [11,-5]

Threads.@threads for Q in Quivers

    HN = all_harder_narasimhan_types(Q, d, theta)
    AmplyStable = is_amply_stable(Q, d, theta)

    # technical condition: it implies AS, but is it equivalent?
    candidates_pseudo_AS = filter(e -> e != ZeroVector(number_of_vertices(Q)) && e != d && slope(e,theta) >= slope(d-e,theta),all_subdimension_vectors(d))
    pseudo_Ample_Stability = all(e -> euler_form(Q,e,d-e) <= -2, candidates_pseudo_AS)
    # technical condition: my condition on ext groups
    MyHypothesis = all(hntype -> euler_form(Q,first(hntype),last(hntype)) < 0, filter(e -> length(e) > 1, HN))

    log = false
    if AmplyStable && !pseudo_Ample_Stability
        log = true
    end
    if !MyHypothesis
        log = true
    end
    if pseudo_Ample_Stability && !MyHypothesis
        log = true
    end
    if log
        println("quiver: ", Q, "Is not ok!")
        println("my hypothesis: ", MyHypothesis)
        println("amply stable: ", AmplyStable)
        println("pseudo ample stability: ", pseudo_Ample_Stability)
        println("\n")
    else
        println("quiver: ", Q, " is ok")
    end
end 