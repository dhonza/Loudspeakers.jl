@recipe function f(a::SampleArray)
    xlabel --> "s"
    ylabel --> "amplitude"
    size --> (1000, 300)
#     ticks --> :native
    domain(a), data(a)
end
