@recipe function f(a::SampleArray)
    xguide --> "s"
    yguide --> "amplitude"
    size --> (1000, 300)
#     ticks --> :native
    domain(a), data(a)
end
