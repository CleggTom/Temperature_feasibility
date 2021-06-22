a = zeros(10)

@show a

Threads.@threads for i = 1:10
    a[i] = Threads.threadid()
end

@show a