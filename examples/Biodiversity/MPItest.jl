using MPI

function main()
    MPI.Init()

    comm = MPI.COMM_WORLD
    Threads.@threads for i in 1:100
        A = rand(1000, 1000)
        A1 = inv(A)
        print("Hello world, I am rank $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm)), thread $i\n")
    end

    return MPI.Finalize()
end

main()
