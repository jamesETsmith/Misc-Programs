def traingular_number_sum(N):
    return int((N*N+N)/2)

def main():
    from mpi4py import MPI
    import sys
    import numpy as np
    import time

    num_proc  = MPI.COMM_WORLD.Get_size()
    my_rank   = MPI.COMM_WORLD.Get_rank()
    node_name = MPI.Get_processor_name()
    comm = MPI.COMM_WORLD

    if (my_rank == 0):
        sys.stdout.write("  %d MPI Processes are now active.\n" %(num_proc))
        sys.stdout.flush()
    comm.Barrier()

    N = 1000
    dN = N/num_proc
    x = np.arange(dN*my_rank,dN*(my_rank+1)) + 1

    local_sum = np.zeros(1)
    local_sum[0] = float(np.sum(x))


    # Print local values
    comm.barrier()
    for i in range(num_proc):
        if (i == my_rank):
            if (i == 0):
                sys.stdout.write("\n\n")
                sys.stdout.write("            Local Results\n")
            sys.stdout.write(
                "  Rank %d reports (SUM): %f.\n"
                % (my_rank, local_sum))
            sys.stdout.flush()
        comm.Barrier()

    # Create global variable to hold the global sum
    global_sum = np.array([0.])

    # Reductions
    comm.Allreduce(local_sum, global_sum, op=MPI.SUM)

    # Print Global Sum
    comm.barrier()
    for i in range(num_proc):
        if (i == my_rank):
            if (i == 0):
                sys.stdout.write("\n\n")
                sys.stdout.write("            Local Results\n")
            sys.stdout.write(
                "  Rank %d reports (SUM): %f.\n"
                % (my_rank, global_sum))
            sys.stdout.flush()
        comm.Barrier()

    # Print Comparison
    comm.barrier()
    time.sleep(0.001)
    if (my_rank==0):
        sys.stdout.write("\nActual Sum: %i\n"%traingular_number_sum(N))

if __name__ == '__main__':
    main()
