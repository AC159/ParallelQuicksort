#include <iostream>
#include <string>
#include <mpi.h>
#include <vector>

#define MASTER 0        // task ID of master task

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "You must specify the number of tasks to spawn!" << std::endl;
        std::cout << "Usage: mpirun -np 4 pool_of_tasks 4" << std::endl;
        return 1;
    }

    int	taskId;	        // task ID - also used as seed number
    int numtasks;       // number of tasks

    // Obtain number of tasks and task ID
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
    std::cout << "MPI task " << taskId << " has started...\n";

    // The master task will dynamically spawn all the other tasks
    const int nbrOfProcessesToStart = std::stoi(argv[1]);

    // Create child processes
     /*
        int MPI_Comm_spawn(
                        const char *command,
                        char *argv[],
                        int maxprocs,
                        MPI_Info info,
                        int root,
                        MPI_Comm comm,
                        MPI_Comm *intercomm,
                        int array_of_errcodes[]
                        )
    */
    std::cout << "Master process spawning workers..." << std::endl;
    MPI_Comm childComm;
    MPI_Comm_spawn("c", argv, nbrOfProcessesToStart, MPI_INFO_NULL, MASTER, MPI_COMM_SELF, &childComm, MPI_ERRCODES_IGNORE);

    // Wait for the children processes to finish?

    // Finalize the MPI environment.
    MPI_Finalize();
}
