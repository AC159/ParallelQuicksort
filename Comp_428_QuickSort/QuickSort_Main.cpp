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
    int numTasks;       // number of tasks

    // Obtain number of tasks and task ID
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
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

    // Collect sorted sequences from the child processes
    //int numOfSequencesToReceive = numTasks;
    //std::vector<std::vector<int>> finalResult(numTasks);

    //while (numOfSequencesToReceive > 0)
    //{
    //    // Probe the incoming message size
    //    MPI_Status receiveStatus;
    //    MPI_Probe(MPI_ANY_SOURCE, 0, childComm, &receiveStatus);

    //    int arraySize;
    //    MPI_Get_count(&receiveStatus, MPI_INT, &arraySize);

    //    int sourceId = receiveStatus.MPI_SOURCE; // Extract the id of the sender
    //    finalResult[sourceId].resize(arraySize);

    //    MPI_Recv((void*) &finalResult[sourceId][0], arraySize, MPI_INT, sourceId, 0, childComm, MPI_STATUS_IGNORE);
    //    --numOfSequencesToReceive;
    //}

    //std::cout << "Master process #" << taskId << " final result:\n[";
    //for (int i = 0; i < finalResult.size(); i++) 
    //{
    //    for ( int j = 0; j < finalResult[i].size(); ++j )
    //    std::cout << finalResult[i][j] << ", ";
    //}
    //std::cout << "]\n";

    //std::cout << "Master process: All child processes have terminated!" << std::endl;

    // Wait for all child processes to terminate
    while (numTasks > 0)
    {
        int terminationCode;
        int res = MPI_Recv((void*) &terminationCode, 1, MPI_INT, MPI_ANY_SOURCE, 0, childComm, MPI_STATUS_IGNORE);
        if (res == -1)
        {
            std::cout << "Error receiving on the master process...\n";
            exit(EXIT_FAILURE);
        }
        --numTasks;
    }

    std::cout << "Master process: All child processes have terminated!" << std::endl;

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
