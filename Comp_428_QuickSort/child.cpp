#include <iostream>
#include <vector>
#include <mpi.h>

int main(int argc, char* argv[])
{
    int	taskId;	     // task ID - also used as seed number
    int numTasks;    // number of tasks

    /* Obtain number of tasks and task ID */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);

    MPI_Comm parentComm;
    std::cout << "Child getting parent communicator..." << std::endl;
    MPI_Comm_get_parent(&parentComm);
    MPI_Status receiveStatus;

    if (taskId == 0)
    {
        std::cout << "Total number of child processes in the child communicator: " << numTasks << std::endl;

        // The master task in this child world will be in charge of selecting the pivots in the quicksort algorithm
        // Before doing that, we must first receive the entire input sequence for the parallel quicksort algorithm

        // First, we will receive the size of the input sequence that we must sort and then the actual sequence
        int sequenceSize;
        MPI_Recv((void*)&sequenceSize, 1, MPI_INT, 0, 0, parentComm, &receiveStatus);
        std::cout << "Child process received sequence size: " << sequenceSize << "\n";

        // Now we will receive the actual sequence
        std::vector<int> sequence(sequenceSize);
        MPI_Recv((void*)&sequence[0], sequenceSize, MPI_INT, 0, 0, parentComm, &receiveStatus);
        std::cout << "Child process received sequence \n";

        // Print the received sequence
        std::cout << "Child #" << taskId << ":\n";
        for (int i = 0; i < sequence.size(); i++) {
            std::cout << sequence[i] << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "\nChild task #" << taskId << " is terminating..." << std::endl;
    MPI_Finalize();

    return 0;
}