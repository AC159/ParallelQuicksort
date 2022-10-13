#include <iostream>
#include <string>
#include <mpi.h>
#include <vector>

#define MASTER 0        // task ID of master task

std::vector<int> generateUnsortedSequence( int sequenceSize )
{
    std::vector<int> sequence(sequenceSize);

    for (int i = 0; i < sequenceSize; ++i)
    {
        sequence[i] = std::rand(); // random number between 0 and RAND_MAX
    }

    return sequence;
}

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

    // set the seed for the random number generator
    srand(time(NULL) + taskId);

    // Generate an unsorted, random sequence
    std::vector<int> sequence = generateUnsortedSequence( 100 );

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

    /*	
        MPI_Scatter(
            void* send_data,
            int send_count,
            MPI_Datatype send_datatype,
            void* recv_data,
            int recv_count,
            MPI_Datatype recv_datatype,
            int root,
            MPI_Comm communicator
            )
    */

    // Chunk size that each process will receive
    // int chunkSizeCount = sequence.size() / nbrOfProcessesToStart;

    // Send the size of the sequence first and then the array
    int sequenceSize = sequence.size();
    MPI_Send((void*)&sequenceSize, 1, MPI_INT, 0, 0, childComm);

    // Send the entire array to the first child process which will be the main process that selects the pivots for the parallel quicksort algorithm
    MPI_Send((void*)&sequence[0], sequence.size(), MPI_INT, 0, 0, childComm);

    // print the sequence
    std::cout << "Sequence from master: \n";
    for (int i = 0; i < sequence.size(); i++) {
        std::cout << sequence[i] << ", ";
    }
    std::cout << "\n";

    //std::vector<int> recvArr(chunkSizeCount);
    //MPI_Scatter(&sequence, chunkSizeCount, MPI_INT, &recvArr, chunkSizeCount, MPI_INT, MASTER, childComm);

    //// Each Process then print their received chunks
    //for (int i = 0; i < chunkSizeCount; i++) {
    //    std::cout << "\nP" << taskId << ": " << recvArr[i] << "\n";
    //}

    // Finalize the MPI environment.
    MPI_Finalize();
}
