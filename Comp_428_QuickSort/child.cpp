#include <iostream>
#include <vector>
#include <mpi.h>

#define CHUNK_SIZE 100

void recursiveHalving(int taskId, int numTasks);

void printVectorContents( std::vector<int> sequence, int taskId )
{
    // Print the received sequence
    std::cout << "Child #" << taskId << ":\n";
    for (int i = 0; i < sequence.size(); i++) {
        std::cout << sequence[i] << ", ";
    }
    std::cout << "\n";
}

std::vector<int> generateUnsortedSequence(int sequenceSize)
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
    int	taskId;	     // task ID - also used as seed number
    int numTasks;    // number of tasks

    /* Obtain number of tasks and task ID */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);

    MPI_Comm parentComm;
    MPI_Comm_get_parent(&parentComm);
    MPI_Status receiveStatus;

    int sequenceSize = numTasks * CHUNK_SIZE;
    std::vector<int> sequence;

    if (taskId == 0)
    {
        std::cout << "Total number of child processes in the child communicator: " << numTasks << std::endl;

        // The master task in this child world will be in charge of selecting the pivots in the quicksort algorithm
        // Before doing that, we must first receive the entire input sequence for the parallel quicksort algorithm

        // set the seed for the random number generator
        srand(time(NULL) + taskId);

        // Generate an unsorted, random sequence
        sequence = generateUnsortedSequence(sequenceSize);
    }

    // The next step is to scatter the input sequence among all the children in this world
    // It is important to place this call outside the if statement because it is a collective
    /*	
        MPI_Scatter(
        void* send_data,
        int send_count,
        MPI_Datatype send_datatype,
        void* recv_data,
        int recv_count,
        MPI_Datatype recv_datatype,
        int root,
        MPI_Comm communicator)
    */
    std::vector<int> subarray( CHUNK_SIZE );
    MPI_Scatter( sequence.data(), CHUNK_SIZE, MPI_INT, subarray.data(), CHUNK_SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    //Each Process then print their received chunks
    std::cout << "\nChild #" << taskId << " received subarray of size " << CHUNK_SIZE << std::endl;
    printVectorContents(subarray, taskId);

    std::cout << "\nChild #" << taskId << " is terminating..." << std::endl;
    MPI_Finalize();

    return 0;
}

void recursiveHalving(int taskId, int numTasks)
{
   
}
