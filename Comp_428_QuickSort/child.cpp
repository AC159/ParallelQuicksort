#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <mpi.h>

#define CHUNK_SIZE 100

void recursiveHalving(int taskId, int numTasks, std::vector<int> sequence);

void printVectorContents( std::vector<int> sequence, int taskId )
{
    // Print the received sequence
    std::cout << "Child #" << taskId << ":\n[";
    for (int i = 0; i < sequence.size(); i++) {
        std::cout << sequence[i] << ", ";
    }
    std::cout << "]\n";
}

std::vector<int> generateUnsortedSequence(int sequenceSize)
{
    std::vector<int> sequence(sequenceSize);

    for (int i = 0; i < sequenceSize; ++i)
    {
        sequence[i] = std::rand() % 100000; // random number between 0 and 100000
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

    // Each Process then print their received chunks
    std::cout << "\nChild #" << taskId << " received subarray of size " << CHUNK_SIZE << std::endl;
    // printVectorContents(subarray, taskId);

    // Now that every process has their respective chunks of the array, we can start the parallel quicksort algorithm
    recursiveHalving( taskId, numTasks, subarray );

    std::cout << "\nChild #" << taskId << " is terminating..." << std::endl;
    MPI_Finalize();

    return 0;
}

int partitionSequence( std::vector<int>& sequenceToPartition, int pivot, int low, int high )
{
    int i = low - 1;

    for (int j = low; j <= high - 1; ++j)
    {
        // If current element is smaller than the pivot
        if (sequenceToPartition[j] < pivot) 
        {
            i++;    // increment index of smaller element
            std::swap(sequenceToPartition[i], sequenceToPartition[j]);
        }
    }
    std::swap(sequenceToPartition[i + 1], sequenceToPartition[high]);
    return i + 1;
}

void recursiveHalving(int taskId, int numTasks, std::vector<int> sequence)
{
    int hyperCubeDimension = log2( numTasks );
    std::cout << "Child #" << taskId << " will iterate " << hyperCubeDimension << " times\n";

    int iterationCount = 1;
    for (int i = hyperCubeDimension; i > 0; --i)
    {
        // Part 1: Choose and broadcast the pivot

        int currentNbrOfGroups = numTasks / pow( 2, i );
        int nbrOfProcessorsPerGroup = numTasks / currentNbrOfGroups;
        int currentGroupId = taskId / nbrOfProcessorsPerGroup; // determine the group id of this process in this iteration

        // Determine the id of the processor that is designated to choose and broadcast the pivot. This will be the first processor in each group
        int processorIdThatWillBroadcastPivot = currentGroupId * nbrOfProcessorsPerGroup; // double the group id to obtain the id the first processor in that group

        int pivot;
        if (taskId == processorIdThatWillBroadcastPivot)
        {
            // The current process has to send the pivot to all the other processes in the current group
            pivot = sequence[sequence.size() -1 ]; // the pivot will be the last element of the sequence

            // Send the pivot to each of the processors in the current group
            for (int i = 1; i < nbrOfProcessorsPerGroup; ++i)
            {
                int receiverId = processorIdThatWillBroadcastPivot + i;
                std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " is sending pivot to child #" << i << std::endl;
                MPI_Send( &pivot, 1, MPI_INT, receiverId, 0, MPI_COMM_WORLD );
            }
        }
        else
        {
            // If the current process is not the designated process to choose and send the pivot, it has to receive the pivot
            MPI_Status receiveStatus;
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " is waiting to receive pivot from child #" << processorIdThatWillBroadcastPivot << std::endl;
            MPI_Recv( &pivot, 1, MPI_INT, processorIdThatWillBroadcastPivot, 0, MPI_COMM_WORLD, &receiveStatus);
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " received pivot: " << pivot << std::endl;
        }

        // Part 2: Partition the local sequence based on the chosen pivot
        int pivotIndex = partitionSequence( sequence, pivot, 0, sequence.size() - 1 );
        // printVectorContents(sequence, taskId);

        // Part 3: Determine the id of the neighbor and determine if we are in the low or high half
        int middleIndexOfGroup = nbrOfProcessorsPerGroup / 2;
        bool isLowHalf = (taskId % nbrOfProcessorsPerGroup) < middleIndexOfGroup;

        std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " # of processors in group: " << nbrOfProcessorsPerGroup << std::endl;

        // Let's find the corresponding neighbor in the other half (upper or lower) of the current processor
        int neighborId = taskId ^ (1 << i - 1);

        // Part 4: First send the upper half and then the lower half is sent
        if (isLowHalf)
        {
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " is in the lower half\n";
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " sending upper half data to neighbor child #" << neighborId << std::endl;

            std::vector<int> upperHalf( sequence.cbegin() + pivotIndex + 1, sequence.cend() );

            // The current processor is in the low half so we send the upper half to the neighbor
            MPI_Send( (void*) &upperHalf[0], upperHalf.size(), MPI_INT, neighborId, 0, MPI_COMM_WORLD);

            // Now we will receive the lower half of the neighbor
            MPI_Status receiveStatus;
            MPI_Probe(neighborId, 0, MPI_COMM_WORLD, &receiveStatus);

            int arraySize;
            MPI_Get_count(&receiveStatus, MPI_INT, &arraySize);

            std::vector<int> otherHalf(arraySize);

            // Now that we know the size of the message that we will receive, we can call the MPI_Recv function
            std::cout << "Iteration #" << iterationCount <<  " Child #" << taskId << " waiting for lower half data from neighbor child #" << neighborId << std::endl;
            MPI_Recv((void*)&otherHalf[0], arraySize, MPI_INT, neighborId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Erase the upper half in the sequence
            sequence.erase(sequence.cbegin() + pivotIndex + 1, sequence.cend());

            // Insert the received half into the current sequence
            sequence.insert(sequence.cend(), otherHalf.begin(), otherHalf.end());
        }
        else
        {
            // The current process is in the high half

            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " is in the upper half\n";
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " waiting for upper half data from neighbor child #" << neighborId << std::endl;
            MPI_Status receiveStatus;
            MPI_Probe( neighborId, 0, MPI_COMM_WORLD, &receiveStatus );

            int arraySize;
            MPI_Get_count(&receiveStatus, MPI_INT, &arraySize);
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " waiting to receive " << arraySize << " data from neighbor child #" << neighborId << std::endl;

            std::vector<int> otherHalf(arraySize);

            // Now that we know the size of the message that we will receive, we can call the MPI_Recv function
            MPI_Recv( (void*) &otherHalf[0], arraySize, MPI_INT, neighborId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Now we can send the lower half
            std::cout << "Iteration #" << iterationCount << " Child #" << taskId << " sending lower half data to neighbor child #" << neighborId << std::endl;
            MPI_Send((void*) &sequence[0], pivotIndex, MPI_INT, neighborId, 0, MPI_COMM_WORLD);
            
            // Erase the low half in the sequence
            sequence.erase( sequence.cbegin(), sequence.cbegin() + pivotIndex + 1 );

            // Insert the received half into the current sequence
            sequence.insert(sequence.cend(), otherHalf.begin(), otherHalf.end());
        }
        ++iterationCount;

        // printVectorContents(sequence, taskId);
    }

    // Perform local quicksort
    std::cout << "Child #" << taskId << " is performing local quicksort\n";
    std::sort( sequence.begin(), sequence.end() );

    // Child with id 0 will gather all the sortied sequences and display the final result
    if (taskId == 0)
    {
        int numOfSequencesToReceive = numTasks - 2;
        std::vector<std::vector<int>> finalResult(numTasks);
        finalResult[0] = sequence; // We already have the sorted sequence of child 0

        while (numOfSequencesToReceive >= 0)
        {
            // Probe the incoming message size
            MPI_Status receiveStatus;
            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &receiveStatus);

            int arraySize;
            MPI_Get_count(&receiveStatus, MPI_INT, &arraySize);

            int sourceId = receiveStatus.MPI_SOURCE; // Extract the id of the sender

            MPI_Recv((void*) &finalResult[sourceId], arraySize, MPI_INT, sourceId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            --numOfSequencesToReceive;
        }

        printVectorContents(sequence, taskId);
    }
    else
    {
        // Every other child will send its sequence to process 0
        MPI_Send((void*)&sequence[0], sequence.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    // printVectorContents( sequence, taskId );

    std::cout << "Child #" << taskId << " has finished!\n";
}
