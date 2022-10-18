#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
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

        std::cout << "Master process generated sequence of size: " << sequence.size() << std::endl;
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

    // Now that every process has their respective chunks of the array, we can start the parallel quicksort algorithm
    recursiveHalving( taskId, numTasks, subarray );

    MPI_Finalize();

    return 0;
}

int partitionSequence( std::vector<int>& sequenceToPartition, int pivot, int low, int high )
{
    if (sequenceToPartition.empty())
    {
        return -1;
    }

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

        int pivot = 5000;
        if (taskId == processorIdThatWillBroadcastPivot)
        {
            // The current process has to send the pivot to all the other processes in the current group

            if (sequence.size() > 0)
            {
                pivot = sequence[sequence.size() - 1]; // the pivot will be the last element of the sequence
            }

            // Send the pivot to each of the processors in the current group
            for (int i = 1; i < nbrOfProcessorsPerGroup; ++i)
            {
                int receiverId = processorIdThatWillBroadcastPivot + i;
                MPI_Send( &pivot, 1, MPI_INT, receiverId, 0, MPI_COMM_WORLD );
            }
        }
        else
        {
            // If the current process is not the designated process to choose and send the pivot, it has to receive the pivot
            MPI_Status receiveStatus;
            MPI_Recv( &pivot, 1, MPI_INT, processorIdThatWillBroadcastPivot, 0, MPI_COMM_WORLD, &receiveStatus);
        }

        // Part 2: Partition the local sequence based on the chosen pivot
        int pivotIndex = partitionSequence( sequence, pivot, 1, sequence.size() - 1 );

        // Let's separate our current sequence into two halves
        std::vector<int> lowerHalf;
        std::vector<int> upperHalf;

        for ( int n : sequence )
        {
            if (n <= pivot) lowerHalf.push_back(n);
            else upperHalf.push_back(n);
        }

        sequence.clear(); // erase the entire sequence as we will insert into it again

        // Part 3: Determine the id of the neighbor and determine if we are in the low or high half
        int middleIndexOfGroup = nbrOfProcessorsPerGroup / 2;
        bool isLowHalf = (taskId % nbrOfProcessorsPerGroup) < middleIndexOfGroup;

        // Let's find the corresponding neighbor in the other half (upper or lower) of the current processor
        int neighborId = taskId ^ (1 << i - 1);

        // Part 4: First send the upper half and then the lower half is sent
        if (isLowHalf)
        {
            // The current processor is in the low half so we send the upper half to the neighbor
            // Send the size of the of the upper half sequence first because we don't want to send an empty sequence
            int upperHalfSequenceSize = upperHalf.size();
            MPI_Send( &upperHalfSequenceSize, 1, MPI_INT, neighborId, 0, MPI_COMM_WORLD );

            if (upperHalfSequenceSize > 0)
            {
               // We only send the upper half sequence to the neighbor if has at least one element
                MPI_Send((void*) &upperHalf[0], upperHalfSequenceSize, MPI_INT, neighborId, 0, MPI_COMM_WORLD);
            }

            // Now we will receive the lower half of the neighbor

            // Get the sequence size first
            int neighborLowerHalfSequenceSize;
            MPI_Recv((void*) &neighborLowerHalfSequenceSize, 1, MPI_INT, neighborId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (neighborLowerHalfSequenceSize > 0)
            {
                std::vector<int> otherHalf(neighborLowerHalfSequenceSize);

                // Now that we know the size of the message that we will receive, we can call the MPI_Recv function
                MPI_Recv((void*) &otherHalf[0], neighborLowerHalfSequenceSize, MPI_INT, neighborId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Insert the received half into the current sequence
                lowerHalf.insert(lowerHalf.cend(), otherHalf.begin(), otherHalf.end());
            }
            sequence.insert(sequence.cend(), lowerHalf.begin(), lowerHalf.end()); // insert the lower half that we kept
        }
        else
        {
            // The current process is in the high half

            // Obtain the size of the neighbor's upper half sequence
            int neighborUpperHalfSequenceSize;
            MPI_Recv( &neighborUpperHalfSequenceSize, 1, MPI_INT, neighborId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::vector<int> otherHalf(neighborUpperHalfSequenceSize);
            if (neighborUpperHalfSequenceSize > 0)
            {
                // Now that we know the size of the message that we will receive, we can call the MPI_Recv function
                MPI_Recv((void*)&otherHalf[0], neighborUpperHalfSequenceSize, MPI_INT, neighborId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Insert the received half into the current sequence
                upperHalf.insert(upperHalf.cend(), otherHalf.begin(), otherHalf.end());
            }
            sequence.insert(sequence.cend(), upperHalf.begin(), upperHalf.end()); // insert our upper half that we kept

            // Now we can send the lower half

            // Send the size of the lower half
            int lowerHalfSequenceSize = lowerHalf.size();
            MPI_Send(&lowerHalfSequenceSize, 1, MPI_INT, neighborId, 0, MPI_COMM_WORLD);

            if (lowerHalfSequenceSize > 0)
            {
               // We only send the sequence if it has values in it
                MPI_Send((void*)&lowerHalf[0], lowerHalfSequenceSize, MPI_INT, neighborId, 0, MPI_COMM_WORLD);
            }
        }
        ++iterationCount;
    }

    // Perform local quicksort
    std::sort( sequence.begin(), sequence.end() );

    // Child with id 0 will gather all the sorted sequences and display the final result
    if (taskId == 0)
    {
        int numOfSequencesToReceive = numTasks - 1;
        std::vector<std::vector<int>> finalResult(numTasks, std::vector<int>(0, 0));
        finalResult[0] = sequence; // We already have the sorted sequence of child 0

        int totalElementsReceived = 0;

        while (numOfSequencesToReceive > 0)
        {
            // Probe the incoming message size
            MPI_Status receiveStatus;
            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &receiveStatus);

            int arraySize;
            MPI_Get_count(&receiveStatus, MPI_INT, &arraySize);

            int sourceId = receiveStatus.MPI_SOURCE; // Extract the id of the sender
            finalResult[sourceId].resize(arraySize);
            totalElementsReceived += arraySize;

            std::cout << "Master process is going to receive a sequence of size: " << arraySize << " from process #" << sourceId << std::endl;
            MPI_Recv((void*)&finalResult[sourceId][0], arraySize, MPI_INT, sourceId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            --numOfSequencesToReceive;
        }

        std::cout << "Master process total elements received (not including master process' sequence): " << totalElementsReceived << std::endl;

        int totalSize = 0;
        std::cout << "Master process #" << taskId << " final result:\n[";
        for (int i = 0; i < finalResult.size(); i++)
        {
            totalSize += finalResult[i].size();
            for (int j = 0; j < finalResult[i].size(); ++j)
            {
                std::cout << finalResult[i][j] << ", ";
            }
        }
        std::cout << "]\n";

        std::cout << "Final result total size: " << totalSize << std::endl;

    }
    else
    {
        // Every other child will send its sequence to process 0
        MPI_Send((void*) &sequence[0], sequence.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}
