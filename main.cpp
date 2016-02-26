#define MPI_MASTERRANK 0
#define MPI_JOBLOADED -1

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::map;

template <typename T>
void transpose(vector<T> &matrix, const unsigned int &size){
	vector<T> buffer;
	for(unsigned int i=0; i<size; ++i){
		for(unsigned int j=0; j<size; ++j){
			unsigned int index = j*size+i;
			if(index<matrix.size()){
				buffer.push_back(matrix[index]);
			}
		}
	}
	matrix = buffer;
}

template <typename T>
void showMatrix(const vector<T> &matrix, const unsigned int &size){
	for(unsigned int i=0; i<size; ++i){
		for(unsigned int j=0; j<size; ++j){
			unsigned int index = i*size+j;
			if(index<matrix.size()){
				cout << matrix[index] << " ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

template <typename T>
void showData(T* data, const unsigned int &dataSize){
	for(unsigned int i=0; i<dataSize; ++i){
		cout<<data[i]<<" ";
	}
	cout<<endl;
	return;
}

void showTime(const clock_t &t){
	double ms = 1000*((double)t)/CLOCKS_PER_SEC;
	cout << (unsigned int)ms << endl;
	return;
}

bool readMatrix(std::istream &input, vector<double> &matrix, const int &size){
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			if(!input.eof()){
				double value;
				input >> value;
				matrix.push_back(value);
			}
			else{
				return false;
			}
		}
	}
	return true;
}

template<typename T>
void generateMatrix(vector<T> &matrix, const int &size){
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			matrix.push_back(T(1));
		}
	}
}

void broadCastData(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator) {
	int world_rank;
	MPI_Comm_rank(communicator, &world_rank);
	int world_size;
	MPI_Comm_size(communicator, &world_size);

	if (world_rank == root) {
		// If we are the root process, send our data to everyone
		for (int i = 0; i < world_size; i++) {
			if (i != world_rank) {
				MPI_Send(data, count, datatype, i, 0, communicator);
			}
		}
	} else {
		// If we are a receiver process, receive the data from the root
		MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
	}
}

template <typename T>
T* packData(const T* data, T* result, const unsigned int dataSize, const int message){
	result[0] = (T)message;
	
	for(unsigned int i=0; i<dataSize; ++i){
		result[i+1] = data[i];
	}

	return result;
} 

template <typename T>
T* unpackData(const T* data, T* result, const unsigned int dataSize, int &message){
	for(unsigned int i=1; i<dataSize; ++i){
		result[i-1] = data[i];
	}

	message = (int)data[0];
	return result;
}

template<typename T>
bool multiplyVector(const vector<T> &data, T &result){
	if(data.size()%2 != 0){
		return false;
	}

	unsigned int halfSize = data.size()/2;
	result = 0;

	for(unsigned int i=0; i<halfSize; ++i){
		result += data[i]*data[halfSize+i];
	}

	return true;
}

bool doMasterJob(unsigned int &matrixSize){
	vector<double> matrixA;
	vector<double> matrixB;
	vector<double> resultMatrix;

	//unsigned int matrixSize;
	//input >> matrixSize;
	resultMatrix.resize(matrixSize*matrixSize);

	//Sendign matrix size to slaves 
	broadCastData(&matrixSize, 1, MPI_UNSIGNED, MPI_MASTERRANK, MPI_COMM_WORLD);

	//Reading matrix A and B
	// if(!readMatrix(input, matrixA, matrixSize) || !readMatrix(input, matrixB, matrixSize)){
	// 	return false;
	// }	

	//Generating matrixes
	generateMatrix(matrixA, matrixSize);
	generateMatrix(matrixB, matrixSize);

	//Transposing matrix B to ease future copmutations
	transpose(matrixB, matrixSize);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	double *dataA = matrixA.data();
	double *dataB = matrixB.data();
	double dataPackage[matrixSize + 1];

	clock_t t1, t2;

	t1 = clock();	
	{
		//Sending task to slaves using the following package structure (taskid, data0, data1,...)
		for(unsigned int i=0; i<matrixSize; ++i){
			for(unsigned int j=0; j<matrixSize; ++j){
				//rank of a slave 
				int source = (i*matrixSize + j) % (world_size - 1) + 1;

				packData(dataA, dataPackage, matrixSize + 1, (int)(i*matrixSize + j));
				MPI_Send(dataPackage, matrixSize + 1, MPI_DOUBLE, source, 0 , MPI_COMM_WORLD);
				packData(dataB, dataPackage, matrixSize + 1, (int)(i*matrixSize + j));
				MPI_Send(dataPackage, matrixSize + 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);

				dataB += matrixSize;
			}
			dataB = matrixB.data();
			dataA += matrixSize;
		}

		//Notifing slaves that the job is destributed 
		packData(dataB, dataPackage, matrixSize + 1, (int)MPI_JOBLOADED);
		broadCastData(dataPackage, matrixSize + 1, MPI_DOUBLE, MPI_MASTERRANK, MPI_COMM_WORLD);

		unsigned int resultRecieved = 0;

		//Resceiving results
		while(resultRecieved < matrixSize*matrixSize){
			MPI_Recv(dataPackage, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			resultMatrix[(unsigned int)dataPackage[0]] = dataPackage[1];
			resultRecieved++;
		}
	}
	t2 = clock();

	//showMatrix(resultMatrix, matrixSize);
	showTime(t2-t1);

	return true;
}

bool doSlaveJob(const int &slaveRank){
	unsigned int dataSize;
	//Receiving size of data 
	broadCastData(&dataSize, 1, MPI_UNSIGNED, MPI_MASTERRANK, MPI_COMM_WORLD);

	int jobID = 0;

	double data[dataSize];
	double dataPackage[dataSize + 1];

	map<int,vector<double> > taskStorage;
	map<int,double> resultStorage;

	//Receiving data packages with tasks 
	while(jobID != MPI_JOBLOADED){
		MPI_Recv(dataPackage, dataSize + 1, MPI_DOUBLE, MPI_MASTERRANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		unpackData(dataPackage, data, dataSize + 1, jobID);
		//Saving tasks
		if(jobID != MPI_JOBLOADED)
			taskStorage[jobID].insert(taskStorage[jobID].end(), data, data + dataSize);
	}

	//Executing multiplication
	for(map<int,vector<double> >::iterator it = taskStorage.begin(); it != taskStorage.end(); ++it){
		double localResult = 0;
		if(multiplyVector(it->second, localResult)){
			resultStorage[it->first] = localResult;
		}
		else{
			cout<< "Data error "<<it->first<<endl;
		}
	}

	//Sending back results
	for(map<int,double>::const_iterator it = resultStorage.begin(); it != resultStorage.end(); ++it){
		packData(&(it->second), dataPackage, (unsigned int)1, (it->first));
		MPI_Send(dataPackage, 2, MPI_DOUBLE, MPI_MASTERRANK, 0 , MPI_COMM_WORLD);
	}

	return true;
}

int main(int argc, char** argv) {
	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	if(argc != 2){
		MPI_Finalize();
		return 0;
	}

	//std::ifstream input(argv[1]);
	unsigned int matrixSize = (unsigned int)atoi(argv[1]);

	if(world_rank == MPI_MASTERRANK){
		doMasterJob(matrixSize);
	}
	else{
		doSlaveJob(world_rank);
	}

	// Print off a bye world message
	// printf("Job is done by processor %s, rank %d"
	//         " out of %d processors\n",
	//         processor_name, world_rank, world_size);

	// Finalize the MPI environment.
	MPI_Finalize();
}
