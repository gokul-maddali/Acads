Enter file contents here#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;

/* Read pairs of key and value from files */
int ReadInput(char *file, vector<int>& key, vector<int>& value);

/* Reduce the map */
int ReduceLocal(vector<int>& key, vector<int>& value);

/* Communication between different processors */
void Communication(vector<int>& key, vector<int>& value, vector< vector<int> >& result, int datasize);

/* Combine the outputs from each processor */
void CombOutFrmProcs(char *file, vector< vector<int> >& result, int datasize);

/* two useful function for sorting */
bool decreasing (int a, int b) {return (a > b); }
bool mySorting (pair<int, int> a, pair<int, int> b) { return (a.first > b.first); }

/* global comm variables */
int commSize;
int commRank;

int main(int argc, char* argv[]){

	/* Analysis arguments */
	if(argc < 3){
		cout << "Clarification: you should specfyed input map file and output map file" << endl;
		exit(-1);
	}

	/* Initialize MPI */
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	double StartTime, EndTime;
	StartTime = MPI_Wtime();

	vector<int> key, value;
	int dataset;
	dataset = ReadInput(argv[1], key, value);

	/* Local map reduce */
	ReduceLocal(key, value);

	/* Communicate with other processor */
	vector< vector<int> > buffer;
	Communication(key, value, buffer, dataset);

	/* Calculate the result */
	CombOutFrmProcs(argv[2], buffer, dataset);

	/* Calculate time */
	EndTime = MPI_Wtime();
	EndTime = EndTime - StartTime;
	MPI_Reduce(&EndTime, &StartTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(commRank == 0)
		cout << "The total time is " << StartTime << " (" << commSize << " cores)"<< endl;

	MPI_Finalize();

}


int ReadInput(char *file, vector<int>& key, vector<int>& value){

	/* Open file */
	fstream in_file;
	in_file.open(file, fstream::in);
	if(!in_file.is_open()){
		cout << "The file " << file << " could not be open." << endl;
		exit(-1);
	}

	int res=0, t1, t2;
	char comma, x, y;
	string inp;
	int count=0;
	while(getline(in_file,inp)){
		char myString[20]; 
		strcpy(myString, inp.c_str());
		char *p = strtok(myString, ",");
		int kv[2];
		int cnt=0;
		while (p) {
			kv[cnt]=atoi(p);
			cnt++;
			p = strtok(NULL, ",");
		}

		bool chk = ((count != 0) && (count%commSize == commRank));
		if(count>=1 && chk){
			res+=1;
			key.resize(res);
			value.resize(res);
			key[res-1]=kv[0];
			value[res-1]=kv[1];
		}
		count++;
	}

	in_file.close();
	return res;
}


int ReduceLocal(vector<int>& key, vector<int>& value) {

	/* Sort the pair (key, value) by increasing key */
	int num = key.size();
	vector<pair <int, int> > kvPair;
	key.resize(num);

	for(int i=0; i<num; i++)
		kvPair.push_back(pair<int, int>(key[i], value[i]));

	std::sort(kvPair.begin(), kvPair.end(), mySorting);

	sort(key.begin(), key.end(), decreasing);
	vector<int>::iterator it;
	it = unique(key.begin(), key.end());
	key.resize(distance(key.begin(), it));

	value.resize(0);
	value.resize(key.size());
	int curloc = 0;						// current location in key array

	/* reduce, in this case, we accumulate the values with same key */
	for(int i=0; i<num; i++){
		if(kvPair[i].first == key[curloc])
			value[curloc] += kvPair[i].second;
		else{
			curloc++;
			value[curloc] = kvPair[i].second;
		}
	}
	
	return curloc;
}

void Communication(vector<int>& key, vector<int>& value, vector< vector<int> >& result, int datasize){

	int length = datasize * 2 + 1;
	int *sendbuf = (int*)malloc(sizeof(int) * length * commSize);
	int *sendcounts = (int*)malloc(sizeof(int) * commSize);
	int *sdispls = (int*)malloc(sizeof(int) * commSize);
	int *recvbuf = (int*)malloc(sizeof(int) * length * commSize);
	int *recvcounts = (int*)malloc(sizeof(int) * commSize);
	int *rdispls = (int*)malloc(sizeof(int) * commSize);

	/* buffer initialize */
	for(int i=0; i<commSize; i++){
		sendcounts[i] = 1;
		sdispls[i] = i * length;
		recvcounts[i] = length;
		rdispls[i] = i * length;
	}

	for(int i=0; i<key.size(); i++){
		int dest = key[i] % commSize;
		int offset = sdispls[dest] + sendcounts[dest];
		sendbuf[offset] = value[i];
		sendbuf[offset+1] = key[i];
		sendcounts[dest] += 2;
	}

	/* set the lengths of each sections in the buffer */
	for(int i=0; i<commSize; i++){
		sendbuf[sdispls[i]] = sendcounts[i];
	}

	/* All to all */
	MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_INT, recvbuf, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);

	/* Assign received buffer into new result vector */
	for(int i=0; i<commSize; i++){
		int len = recvbuf[rdispls[i]];
		vector<int> temp(&recvbuf[rdispls[i]+1], &recvbuf[rdispls[i]+len]);
		result.push_back(temp);
	}
	return;
}

void CombOutFrmProcs(char *file, vector< vector<int> >& result, int datasize){
	string name(file);
	
	/* write files */
	fstream out_file;
	out_file.open(name.c_str(), fstream::out);
	if(!out_file.is_open()){
		cout << "Could not open the output file " << name << endl;
		exit(-1);
	}
	/* each vector still has sth */
	int *file_status = (int*)malloc(sizeof(int) * commSize);
	for(int i=0; i<commSize; i++)
		file_status[i] = 1;

	/* reduce from multiple vectors<int> */
	int curkey = commRank;
	int file_alive = commSize;
	int output, value;
	vector<int> temp;



	while(file_alive > 0) {
		output = 0;
		value = 0;
		for(int i=0; i<commSize; i++){
			if(!file_status[i])
				continue;
			if(result[i].back() > curkey)
				continue;
			if(result[i].back() == curkey){
				result[i].pop_back();
				value += result[i].back();
				result[i].pop_back();
				output = 1;
				if(result[i].empty()){
					file_status[i] = 0;
					file_alive--;
				}
			}
		}
		if(output){
			temp.push_back(curkey);
			temp.push_back(value);
		}
		curkey += commSize;
	}

	if(commRank!=0){
		MPI_Send(&temp[0],temp.size(),MPI_INT,0,0,MPI_COMM_WORLD);
		cout << "Exited Send: " << commRank << endl;
	}

	MPI_Status status;
	if(commRank == 0){
		int * buf;
		map<int, int> finalTemp;
		int count1, count2;
		for(int j=1; j < commSize; j++){
		       	buf = new int[50000];
		        cout << commRank << " waiting for " << j << endl;
		        MPI_Recv(&buf[0], 50000, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
		        MPI_Get_count(&status, MPI_INT, &count1);
		        cout << commRank << " has received messages of length " << count1 << " from " << j << endl;
			for(int i = 0; i < count1; i++){
				temp.push_back(buf[i]);
			}
		}
		for(int i = 0; i < temp.size(); i+=2){
			finalTemp.insert(pair<int, int>(temp[i], temp[i+1]));
		}
		for(map<int, int>::iterator it = finalTemp.begin(); it != finalTemp.end(); it++){
			out_file << "Key, Value: " << it->first << ", " << it->second << endl;
		}
	}


	out_file.close();
	return;
}

