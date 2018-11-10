#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <curand_kernel.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>
#include <sys/time.h>
using namespace std;

__global__ void partition_step (curandState * state, unsigned long seed )
{
    int i= blockDim.x * blockIdx.x + threadIdx.x;
    curand_init (seed, i, 0, &state[i]);
} 

__global__ void randomColouring (curandState* globalState, int *degreeCount, int n, int limit){

	int i= blockDim.x * blockIdx.x + threadIdx.x;
	
	curandState localState = globalState[i];
    	float RANDOM = curand_uniform( &localState );
    	globalState[i] = localState;
    	
    	RANDOM *= (limit - 1 + 0.999999);
    	RANDOM += 1;
	
	degreeCount[i] = (int) RANDOM;
}

__global__ void conflictDetection (int *vertexArray, int *neighbourArray, int *degreeCount, int n, int m, int *detectConflict){

	int i= blockDim.x * blockIdx.x + threadIdx.x;
	
	if (i>=n){
		return;
	}
	
	int myColour = degreeCount[i];
	
	int incoming = -1, stop = -1;
	
	incoming = vertexArray[i];
	
	if (i==n-1){	
		stop = m;
	}
	
	else{
		stop = vertexArray[i+1];
	}
	
	for (int j=incoming; j<stop; j++){
		if (degreeCount[neighbourArray[j]-1] == myColour){

			detectConflict[i]=1;
			break;
		}
	}
}


__global__ void degreeCalc (int *vertexArray, int *neighbourArray, int *degreeCount, int n, int m){

	int i= blockDim.x * blockIdx.x + threadIdx.x;
	
	if (i>=n){
		return;
	}
	
	
	int incoming = -1, stop = -1;
	int diff=0;
	
	incoming = vertexArray[i];
	
	if (i==n-1){	
		stop = m;
	}
	
	else{
		stop = vertexArray[i+1];
	}

	diff = stop-incoming;
		
	atomicAdd(&degreeCount[i], diff);
	
	for (int j=incoming; j<stop; j++){
		atomicAdd(&degreeCount[neighbourArray[j]-1], 1);
	}
}

int main(int argc, char const *argv[])
{

	int n, m;
	// Enter number of vertices and edges
	cin>>n>>m;

	int h_vertexArray[n];
	int h_neighbourArray[m];
	int h_degreeCount[n];
	int h_detectConflict[n];
	
	
	// Cuda memory allocation
	size_t bytes = n*sizeof(int);
    int *d_vertexArray = NULL;
    cudaMalloc((void **)&d_vertexArray, bytes);
    
    int *d_neighbourArray = NULL;
    cudaMalloc((void **)&d_neighbourArray, m*sizeof(int));
    
    int *d_detectConflict = NULL;
    cudaMalloc((void **)&d_detectConflict, bytes);
    cudaMemset((void *)d_detectConflict, 0,bytes);
    
    int *d_degreeCount = NULL;
    cudaMalloc((void **)&d_degreeCount, bytes);
    cudaMemset((void *)d_degreeCount, 0, bytes);
    
    curandState* partition_states;
    cudaMalloc ( &partition_states, n*sizeof( curandState ) );
    	
	for (int i = 0; i < n; ++i)
	{
		/* code */
		h_vertexArray[i]=m;
	}

	int temp = 0;

	int current = 0;
	int mark = 1;
// Add the graph based on input file
	for (int i = 0; i < m; ++i)
	{
		/* code */

		int incoming;
		int end;

		cin>>incoming>>end;
		incoming++;
		end++;

		if (incoming!=mark){ 

			if (incoming == mark+1 && h_vertexArray[mark-1]!=m){ 

			}

			else{

				for (int j = mark; j<incoming; j++){ 
					h_vertexArray[j-1]=temp;
					
				}
			}
			mark = incoming;

		}

		if (incoming==current){ 
			h_neighbourArray[temp]=end;
			temp++;
		}

		else { 
			current = incoming;

			h_vertexArray[current-1]=temp;

			h_neighbourArray[temp]=end;
			temp++;
		}
	}

	
	cudaMemcpy(d_vertexArray, h_vertexArray, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_neighbourArray, h_neighbourArray, m*sizeof(int), cudaMemcpyHostToDevice);
	
	int threadsPerBlock = 512;
	int blocksPerGrid = (n + threadsPerBlock -1)/threadsPerBlock;
	struct timeval startTime;
	struct timeval endTime;
	struct timezone startZone;
	struct timezone endZone;
	long startt,endt;
	double overhead;
	cout<<threadsPerBlock<<" "<<blocksPerGrid<<endl;
	gettimeofday(&startTime,&startZone);
	// Step 0 : Calculate degree of each vertex
	degreeCalc<<<blocksPerGrid, threadsPerBlock>>>(d_vertexArray, d_neighbourArray, d_degreeCount, n, m);

	thrust::device_ptr<int> d_ptr = thrust::device_pointer_cast(d_degreeCount);
  	int max = *(thrust::max_element(d_ptr, d_ptr + n));
	
	cout<<"Max number of colours = "<<max<<endl;


	partition_step <<<blocksPerGrid, threadsPerBlock>>> ( partition_states, time(NULL) );
	
	// Step 1 - Randomly assign colours
	randomColouring<<<blocksPerGrid, threadsPerBlock>>>(partition_states, d_degreeCount, n, max);

	cudaMemcpy(h_degreeCount, d_degreeCount, n*sizeof(int), cudaMemcpyDeviceToHost);
    cout<<"randomColouring"<<endl;
	for (int i=0; i<n; i++){
		cout<<"Color of"<<i<<": "<<h_degreeCount[i]<<endl;
	}
	cout<<endl;
	conflictDetection<<<blocksPerGrid, threadsPerBlock>>>(d_vertexArray, d_neighbourArray, d_degreeCount, n, m, d_detectConflict);
	
	thrust::device_ptr<int> d_detectConflict_ptr = thrust::device_pointer_cast(d_detectConflict);
  	int count1 = thrust::reduce(d_detectConflict_ptr, d_detectConflict_ptr + n);
  	
  	cudaMemcpy(h_detectConflict, d_detectConflict, n*sizeof(int), cudaMemcpyDeviceToHost);
	
	int countnew=0;
	int old_colors[n];
	for (int i = 0; i < n; ++i)
	{
		/* code */
		old_colors[i] = -1;
	}
	for (int i=0; i<n-1; i++){
		
		if (h_detectConflict[i]==0){
			continue;
		}
		
		countnew++;
		
		bool usedColours[n];
		
		fill(usedColours, usedColours+n, false);
		

		
		int incoming = -1, stop = -1;
	
		incoming = h_vertexArray[i];
		
		stop = h_vertexArray[i+1];
		old_colors[i] = h_degreeCount[i];
		
		
		for (int j=incoming; j<stop; j++){
		

			usedColours[h_degreeCount[h_neighbourArray[j]-1]-1] = true;
		}

		
		for (int j=0; j<n; j++){
			if (usedColours[j]==false){
				h_degreeCount[i]=j+1;
				break;
			}
		}
	}
	
	
	
	if (h_detectConflict[n-1]!=0){

		bool usedColours[n];
		
		countnew++;
		
		fill(usedColours, usedColours+n, false);
		
		int incoming = -1, stop = -1;
	
		incoming = h_vertexArray[n-1];
	
		stop = m;
		
	
		for (int j=incoming; j<stop; j++){
			usedColours[h_degreeCount[h_neighbourArray[j]-1]-1] = true;
		}
		
		for (int j=0; j<n; j++){
			if (usedColours[j]==false){
				h_degreeCount[n-1]=j+1;
				break;
			}
		}
	}
	for (int i = 0; i < n; ++i)
	{
		cout<<"Colour of i from" <<i <<" "<<old_colors[i]<<":"<<h_degreeCount[i]<<endl;
	}
	

	cudaMemset((void *)d_detectConflict, 0, (n)*sizeof(int));
	
	cudaMemcpy(d_degreeCount, h_degreeCount, n*sizeof(int), cudaMemcpyHostToDevice);



	conflictDetection<<<blocksPerGrid, threadsPerBlock>>>(d_vertexArray, d_neighbourArray, d_degreeCount, n, m, d_detectConflict);
	gettimeofday(&endTime,&endZone);
	startt = startTime.tv_sec*1000000+startTime.tv_usec;
	endt = endTime.tv_sec*1000000+endTime.tv_usec;
	overhead = (endt-startt)/1000000.0;
	count1 = thrust::reduce(d_detectConflict_ptr, d_detectConflict_ptr + n);
  		
	cout<<"Count: "<<count1<<"    "<<countnew<<endl;
	cout<<"time taken is"<<overhead<<endl;

	cudaFree(d_neighbourArray);
	cudaFree(d_vertexArray);
	cudaFree(d_degreeCount);
	cudaFree(d_detectConflict);
	
	cudaDeviceReset();
	return 0; 

}
