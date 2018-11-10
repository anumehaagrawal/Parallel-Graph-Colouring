#include <iostream> 
#include <list> 
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include <bits/stdc++.h>
using namespace std; 
  
// A class that represents an undirected graph 
class Graph 
{
    struct timeval TimeValue_Start;
    struct timezone TimeZone_Start;
    struct timeval TimeValue_Final;
    struct timezone TimeZone_Final;
    long time_start, time_end;
    double time_overhead; 
    int V;    // No. of vertices 
    list<int> *adj;    // A dynamic array of adjacency lists 
public: 
    // Constructor and destructor 
    Graph(int V, list<int> *adj)   { this->V = V; this->adj = adj; } 
    ~Graph()       { delete [] adj; } 
  
  
    // Prints greedy coloring of the vertices 
    void greedyColoring(); 
}; 

void Graph::greedyColoring() 
{   

    // Array that conatins the color values for the vertices 0 -> V
    int color[V];
    gettimeofday(&TimeValue_Start, &TimeZone_Start);
    #pragma omp parallel num_threads(16)
    { 
    // Assign the first color to first vertex 
    color[0] = 0; 
  
    // Initialize remaining V-1 vertices as unassigned 
    #pragma omp for
    for (int u = 1; u < V; u++) 
        color[u] = -1;  // no color is assigned to u 
    #pragma omp barrier 
  
    // A temporary array to store the available colors.  
    // Value > -1 of colorMask[cr] would mean that the color cr is 
    // assigned to one of its adjacent vertices 
    int colorMask[V]; 
    #pragma omp for
    for (int cr = 0; cr < V; cr++) 
        colorMask[cr] = -1;
    #pragma omp barrier 
  
    // Assign colors to remaining V-1 vertices 
    #pragma omp parallel for schedule(dynamic)
    for (int u = 1; u < V; u++) 
    { 
        // Process all adjacent vertices and flag their colors 
        // as unavailable 
        list<int>::iterator i; 
        for (i = adj[u].begin(); i != adj[u].end(); ++i) 
            if (color[*i] != -1) 
                colorMask[color[*i]] = u; 
  
        // Find the first available color 
        int cr; 
        for (cr = 0; cr < V; cr++) 
            if (colorMask[cr] != u) 
                break; 
  
        color[u] = cr; // Assign the found color 
    }
    #pragma omp barrier 
    }
    gettimeofday(&TimeValue_Final, &TimeZone_Final);
    // print the result 
    for (int u = 0; u < V; u++) 
        cout << "Vertex " << u<< "  --->  Color "
             << color[u] << endl; 
    time_start = TimeValue_Start.tv_sec * 1000000 + TimeValue_Start.tv_usec;
    time_end = TimeValue_Final.tv_sec * 1000000 + TimeValue_Final.tv_usec;
    time_overhead = (time_end - time_start)/1000000.0;
    cout<<"\n\nTime in Seconds (T) :"<<time_overhead;
} 

// A function to generate random graph.
void GenerateRandGraphs(int NOE, int NOV)
{
    //NOE=2506;
    //NOV=100;
    int i, j, edge[NOE][2], count;
    i = 0;
    // Build a connection between two random vertex.
    while(i < NOE)
    {
        edge[i][0] = rand()%NOV+1;
        edge[i][1] = rand()%NOV+1;
 
        if(edge[i][0] == edge[i][1])
            continue;
        else
        {
            for(j = 0; j < i; j++)
            {
                if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) || (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0]))
                    i--;
            }
        }
        i++;
    }
    //adj = new list<int>[NOV];
    // Print the random graph.
    //cout<<"\nThe generated random random graph is: ";
    list <int> *adj;
    adj = new list<int>[NOV];
    for(i = 0; i < NOV; i++)
    {
        count = 0;
        //cout<<"\n\t"<<i<<"-> { ";
        for(j = 0; j < NOE; j++)
        {
            if(edge[j][0] == i+1)
            {
                int temp = edge[j][1];
                temp = temp - 1;
                //cout<<temp<<"   ";
                count++;
                adj[i].push_back(temp);
            }
            else if(edge[j][1] == i+1)
            {
                int temp = edge[j][0];
                temp = temp - 1;
                //cout<<temp<<"   ";
                count++;
                adj[i].push_back(temp);
            }
            else if(j == NOE-1 && count == 0)
            {
                //cout<<"Isolated Vertex!";
            }
        }
        //cout<<" }";
    }


    list<int>::iterator a;
    //cout<<"\n\nAdjacency List \n";
    for(int j = 0 ; j < 100 ; ++j){

        cout <<"\n";
        for (a = adj[j].begin(); a != adj[j].end(); ++a) 
            cout << *a << " ";
    }
    cout<<"\n\n";
    /*fstream file; 
    string word, t, q, filename,word1; 
    int i,j;
    // filename of the file 
    filename = "coloring.txt"; 
    list <int> *adj;
    adj = new list<int>[100]; 
    // opening file 
    file.open(filename.c_str()); 
  
    // extracting words form the file 
    while (file >> word >> word1) 
    { 
        // displaying content 
        //cout << word << word1 << endl;
        i=stoi(word);
        j=stoi(word1); 
        //cout<<i<<" "<<j<<endl;
        adj[i].push_back(j);
    } 
    NOV=100;*/
    //Creates graph using the random vertice number
    Graph g(NOV,adj);
    //Calls the function to color the graph
    g.greedyColoring();
}

  
// Driver program to test above function 
int main() 
{ 
    
    int e, v;
    e=2056;
    v=100;
    GenerateRandGraphs(e, v);
    cout << endl;
  
    return 0; 
} 
