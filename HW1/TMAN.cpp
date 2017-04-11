#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <sstream>
#include "gnuplot_i.hpp"


using namespace std;

class node{
	public:
	int nid; // nodeID
	float theta ; //theta
	//int r; // radius of circle
	int *neighbors;
};

class cnode{
	public:
	int nid; // nodeID
	float theta ; //theta
	int radius; // radius of circle node belongs
	int *neighbors;
};

class bNode {
	public: 
	int nid;
	int level;
	int *neighbors;
};

class DynamicRing{
	
	public:
	int DynamicRingFun(int N, int k, int n, int r[]);
	void Initialize(int N, int k, int r[],node nodes[]);
	void Evolution(int N, int k, int r[],int n, node nodes[]);
	void UpdateNeighbors(int snid, int rnid, int k, int r, node nodes[]); //snid = sender nid, rnid = receiver nid
	float Distance(int n1, int n2, int r, node nodes[]);
	void PrintNeighbors(int N, int k, int r, string filename1, string filename2, node nodes[]);
	double SumofDistance(int N, int k, int r, node nodes[]);
	void PrintDistance(int cycleno, double sumofdistance, string filename1, string filename2);
	void wait_for_key ();
	int SumofDistGraph(int N, int k);
	int RingGraph(int N, int k, int cycle);
};

class BinaryTree{
	public:
	int BinaryTreeFun(int N, int k);
	int binaryInitialize(int N, int k, bNode bnodes[]);
	int binaryEvolution(int N, int k, bNode bnodes[]);
	int getLevel(int index);
	int Distance(int snid, int rnid); //snid - sending nodeID, rnid - Receiving nodeID
	int UpdateNeighbors(int snid, int rnid, int k, bNode bnodes[]);
	int SumofDistance(int N, int k, bNode bnodes[]);
	void PrintNeighbors(int N, int k, string filename, bNode bnodes[]);
	void PrintDistance(int cycleno, int sumofdistance, string filename1, string filename2);
	void wait_for_key ();
	int SumofDistGraph(int N, int k);
};

class CrescentMoon{
	public:
	float r1;// = 10; //assume radius of first circle is 10
	float r2; //= 15; // assume radius of second circle as 15
	
	float centre_x1;// = 0, //centre_y1 = 0;//centre of circle 1
	float centre_y2;// = 0, //centre_x2 = sqrt((r2*r2) - (r1*r1));
	float centre_y1;
	float centre_x2;
	CrescentMoon();
	
	int CrescentMoonFun(int N, int k); 
	int crescentInitialize(int N,int k, cnode cnodes[]);
	int crescentEvolution(int N, int k, cnode cnodes[]);
	int UpdateNeighbors(int snid, int rnid,int N, int k, cnode cnodes[]);
	float Distance(int snid, int rnid,  cnode cnodes[]);
	float EuDistance(float x1, float y1, float x2, float y2);
	void PrintNeighbors(int N, int k, string filename1, string filename2, cnode cnodes[]);
	double SumofDistance(int N, int k, cnode cnodes[]);
	void PrintDistance(int cycleno, int sumofdistance, string filename1, string filename2);
	void wait_for_key ();
	int SumofDistGraph(int N, int k);
	float getX(int nid, cnode cnodes[]);
	float getY(int nid, cnode cnodes[]);
	int CrescentGraph(int N, int k, int cycle);
};

int CrescentMoon :: CrescentMoonFun(int N, int k){
		cout<<"Inside CrescentMoon";	//testprint
		cnode cnodes[N];
		crescentInitialize(N,k,cnodes);
		crescentEvolution(N,k,cnodes);
		SumofDistGraph( N,  k);
}

CrescentMoon::CrescentMoon(){
	float r1 = 10; //assume radius of first circle is 10
	float r2 = 15; // assume radius of second circle as 15
	float rad_ratio = r2/r1;
	float centre_x1 = 0, centre_y1 = 0;//centre of circle 1
	float centre_y2 = 0, centre_x2 = sqrt((r2*r2) - (r1*r1));//centre of circle 2
}
int CrescentMoon :: crescentInitialize(int N,int k, cnode cnodes[]){
	
	float r1,r2;
	float theta1, theta2, rad_ratio;
	
	r1 = 10; //assume radius of first circle is 10
	r2 = 15; // assume radius of second circle as 15
	rad_ratio = r2/r1;
	float centre_x1 = 0, centre_y1 = 0;//centre of circle 1
	float centre_y2 = 0, centre_x2 = sqrt((r2*r2) - (r1*r1));//centre of circle 2
	
	cout <<"Centre of Circle two is ("<<centre_x2<<","<<centre_y2<<")"<<endl;
	
	float angle1_arc1 = M_PI/2, angle2_arc1 = 3*M_PI/2;
	cout<<"Arc angle 1 of circle 1 is "<<angle1_arc1<<" and "<<angle2_arc1<<endl;
	
	float angle1_arc2 = M_PI - asin(r1/r2), angle2_arc2 = (M_PI) + asin(r1/r2); 
	cout<<"Arc angle 1 of circle 2 is "<<angle1_arc2<<" and "<<angle2_arc2<<endl;
	
	float arc_len1 = r1*(angle2_arc1 - angle1_arc1), arc_len2 = r2*(angle2_arc2 - angle1_arc2);
	cout<<"Arc length of arc 1 is "<<arc_len1<<" arc 2 is"<<arc_len2<<endl;
	
	float len_ratio = arc_len1/arc_len2;
	cout<<"ratio of arc lengths is "<<len_ratio<<"\n";
	
	int numNodesArc1 = floor((arc_len1/(arc_len2 + arc_len1))*N);
	int numNodesArc2 = N - numNodesArc1;
	cout<<"No of nodes in arc1: "<<numNodesArc1<<" No of nodes in arc2: "<<numNodesArc2;
	
	theta1 = (angle2_arc1 - angle1_arc1) / numNodesArc1;
	theta2 = (angle2_arc2 - angle1_arc2) / numNodesArc2;
	cout<<"\ntheta1 is "<<theta1<<" theta2 is "<<theta2<<"\n";
	for (int i=0;i<N;i++){
		cnodes[i].nid = i;
		cnodes[i].neighbors = new int[k];
	}
	//cnodes[0].theta = M_PI/2;
	cout<<"numNodesArc1 "<<numNodesArc1<<endl;
	for (int i=0;i<numNodesArc1;i++){
		cnodes[i].theta = angle1_arc1 + i*theta1;
		cnodes[i].radius = r1;
	}
	int j = numNodesArc1;
	for (int i=1;i<=numNodesArc2;i++){
		cnodes[j].theta = angle1_arc2 + i*theta2;
		cnodes[j].radius = r2;
		j++;
	}
	
	
	int ranarray[N];
	
	int ran;
	for (int m=0;m<N;m++){
		ranarray[m] = m;
	}
	
	for (int i=0;i<N;i++){
		
		random_shuffle(&ranarray[0],&ranarray[N]);
		ran = rand()%(N-k+1);
		
		
		for (int j=0; j<k; j++){
			
			if(ranarray[ran] != i){
				cnodes[i].neighbors[j] = ranarray[ran];
				ran++;		
				//ranarray[ran] = i;
			}
			else{
				ran++;		
				cnodes[i].neighbors[j] = ranarray[ran];	
				ran++;
			}
			
			
		}
		
	}
}

int CrescentMoon :: crescentEvolution(int N, int k, cnode cnodes[]){
	int index;
	int randselnid;
	int sumofdistance;
	
	
	
	for (int j=1;j<=50;j++){
		
		for (int i=0;i<N;i++){
						
			index = rand()%k; 
			
			randselnid = cnodes[i].neighbors[index];//randomly selected nodeid
			
			
			UpdateNeighbors(i,randselnid,N,k, cnodes);//fun to send/receive and update neighbors list for two nodes 
				
		}
		
		cout <<j<<"th Cycle Completed..."<<endl;
		
		if (j == 1 || (j%5 == 0 && j<=15) ){
			
			ostringstream oss;
			oss << "C_N"<<N<<"_k"<<k<<"_"<<j<<".txt";
			string filename1 = oss.str();
			
			ostringstream oss1;
			oss1 << "C_N"<<N<<"_k"<<k<<"_"<<j<<"plot.txt";
			string filename2 = oss1.str();
			
			PrintNeighbors( N,  k, string(filename1), string(filename2), cnodes);
			CrescentGraph( N,  k,  j);
		}
		sumofdistance = SumofDistance(N, k,  cnodes);
	
		ostringstream oss;
		oss << "C_N"<<N<<"_k"<<k<<".txt";
		string filename1 = oss.str();
		
		ostringstream oss1;
		oss1 << "C_N"<<N<<"_k"<<k<<"_SumofDistance.txt";
		string filename2 = oss1.str();
		PrintDistance(j, sumofdistance,filename1,filename2);
	}
}


int CrescentMoon :: UpdateNeighbors(int snid, int rnid,int N, int k, cnode cnodes[]){
	
	int newlists[(2*k)+1];
	int newlistr[(2*k)+1];
	float distlists[(2*k)+1];
	float distlistr[(2*k)+1];
	float swap;
	int a=0; int i=0; int flags = 0, flagr =0;
	
	for (int i=0;i<k;i++){
		newlists[i] = cnodes[rnid].neighbors[i];
		newlistr[i] = cnodes[rnid].neighbors[i];
	}
	for (int i=k,j=0;i<(2*k);i++,j++){
		newlists[i] = cnodes[snid].neighbors[j];
		newlistr[i] = cnodes[snid].neighbors[j];
	}
	newlistr[2*k] = cnodes[snid].nid;
	newlists[2*k] = cnodes[rnid].nid;
	

	for(int i=0;i<=(2*k);i++){
		distlists[i] = Distance(snid,newlists[i], cnodes);
		distlistr[i] = Distance(rnid,newlistr[i], cnodes);
	} 
	

	for(int c=0;c<=((2*k));c++){
		for (int d=0;d<((2*k)-c);d++){
			if(distlists[d]>distlists[d+1]){
				swap = distlists[d];
				distlists[d] = distlists[d+1];
				distlists[d+1] = swap;
				
				swap = newlists[d];
				newlists[d] = newlists[d+1];
				newlists[d+1] = swap;
			}	

			if(distlistr[d]>distlistr[d+1]){
				swap = distlistr[d];
				distlistr[d] = distlistr[d+1];
				distlistr[d+1] = swap;
				
				swap = newlistr[d];
				newlistr[d] = newlistr[d+1];
				newlistr[d+1] = swap;
							
			}	
		}
	}
	
	for (int i=0;i<2*k;i++){
		//cout << "i in for is***************** "<<i<<endl;
		if(distlists[i] == distlists[i+1] && newlists[i] > newlists[i+1]){
			swap = newlists[i];
			newlists[i] = newlists[i+1];
			newlists[i+1] = swap;
		}
		if(distlistr[i] == distlistr[i+1] && newlistr[i] > newlistr[i+1]){
			swap = newlistr[i];
			newlistr[i] = newlistr[i+1];
			newlistr[i+1] = swap;
		}
	}
	
    
	
	a = 0; i = 0;
	while((a<k) && (i<2*k)){
		
		if(a==0 && (newlists[i] != snid) ){
			cnodes[snid].neighbors[a] = newlists[i];
			i++; a++;
			
		}
		for (int p=0;p<a;p++){
			if (newlists[i] == cnodes[snid].neighbors[p] && a!=1){
				flags = 1;
				//cout<<"Flags is Set";	
			}
		}
		if (flags == 1 && i<(2*k)){
			i++;
		}
		if((newlists[i] != snid) && (newlists[i] != cnodes[snid].neighbors[a-1])){
			cnodes[snid].neighbors[a] = newlists[i];	
			a++;
			i++;		
		}	
		else{i++;}
		
	}
	
	a = 0;
	i = 0;
	while((a<k) && (i<2*k)){
		
		if(a==0 && (newlistr[i] != rnid) ){
			cnodes[rnid].neighbors[a] = newlistr[i];
			i++; a++;
			
		}
		for (int p=0;p<a;p++){
			if (newlistr[i] == cnodes[rnid].neighbors[p] && a!=1){
				flagr = 1;
				//cout<<"Flagr is Set";
				
			}
		}
		if (flagr == 1 && i<(2*k)){
			i++;
		}
		if((newlistr[i] != rnid) && (newlistr[i] != cnodes[rnid].neighbors[a-1])){
			cnodes[rnid].neighbors[a] = newlistr[i];	
			a++;
			i++;		
		}	
		else{i++;}
		
	}

}

float CrescentMoon :: Distance(int snid, int rnid, cnode cnodes[]){
	
	
	float theta,x1,y1,x2,y2,top_x,top_y,bot_x,bot_y,dist_top,dist_bot,D1,D2;
	x1 = centre_x1 + cnodes[snid].radius * cos(cnodes[snid].theta);
	x2 = centre_x2 + cnodes[rnid].radius * cos(cnodes[rnid].theta);
	y1 = centre_y1 + cnodes[snid].radius * cos(cnodes[snid].theta);
	y2 = centre_y2 + cnodes[rnid].radius * cos(cnodes[rnid].theta);
	
	if(x1 == top_x && y1 == top_y){
		return EuDistance(top_x,top_y,x2,y2);
	}
	else if(x1 == bot_x && y1 == bot_y){
		return EuDistance(bot_x,bot_y,x2,y2);
	}
	else if(x2 == top_x && y2 == top_y){
		return EuDistance(top_x,top_y,x1,y1);
	}
	else if(x2 == bot_x && y2 == bot_y){
		return EuDistance(bot_x,bot_y,x1,y1);
	}
	else if(cnodes[snid].radius == cnodes[rnid].radius){
		return EuDistance(x1,y1,x2,y2);	
	}
	
	else{
		D1 = EuDistance(x1,y1,top_x,top_y) +  EuDistance(x2,y2,top_x,top_y);
		D2 = EuDistance(x1,y1,bot_x,bot_y) + EuDistance(x2,y2,bot_x,bot_y);
		
		return min(D1,D2);
	}
		
}
float CrescentMoon :: EuDistance(float x1, float y1, float x2, float y2){
	return sqrt(pow((x1-x1),2) + pow((y1-y2),2));
}

void CrescentMoon :: PrintNeighbors(int N, int k, string filename1, string filename2, cnode cnodes[])
{
	//string file = filename;
	int Nnid;
	ofstream myfile (filename1.c_str());
	
	ofstream myfile2(filename2.c_str());
	if (myfile.is_open() && myfile2.is_open()){
  		for (int i=0;i<N;i++){
			myfile << "Neighbor List of "<<i<<"th node is: \n";
			
			for (int j=0;j<k;j++){
				myfile<<cnodes[i].neighbors[j]<<"\t";
				
				myfile2<<getX(i,cnodes)<<" "<<getY(i,cnodes)<<"\n";
				Nnid = cnodes[i].neighbors[j];
				
				myfile2<<getX(Nnid,cnodes)<<" "<<getY(Nnid,cnodes)<<"\n\n";
				
			}
			myfile<<"\n";
		}
    myfile.close();
	myfile2.close();
  	}
}		
float CrescentMoon :: getX(int nid, cnode cnodes[]){
	float ret;
	
	float centre_x1 = 0, centre_x2 = 11.1803;
	float r1 = 10;
	float r2 = 15;
	if (cnodes[nid].radius == r1){
		
		ret =  centre_x1 + r1*cos(cnodes[nid].theta);
	}
	else if (cnodes[nid].radius ==r2){
		
		ret =  centre_x2 + r2*cos(cnodes[nid].theta);
	}
	return ret;
}

float CrescentMoon :: getY(int nid, cnode cnodes[]){
	float centre_y1 = 0, centre_y2 = 0, ret;
	float r1 = 10;
	float r2 = 15;
	if (cnodes[nid].radius == r1){
		ret =  centre_y1 + r1*sin(cnodes[nid].theta);
	}
	else if (cnodes[nid].radius ==r2){
		ret =  centre_y2 + r2*sin(cnodes[nid].theta);
	}
	return ret;
}

double CrescentMoon :: SumofDistance(int N, int k, cnode cnodes[]){
	
	double NeighborDistSum;
	NeighborDistSum = 0;
	for (int i=0;i<N;i++){
		
		for (int j=0;j<k;j++){
			
			NeighborDistSum = NeighborDistSum + Distance(i,cnodes[i].neighbors[j],cnodes);
			
		}
	}
	
	return NeighborDistSum;
}



int BinaryTree :: BinaryTreeFun(int N, int k){
	
	bNode bnodes[N];
	binaryInitialize(N,k,bnodes);
	binaryEvolution(N,k,bnodes);
	SumofDistGraph( N,  k);
}

int BinaryTree :: binaryInitialize(int N, int k, bNode bnodes[]){
	
	for (int i=0;i<N;i++){
		bnodes[i].nid = i;
		bnodes[i].level = floor(log2(i)) + 1;
		bnodes[i].neighbors = new int[k];
	}
	int ranarray[N];
	
	int ran;
	for (int m=0;m<N;m++){
		ranarray[m] = m;
	}
	
	for (int i=0;i<N;i++){
		
		random_shuffle(&ranarray[0],&ranarray[N]);
		ran = rand()%(N-k+1);
		
		
		for (int j=0; j<k; j++){
			
			if(ranarray[ran] != i){
				bnodes[i].neighbors[j] = ranarray[ran];
				ran++;		
				//ranarray[ran] = i;
			}
			else{
				ran++;		
				bnodes[i].neighbors[j] = ranarray[ran];	
				ran++;
			}
			
			
		}
		cout<<"\n\n";
	}
}

int BinaryTree :: getLevel(int index){
	return 	floor(log2(index+1)) + 1;
}

int BinaryTree :: binaryEvolution(int N, int k, bNode bnodes[]){
	int index;
	int randselnid;
	int sumofdistance;
	for (int cycle = 1;cycle<=50;	cycle++){
		cout<<"This is "<<cycle<<"th cycle "<<"\n";
		
		for (int i=0;i<N;i++){
						
			index = rand()%k; 
			
			randselnid = bnodes[i].neighbors[index];//randomly selected nodeid
			
			
			UpdateNeighbors(i,randselnid,k,bnodes);//fun to send/receive and update neighbors list for two nodes 
				
		}
		
		
	
		if (cycle == 1 || (cycle%5 == 0 && cycle<=15)){
			
				ostringstream oss;
				oss << "B_N"<<N<<"_k"<<k<<"_"<<cycle<<".txt";
				string filename = oss.str();
				PrintNeighbors( N,  k, string(filename), bnodes);
			}
			sumofdistance = SumofDistance(N, k,  bnodes);
		
			ostringstream oss;
			oss << "B_N"<<N<<"_k"<<k<<".txt";
			string filename1 = oss.str();
		
			ostringstream oss1;
			oss1 << "B_N"<<N<<"_k"<<k<<"_SumofDistance.txt";
			string filename2 = oss1.str();
		
			PrintDistance(cycle, sumofdistance,filename1,filename2);
		}
	}
	
	
		


int BinaryTree	::	UpdateNeighbors(int snid, int rnid, int k, bNode bnodes[]){
	
	int newlists[(2*k)+1];
	int newlistr[(2*k)+1];
	int distlists[(2*k)+1];
	int distlistr[(2*k)+1];
	int swap;
	int a=0; int i=0; int flags = 0, flagr =0;
	
	for (int i=0;i<k;i++){
		newlists[i] = bnodes[rnid].neighbors[i];
		newlistr[i] = bnodes[rnid].neighbors[i];
	}
	
	for (int i=k,j=0;i<(2*k);i++,j++){
		newlists[i] = bnodes[snid].neighbors[j];
		newlistr[i] = bnodes[snid].neighbors[j];
	}
	
	newlistr[2*k] = bnodes[snid].nid;
	newlists[2*k] = bnodes[rnid].nid;
	
	
	for (int i=0; i<=(2*k);i++){
		
		distlists[i] = Distance(snid+1,newlists[i]+1);
		distlistr[i] = Distance(rnid+1,newlistr[i]+1);
	}
	
	//cout<<"After making distanceList";
	for(int c=0;c<=((2*k));c++){
		for (int d=0;d<((2*k)-c);d++){
			if(distlists[d]>distlists[d+1]){
				swap = distlists[d];
				distlists[d] = distlists[d+1];
				distlists[d+1] = swap;
				
				swap = newlists[d];
				newlists[d] = newlists[d+1];
				newlists[d+1] = swap;
			}	

			if(distlistr[d]>distlistr[d+1]){
				swap = distlistr[d];
				distlistr[d] = distlistr[d+1];
				distlistr[d+1] = swap;
				
				swap = newlistr[d];
				newlistr[d] = newlistr[d+1];
				newlistr[d+1] = swap;
							
			}	
		}
	}
	
	for (int i=0;i<2*k;i++){
		//cout << "i in for is***************** "<<i<<endl;
		if(distlists[i] == distlists[i+1] && newlists[i] > newlists[i+1]){
			swap = newlists[i];
			newlists[i] = newlists[i+1];
			newlists[i+1] = swap;
		}
		if(distlistr[i] == distlistr[i+1] && newlistr[i] > newlistr[i+1]){
			swap = newlistr[i];
			newlistr[i] = newlistr[i+1];
			newlistr[i+1] = swap;
		}
	}
	
    
	
	a = 0; i = 0;
	while((a<k) && (i<2*k)){
		
		if(a==0 && (newlists[i] != snid) ){
			bnodes[snid].neighbors[a] = newlists[i];
			i++; a++;
			
		}
		for (int p=0;p<a;p++){
			if (newlists[i] == bnodes[snid].neighbors[p] && a!=1){
				flags = 1;
					
			}
		}
		if (flags == 1 && i<(2*k)){
			i++;
		}
		if((newlists[i] != snid) && (newlists[i] != bnodes[snid].neighbors[a-1])){
			bnodes[snid].neighbors[a] = newlists[i];	
			a++;
			i++;		
		}	
		else{i++;}
		
	}
	
	a = 0;
	i = 0;
	while((a<k) && (i<2*k)){
		
		if(a==0 && (newlistr[i] != rnid) ){
			bnodes[rnid].neighbors[a] = newlistr[i];
			i++; a++;
			
		}
		for (int p=0;p<a;p++){
			if (newlistr[i] == bnodes[rnid].neighbors[p] && a!=1){
				flagr = 1;
				
				
			}
		}
		if (flagr == 1 && i<(2*k)){
			i++;
		}
		if((newlistr[i] != rnid) && (newlistr[i] != bnodes[rnid].neighbors[a-1])){
			bnodes[rnid].neighbors[a] = newlistr[i];	
			a++;
			i++;		
		}	
		else{i++;}
		
	}
	
}

int BinaryTree	:: Distance(int a, int b){
	//cout<<"Inside Distance function";
	int bits = 32;
	int alevel=bits;
	int blevel=bits;
	int commonprefix=0;
	int mask = 1 << bits-1;
	// find the level of node a
	while( (mask & a) == 0 )
	{
		a <<= 1;
		alevel--;
	}
	// find the level of node b
	while( (mask & b) == 0 )
	{
		b <<= 1;
		blevel--;
	}
	int length = min(alevel,blevel);
	while( (mask & ~(a ^ b)) != 0 && length>0)
	{
		b <<= 1;
		a <<= 1;
		commonprefix++;
		length--;
	}
	return alevel - commonprefix + blevel - commonprefix;
}

void BinaryTree :: PrintNeighbors(int N, int k, string filename, bNode bnodes[])
{
	//string file = filename;
	ofstream myfile (filename.c_str());
	if (myfile.is_open()){
  		for (int i=0;i<N;i++){
			myfile << "Neighbor List of "<<i+1<<"th node is: \n";
			for (int j=0;j<k;j++){
				myfile<<bnodes[i].neighbors[j]+1<<"\t";
			}
			myfile<<"\n";
		}
    myfile.close();
  	}
}	

int BinaryTree :: SumofDistance(int N, int k, bNode bnodes[]){
	
	int NeighborDistSum;
	NeighborDistSum = 0;
	for (int i=0;i<N;i++){
		
		for (int j=0;j<k;j++){
			
			NeighborDistSum = NeighborDistSum + Distance(i+1,bnodes[i].neighbors[j]+1);
			
		}
	}
	
	return NeighborDistSum;
}

void BinaryTree ::	PrintDistance(int cycleno, int sumofdistance, string filename1, string filename2){
	ofstream myfile (filename1.c_str(),ofstream::out|ofstream::app);
		myfile<<"Sum of Distance in "<<cycleno<<"th cycle is: "<<sumofdistance<<"\n";
		
		if (cycleno == 50){
			myfile.close();
		}
	ofstream myfile2 (filename2.c_str(),ofstream::out|ofstream::app);
		myfile2<<cycleno<<" "<<sumofdistance<<"\n";
		
		if (cycleno == 50){
			myfile2.close();
		}
}

void CrescentMoon ::	PrintDistance(int cycleno, int sumofdistance, string filename1, string filename2){
	ofstream myfile (filename1.c_str(),ofstream::out|ofstream::app);
		myfile<<"Sum of Distance in "<<cycleno<<"th cycle is: "<<sumofdistance<<"\n";
		
		if (cycleno == 50){
			myfile.close();
		}
	ofstream myfile2 (filename2.c_str(),ofstream::out|ofstream::app);
		myfile2<<cycleno<<" "<<sumofdistance<<"\n";
		
		if (cycleno == 50){
			myfile2.close();
		}
}

int DynamicRing :: DynamicRingFun(int N, int k, int n, int r[]){
		node nodes[N];
		
		
		
		Initialize(N,k,r,nodes); //Initialisation of network
		
		Evolution(N,k,r,n,nodes);//Network Evolution
		SumofDistGraph(N,k);
		
		
	}
void DynamicRing :: Initialize(int N, int k, int r[],node nodes[]){
	
	int ranarray[N];
	
	int ran;
	for (int m=0;m<N;m++){
		ranarray[m] = m;
	}
	
	

	for (int i=0;i<N;i++){
		nodes[i].nid = i; 
		nodes[i].theta = ((2*M_PI)/N)*i;
		
		nodes[i].neighbors = new int[k];
	}	
	
	for (int i=0;i<N;i++){
		
		random_shuffle(&ranarray[0],&ranarray[N]);
		ran = rand()%(N-k+1);
		for (int j=0; j<k; j++){
			if(ranarray[ran] != i){
				nodes[i].neighbors[j] = ranarray[ran];
				ran++;		
			}
			else{
				ran++;		
				nodes[i].neighbors[j] = ranarray[ran];	
				ran++;
			}
			
		}
		
	}
}

void DynamicRing :: Evolution(int N, int k, int r[],int n, node nodes[]){
	cout<<"Inside the Evolution Function"; //testprint
	
	int r_target;
	int radius;
	int temp = 0;
	r_target = r[temp++];
	radius = r_target;
	int randselnid,index;
	double sumofdistance;
	
	
	for (int j=1;j<=50;j++){ // 50 cycles for evolution of the network
		
		if ((j%3==0) && (radius < r_target)){
			
			radius++;		
		}
		if ((j%5 == 0) && (temp<n)){
			
			r_target = r[temp++];
			 
		}
		
		for (int i=0;i<N;i++){
						
			index = rand()%k; 
			
			randselnid = nodes[i].neighbors[index];//randomly selected nodeid
			UpdateNeighbors(i,randselnid,k,radius,nodes);//fun to send/receive and update neighbors list for two nodes 
				
		}
		cout<<j<<"th Cycle of Evolution Completed ...\n";
		if (j == 1 || (j%5 == 0 && j<=15)){
			
			ostringstream oss;
			oss << "D_N"<<N<<"_k"<<k<<"_"<<j<<".txt";
			string filename = oss.str();
			ostringstream oss1;
			oss1 << "D_N"<<N<<"_k"<<k<<"_"<<j<<"plot.txt";
			string filename2 = oss1.str();
			PrintNeighbors( N,  k, radius, string(filename),string(filename2) ,nodes);
			RingGraph(N,k,j);
		}
		sumofdistance = SumofDistance(N, k, radius,  nodes);
		
		ostringstream oss;
		oss << "D_N"<<N<<"_k"<<k<<".txt";
		string filename1 = oss.str();
		
		ostringstream oss1;
		oss1 << "D_N"<<N<<"_k"<<k<<"_SumofDistance.txt";
		string filename2 = oss1.str();
		PrintDistance(j, sumofdistance,filename1,filename2);
		
		
	}
	
	
}

void DynamicRing :: UpdateNeighbors(int snid, int rnid, int k, int r, node nodes[]){
	int newlists[(2*k)+1];
	int newlistr[(2*k)+1];
	float distlists[(2*k)+1];
	float distlistr[(2*k)+1];
	float swap;
	int a = 0,i=0, flags=0,flagr=0;
	for (int i=0;i<k;i++){
		newlists[i] = nodes[rnid].neighbors[i];
		newlistr[i] = nodes[rnid].neighbors[i];
	}
	for (int i=k,j=0;i<(2*k);i++,j++){
		newlists[i] = nodes[snid].neighbors[j];
		newlistr[i] = nodes[snid].neighbors[j];
	}
	newlistr[2*k] = nodes[snid].nid;
	newlists[2*k] = nodes[rnid].nid;
	

	for(int i=0;i<=(2*k);i++){
		distlists[i] = Distance(snid,newlists[i],r, nodes);
		distlistr[i] = Distance(rnid,newlistr[i],r, nodes);
	} 
	

	for(int c=0;c<=((2*k));c++){
		for (int d=0;d<((2*k)-c);d++){
			if(distlists[d]>distlists[d+1]){
				swap = distlists[d];
				distlists[d] = distlists[d+1];
				distlists[d+1] = swap;
				
				swap = newlists[d];
				newlists[d] = newlists[d+1];
				newlists[d+1] = swap;
			}	

			if(distlistr[d]>distlistr[d+1]){
				swap = distlistr[d];
				distlistr[d] = distlistr[d+1];
				distlistr[d+1] = swap;
				
				swap = newlistr[d];
				newlistr[d] = newlistr[d+1];
				newlistr[d+1] = swap;
							
			}	
		}
	}
	
	a = 0; i = 0;
	while((a<k) && (i<2*k)){
		
		if(a==0 && (newlists[i] != snid) ){
			nodes[snid].neighbors[a] = newlists[i];
			i++; a++;
			
		}
		for (int p=0;p<a;p++){
			if (newlists[i] == nodes[snid].neighbors[p] && a!=1){
				flags = 1;
			}
		}
		if (flags == 1 && i<(2*k)){
			i++;
		}
		if((newlists[i] != snid) && (newlists[i] != nodes[snid].neighbors[a-1])){
			nodes[snid].neighbors[a] = newlists[i];	
			a++;
			i++;		
		}	
		else{i++;}
		
	}
	
	a = 0;
	i = 0;
	while((a<k) && (i<2*k)){
		
		if(a==0 && (newlistr[i] != rnid) ){
			nodes[rnid].neighbors[a] = newlistr[i];
		}
		for (int p=0;p<a;p++){
			if (newlistr[i] == nodes[rnid].neighbors[p] && a!=1){
				flagr = 1;
			}
		}
		if (flagr == 1 && i<(2*k)){
			i++;
		}
		if((newlistr[i] != rnid) && (newlistr[i] != nodes[rnid].neighbors[a-1])){
			nodes[rnid].neighbors[a] = newlistr[i];	
			a++;
			i++;		
		}	
		else{i++;}
	}
		
	}



float DynamicRing :: Distance(int n1, int n2, int r, node nodes[])
{
	
	float x1 = r*(cos(nodes[n1].theta));
	float x2 = r*(cos(nodes[n2].theta));
	float y1 = r*(sin(nodes[n1].theta));
	float y2 = r*(sin(nodes[n2].theta));
	
	return sqrt(((x1-x2)*(x1-x2))+((y1-y2)*(y1-y2)));
	
}

void DynamicRing :: PrintNeighbors(int N, int k, int r, string filename1, string filename2, node nodes[])
{
	//string file = filename;
	int Nnid;
	ofstream myfile (filename1.c_str());
	
	ofstream myfile2(filename2.c_str());
	if (myfile.is_open() && myfile2.is_open()){
  		for (int i=0;i<N;i++){
			myfile << "Neighbor List of "<<i<<"th node is: \n";
			
			for (int j=0;j<k;j++){
				myfile<<nodes[i].neighbors[j]<<"\t";
				myfile2<<r*cos(nodes[i].theta)<<" "<<r*sin(nodes[i].theta)<<"\n";
				Nnid = nodes[i].neighbors[j];
				//cout<<"Nnid is "<<Nnid<<"\t";
				myfile2<<r*cos(nodes[Nnid].theta)<<" "<<r*sin(nodes[Nnid].theta)<<"\n\n";
			}
			myfile<<"\n";
		}
    myfile.close();
	myfile2.close();
  	}
}	

double DynamicRing :: SumofDistance(int N, int k, int r, node nodes[]){
	
	double NeighborDistSum;
	NeighborDistSum = 0;
	for (int i=0;i<N;i++){
		
		for (int j=0;j<k;j++){
			
			NeighborDistSum = NeighborDistSum + Distance(i,nodes[i].neighbors[j],r, nodes);
			
		}
	}
	
	return NeighborDistSum;
}

int DynamicRing :: RingGraph(int N, int k, int cycle){
	try{
		cout<<"Sum of Distance: \n";
		ostringstream fname, gname;
		
		fname<<"D_N"<<N<<"_"<<"k"<<k<<"_"<<cycle<<"plot.txt";
		string png;
		Gnuplot gp("lines");
		//gp.set_ylabel("Sum of Distances");
		//gp.set_xlabel("Cycles");
		gp.set_title("Dynamic Ring Topology ");
		gp.plotfile_xy(fname.str().c_str(),1,2);
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)		
		cout<<"Press Enter to Continue..."<<endl;
		
		cin.clear();
		cin.ignore(std::cin.rdbuf() -> in_avail());
		cin.get();
#endif
		
		gname<<"D_N"<<N<<"_"<<"k"<<k<<"_"<<cycle<<".png";
		gp.cmd("set terminal png\n");
		png = string("set output") + "'" + gname.str() + "'\n";
		gp.cmd(png);
		gp.replot();
	}
	catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
}

int CrescentMoon :: CrescentGraph(int N, int k, int cycle){
	try{
		//cout<<"Sum of Distance: \n";
		ostringstream fname, gname;
		
		fname<<"C_N"<<N<<"_"<<"k"<<k<<"_"<<cycle<<"plot.txt";
		string png;
		Gnuplot gp("lines");
		//gp.set_ylabel("Sum of Distances");
		//gp.set_xlabel("Cycles");
		gp.set_title("Crescent Moon Topology ");
		gp.plotfile_xy(fname.str().c_str(),1,2);
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)		
		cout<<"Press Enter to Continue..."<<endl;
		
		cin.clear();
		cin.ignore(std::cin.rdbuf() -> in_avail());
		cin.get();
#endif
		
		gname<<"C_N"<<N<<"_"<<"k"<<k<<"_"<<cycle<<".png";
		gp.cmd("set terminal png\n");
		png = string("set output") + "'" + gname.str() + "'\n";
		gp.cmd(png);
		gp.replot();
	}
	catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
}

int DynamicRing :: SumofDistGraph(int N, int k){
	try{
		cout<<"Sum of Distance: \n";
		ostringstream fname, gname;
		
		fname<<"D_N"<<N<<"_"<<"k"<<k<<"_SumofDistance.txt";
		string png;
		Gnuplot gp("lines");
		gp.set_ylabel("Sum of Distances");
		gp.set_xlabel("Cycles");
		gp.set_title("Dynamic Ring Topology - Sum of Distances per Cycle");
		gp.plotfile_xy(fname.str().c_str(),1,2);
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)		
		cout<<"Press Enter to Continue..."<<endl;
		
		cin.clear();
		cin.ignore(std::cin.rdbuf() -> in_avail());
		cin.get();
#endif
		
		gname<<"C_N"<<N<<"_"<<"k"<<k<<"_"<<".png";
		gp.cmd("set terminal png\n");
		png = string("set output") + "'" + gname.str() + "'\n";
		gp.cmd(png);
		gp.replot();
	}
	catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
}

int BinaryTree :: SumofDistGraph(int N, int k){
	try{
		cout<<"Sum of Distance: \n";
		ostringstream fname, gname;
		
		fname<<"B_N"<<N<<"_"<<"k"<<k<<"_SumofDistance.txt";
		string png;
		Gnuplot gp("lines");
		gp.set_ylabel("Sum of Distances");
		gp.set_xlabel("Cycles");
		gp.set_title("Binary Tree Topology - Sum of Distances per Cycle");
		gp.plotfile_xy(fname.str().c_str(),1,2);
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)		
		cout<<"Press Enter to Continue..."<<endl;
		
		cin.clear();
		cin.ignore(std::cin.rdbuf() -> in_avail());
		cin.get();
#endif
		
		gname<<"B_N"<<N<<"_"<<"k"<<k<<".png";
		gp.cmd("set terminal png\n");
		png = string("set output") + "'" + gname.str() + "'\n";
		gp.cmd(png);
		gp.replot();
	}
	catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
}

int CrescentMoon :: SumofDistGraph(int N, int k){
	try{
		cout<<"Sum of Distance: \n";
		ostringstream fname, gname;
		
		fname<<"C_N"<<N<<"_"<<"k"<<k<<"_SumofDistance.txt";
		string png;
		Gnuplot gp("lines");
		gp.set_ylabel("Sum of Distances");
		gp.set_xlabel("Cycles");
		gp.set_title("Crescent Moon Topology - Sum of Distances per Cycle");
		gp.plotfile_xy(fname.str().c_str(),1,2);
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)		
		cout<<"Press Enter to Continue..."<<endl;
		
		cin.clear();
		cin.ignore(std::cin.rdbuf() -> in_avail());
		cin.get();
#endif
		
		gname<<"C_N"<<N<<"_"<<"k"<<k<<".png";
		gp.cmd("set terminal png\n");
		png = string("set output") + "'" + gname.str() + "'\n";
		gp.cmd(png);
		gp.replot();
	}
	catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
}
void DynamicRing ::	PrintDistance(int cycleno, double sumofdistance, string filename1, string filename2){
	ofstream myfile (filename1.c_str(),ofstream::out|ofstream::app);
		myfile<<"Sum of Distance in "<<cycleno<<"th cycle is: "<<sumofdistance<<"\n";
		
		if (cycleno == 50){
			myfile.close();
		}
	  
	ofstream myfile2 (filename2.c_str(),ofstream::out|ofstream::app);
		myfile2<<cycleno<<" "<<sumofdistance<<"\n";
		
		if (cycleno == 50){
			myfile2.close();
		}
}

void DynamicRing :: wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void BinaryTree :: wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void CrescentMoon :: wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

int main(int argc,char *argv[])
{  
	int N = atoi(argv[1]); //Number of nodes in the system
	int k = atoi(argv[2]); //number of neighbors to be maintained by each node
	char topology = *argv[3]; //topology
	
	DynamicRing top1; //create object of Dynamic Ring class
	BinaryTree  top2;
	CrescentMoon top3;
	int status;
	
	if(topology == 'D'){
	int n = atoi(argv[4]); //Number of radius
	int r[n]; // array containg all radius
	
	int j = 0;
	char ch;
	
	char *tok = strtok(argv[5],",");

	while(tok !=NULL){
		r[j] = atoi(tok);
		tok = strtok(NULL,",");
		j++;
	}
		
	cout<<"Number of nodes in the system: "<<N<<endl;
	cout<<"Number of neighbors maintained by each node: "<<k;
	cout<<"Network Topology: "<<topology<<endl;
	cout<< "Number of radii: "<<n<<endl;
	
	
		status = top1.DynamicRingFun(N,k,n,r);
	}
	
	
	else if(topology == 'B') {
		status = top2.BinaryTreeFun(N,k);
	}
	else if(topology == 'C'){
		status = top3.CrescentMoonFun(N,k);
	}

	return 0;	
}
