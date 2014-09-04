#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string>
#include <hash_map>
#include <utility>
#include <fstream>
#include <list>

//included GNU Name Space for using Hash Map collection
using namespace __gnu_cxx;

using namespace std;

//
typedef int element;

//Structure of Vertex containg a value if used. vVal is used for storing edge weights incident on them
struct vertex{
	element vVal;
};

//temporary list used to print out to the std out, the vertex values and the corresponding edge weights 
list<int> tempList;

//Fibonacci Heap Node containg all the 8 items specified in the class used as a node fibonacci heap
class minFibNode{
	private:
		int degree, nVal; 					//used for storing degree of the node and edge weight respectively
		bool isRChildCut; 					//used for child cut check
		int refArrayInd;					//used for reference index from the array
		minFibNode* nLeft;					//used as pointer to left node of current node
		minFibNode* nRight;					//used as pointer to right node of current node
		minFibNode* nChild;					//used as pointer to child node
		minFibNode* nParent;				//used as pointer parent node
	public:
		void setDegree(int);				//used to set the degree	
		void setValue(int);					//used to set the value
		void setChildCut(bool);				//used to set the child cut
		void setIndex(int);					//used to set index
		void setLeft(minFibNode*);			//used to set left node pointer
		void setRight(minFibNode*);			//used to set right node pointer
		void setChild(minFibNode*);			//used to set child node pointer
		void setParent(minFibNode*);		//used to set parent node pointer
		minFibNode* createDefNode();		//used to create a default node
		int getDegree();					//used to get degree of the node in the heap
		int getValue();						//used to get value of the node in the heap
		int getIndex();						//used to get index of the node in the reference array of heap
		bool getChildCut();					//used to get degree of the node in the heap
		minFibNode* getLeft();				//used to get left node pointer of the node in the heap
		minFibNode* getRight();				//used to get right node pointer of the node in the heap
		minFibNode* getChild();				//used to get child node pointer of the node in the heap
		minFibNode* getParent();			//used to get parent node pointer of the node in the heap
};

//random function to generate random numbers
int random(int range)
{
	int i = rand() % range;
	return i;
}

//Graph class containg vertex array and map array corresponding to each vertex
class CGraph{
public:
	vertex *Vertices;   																//array of vetices of the graph
	
	hash_map<int,int> *edgeMap;															//array of Hash Maps corresponding to each vertex to maintain adjecency list
	
	int numEdge;																		//number of edges of the graph
				
	int numVer;																			//number of vertices of the graph
	
	void initializeGraph(int n, int m, int ver1Ind[], int ver2Ind[], int uEdgeCost[]){  //init method to create user mode graph 
	
		numEdge = m;
		
		numVer = n;
		
		this->edgeMap = new hash_map<int,int>[numVer];
		
		this->Vertices = new vertex[numVer];
		
		for(int verIter = 0; verIter < numVer; verIter++){								//for loop setting all vVal of vertex to 0
			this->Vertices[verIter].vVal = 0;
		}
		
		for(int edgeIter = 0; edgeIter < numEdge; edgeIter++){							//edge iterator to include all edges in adjecency list of each of the vertex
		
			int l1 = ver1Ind[edgeIter];
			
			int l2 = ver2Ind[edgeIter];
			
			int edgeCost = uEdgeCost[edgeIter];
			
			this->edgeMap[l1].insert(pair<int, int>(l2,edgeCost));
			this->edgeMap[l2].insert(pair<int, int>(l1,edgeCost));			
		}
				
	}
	
	
	void initializeGraph(int n, int maxEdgeCount, int density){							//init method overloaded to create graph in random mode
	
		numEdge = maxEdgeCount * density/100;
		
		numVer = n;
		
		this->edgeMap = new hash_map<int,int>[numVer];
		this->Vertices = new vertex[numVer];
		
		for(int verIter = 0; verIter < numVer; verIter++){								//for loop setting all vVal of vertex to 0
			this->Vertices[verIter].vVal = 0;
		}			
		
		srand (time(NULL));																//seek random to generate different values each time a rand() function is called
		
		for(int edgeIter = 0; edgeIter < numEdge; edgeIter++){							//for loop used to generate random edges
		
			int l1,l2;
			
			int edgeCost = random(1000)+1;
			
			do{																			
				l1 = random(n);
				l2 = 0;
				do{
					l2 = random(n);
				}while (l2 == l1);														//while loop used to check if random values generated are not equal
			}while(edgeMap[l1].find(l2)!=edgeMap[l1].end());							//check made if randomly generated edge is already present
			
			this->edgeMap[l1].insert(pair<int, int>(l2,edgeCost));
			this->edgeMap[l2].insert(pair<int, int>(l1,edgeCost));	
		}
		}
};


//Depth First Search Method performed to check if a graph is connected
int DFS(CGraph G, int verIndex){
	
	int nextIndex;
	
	int verCount = 0;
	
	hash_map<int,int>::iterator it;
	
	//if no edges are generated
	if(G.edgeMap[verIndex].begin() == G.edgeMap[verIndex].end()){
		return 0;
	}
	
	//marking the edge visited by setting vVal of the corresponding edge to 1
	G.Vertices[verIndex].vVal = 1;
	
	//iterator to perform DFS
	it = G.edgeMap[verIndex].begin();
	
	//while loop to do a recursive DFS to get the vertex count
	while (it != G.edgeMap[verIndex].end()) {
		int l2 = it->first;
		if(G.Vertices[l2].vVal != 1){
				verCount = verCount + DFS(G,l2);
		}
		++it;
	};
	
	//returning verCount in each DFS call to get the total vetrtex count
	return (verCount+1);
}


//to check if the random generated graph is connected or not by checking if all the vertices are visited
int isConnected(CGraph G){
	
	int l1 = random(G.numVer);	
	
	int verCount = DFS(G, l1);	
	
	if (verCount < G.numVer)
		return 0;
	return 1;
}


//create method to create random graph
CGraph createGraph(int numVert, int maxEdgeCount, int density){
	
	CGraph rG;
	
	rG.initializeGraph(numVert, maxEdgeCount, density);
	
	while (!isConnected(rG)){
		delete [] rG.edgeMap;
		delete [] rG.Vertices;
		rG.initializeGraph(numVert, maxEdgeCount, density);	
	};
	
	return rG;	
}


//Create user Graph function to create user graph
CGraph createUserGraph(int numVert,int numEdg, int ver1Ind[], int ver2Ind[], int uEdgeCost[]){
	
	CGraph uG;
	
	uG.initializeGraph(numVert,numEdg,ver1Ind,ver2Ind,uEdgeCost);
	
	return uG;
}

//Used in simple scheme to pick the minimum element by traversing the array of vertices
int simpleMinPick(CGraph G){
	
	int verInd;
	
	int runMin = 1001;																		//Marking values 1001 as the maximum values of the edge is 1000
	
	for(int i = 0; i < G.numVer; i++){
		if(runMin > G.Vertices[i].vVal && G.Vertices[i].vVal > 0){
			runMin = G.Vertices[i].vVal;
			verInd = i;
		}
	}
	
	return verInd;
	
}

//Generates MST using simple array scheme
int MST(CGraph G,int verIndex,int frmVerInd[],int mstCount){
	
	int edgeValSum = 0;
	
	G.Vertices[verIndex].vVal = 0;
	
	hash_map<int,int>::iterator it;
	
	//check made to see if vertex count exceeded number of vertices
	if (mstCount < G.numVer-1){
		int nextIndex;		
		
		it = G.edgeMap[verIndex].begin();	
		
		if(G.edgeMap[verIndex].begin() == G.edgeMap[verIndex].end()){
			return 0;
		}
		
		vertex* vTemp;
		
		//update all the vertex values if their current min is less than old min
		while(it != G.edgeMap[verIndex].end()){
			
			if(it->second < G.Vertices[it->first].vVal){
				G.Vertices[it->first].vVal = it->second;
				frmVerInd[it->first] = verIndex;
			}	
			
			++it;
		}
		
		//next Index to add to the MST
		nextIndex = simpleMinPick(G);
		
		//List collection used to send output to the standard output
		tempList.push_back(nextIndex);
		tempList.push_back(frmVerInd[nextIndex]);
		tempList.push_back(G.Vertices[nextIndex].vVal);
		
		//running count used to check vertex count
		mstCount++;
		
		//running edge val sum and call to recursive MST generation
		edgeValSum = G.Vertices[nextIndex].vVal + MST(G,nextIndex,frmVerInd,mstCount);
	}
	
	return edgeValSum;
}

//prototype for fibonacci MST generation as it is used in the Prim's algorithm
int fibMST(CGraph,int*);

//Prim's Algorithm
void primMst(CGraph G,char algo){
	
	int l1 = random(G.numVer);	     								//random vertex to start
	int frmVerInd[G.numVer];										//array of vertex indices to track the minimum edge value source vertex for a corresponding vertex
	
	list<int>::iterator it;
	
	//setting all vVals of vertices to maximum as they track the edge values in MST function
	for(int i = 0; i < G.numVer; i++){
		G.Vertices[i].vVal = 1001;  
		frmVerInd[i] = G.numVer; 
	}
	
	frmVerInd[l1] = l1;
	
	//setting vVal of random vetex to 0 so that it is added 1st to MST
	G.Vertices[l1].vVal = 0;	
	
	int edgeCount;
	int tempCount = 1;
	
	if(algo == 's' || algo == 'r'){
		
		//clock_t Start1, Time1;
		//Start1 = clock();
		edgeCount = MST(G,l1,frmVerInd,0);
		//Time1 = clock() - Start1;
		
		//cout << Time1 << " s Scheme" << endl;
		//Output if user simple scheme is requested
		if(algo == 's'){
		    cout << edgeCount << endl;		
		    
			for(it = tempList.begin(); it!=tempList.end(); it++){
				cout << *it;
				
				if(tempCount%3 != 0)
					cout << " ";
				else
					cout << endl;
				tempCount++;
			}
		
		}
    }
    
	if(algo == 'f' || algo == 'r'){
		
		//clock_t Start2, Time2;
		//Start2 = clock();
		edgeCount = fibMST(G, frmVerInd);
		//Time2 = clock() - Start2;
		//cout << Time2 << " f Scheme" << endl;
		//Output to std out from temp List if user fibonacci is called
		
		if(algo == 'f'){			
	    cout << edgeCount << endl;	
			
		for(it = tempList.begin(); it!=tempList.end(); it++){
			cout << *it;
			if(tempCount%3 != 0)
				cout << " ";
			else
				cout << endl;
			tempCount++;
			}
		}
	}

}

//Default node function to initialize minim fibonacci Heap Node
minFibNode* minFibNode::createDefNode()
{
		minFibNode* nDef;
		nDef = new minFibNode;
		nDef->setDegree(0);
		nDef->setValue(-1);
		nDef->setChildCut(false);
		nDef->setChild(NULL);
		nDef->setIndex(0);
		nDef->setLeft(NULL);
		nDef->setParent(NULL);
		nDef->setRight(NULL);
		return nDef;
}

//following funcions are explained above for their respective uses in the node class
void minFibNode::setChildCut(bool b){
		isRChildCut = b;
}

void minFibNode::setChild(minFibNode* child){
	nChild = child;
}

void minFibNode::setDegree(int num){
	degree = num;
}

void minFibNode::setIndex(int r){
	refArrayInd = r;
}

void minFibNode::setLeft(minFibNode* left){
	nLeft = left;
}

void minFibNode::setRight(minFibNode* right){
	nRight = right;
}

void minFibNode::setParent(minFibNode* parent){
	nParent = parent;
}

void minFibNode::setValue(int num){
	nVal = num;
}

int minFibNode::getDegree(){
	return degree;
}

int minFibNode::getValue(){
	return nVal;
}

int minFibNode::getIndex(){
	return refArrayInd;
}

bool minFibNode::getChildCut(){
	return isRChildCut;
}

minFibNode* minFibNode::getLeft(){
	return nLeft;
}

minFibNode* minFibNode::getRight(){
	return nRight;
}

minFibNode* minFibNode::getChild(){
	return nChild;
}

minFibNode* minFibNode::getParent(){
	return nParent;
}


//Minimum Fibonacci Heap Class to do the fibonacci min pick in Prmi's algorithm
class minFibHeap
{
   private : 
   
		 minFibNode *min, *root;														//Min and Root pointer for heap and each min tree respectively
		 
		 minFibNode *dummy;																//Dummy node to initialize the heap
		 
		 minFibNode **refArray, **degT;													//Reference array storing the vertex index
		 
		 int refArrayInd, maxRefArrayInd;												//Reference array index and maximum value of reference array index
		 

   public :
   	
		 // Init function constructor to create min Fib Heap
		 minFibHeap(int num)
		 {
			min = dummy->createDefNode();
		        root = dummy->createDefNode();
			this->createAndInitRefArray(num);
			refArrayInd = 0;
			maxRefArrayInd = num;
		 }		
		 
		 void createAndInitRefArray(int);												//used to initialize reference array
		 
		 void pairwiseCombine(minFibNode*,minFibNode*);									//pairwise Combine function to combine 2 min trees
		 
		 minFibNode* insert(minFibNode*,int);											//insert function to insert into the top level circular list of nodes
		 
	   	 int removeMin(minFibNode*);													//remove Min function to remove the minimum pointer Node
	   	 
		 int arbitraryremove(minFibNode*,int);											//arbitrary remove function implemented not used in Prim's algorithm
		 
		 int decreasekey(minFibNode*, int, int);										//Decrease Key function to run decrease the edge value of a vertex
		 
		 minFibNode* getMinPtr();														//get minimum fib heap pointer to 
		 
		 minFibNode* getMinFibHeap();
		 
		 minFibNode* getRefArrayElement(int);											//get Reference array element that is the node for a particular vertex index
		 
		 minFibNode* setRefArrayIndex(int);												//set Reference array element that is the node for a particular vertex index
};

// Following functions are explained above
void minFibHeap :: createAndInitRefArray(int num)
{
	refArray = new minFibNode*[num];
	
	//initialize refArray to null
	for(int i=0 ; i < num ; i++)
	{
		refArray[i] = NULL;
	}
}

// Returns pointer to the min Fibonacci heap
minFibNode* minFibHeap :: getMinFibHeap()
{
	return root;
}

minFibNode* minFibHeap :: setRefArrayIndex(int x)
{
	refArrayInd = x;
}

minFibNode* minFibHeap :: getMinPtr()
{
	return min;
}

minFibNode* minFibHeap :: getRefArrayElement(int n)
{
	return refArray[n];
}

minFibNode* minFibHeap :: insert(minFibNode* rNode,int val)
{
	minFibNode *currNode;

	//Adding element to an empty heap
	if(rNode->getValue() == -1)
	{
		// initializing heap and setting the fields
		currNode= new minFibNode;
		currNode->setLeft(currNode);
		currNode->setRight(currNode);		
		currNode->setChild(NULL);
		currNode->setParent(NULL);			
		currNode->setDegree(0);
		currNode->setValue(val);
		currNode->setChildCut(false);
		currNode->setIndex(refArrayInd);
		
		// Maintaining min pointer and root node the same
		rNode = currNode;
		min = currNode;

/*		// Check insert is not causing reference array bounds exceed.. If not add node to reference array
		if(refArrayInd < maxRefArrayInd)
		{
			refArray[refArrayInd] = rNode;
			refArrayInd++;
		}
		else
		{
			cout<<"Out of Reference Array Bounds!!\n\n";
		}*/
	}

	// If heap already contains more than zero elements
	else
	{
		// Allocate memory and set the fields		 
		currNode= new minFibNode;
		currNode->setLeft(NULL);
		currNode->setRight(NULL);		
		currNode->setChild(NULL);
		currNode->setParent(NULL);
		//cout << "Set degree for non root pointers" << endl;
		currNode->setDegree(0);
		currNode->setValue(val);
		currNode->setChildCut(false);
		currNode->setIndex(refArrayInd);		
		
		
		// Update min pointer if the value of newly added node is less than value pointed by min pointer
		if(currNode->getValue() < min->getValue())
		{
			min = currNode;
		}

		// Meld newly created node with heap
		rNode->getRight()->setLeft(currNode);
		currNode->setLeft(rNode);
		currNode->setRight(rNode->getRight());
		rNode->setRight(currNode);


		// Check insert is not causing reference array bounds exceed.. If not add node to reference array
		if(refArrayInd < maxRefArrayInd)
		{
			refArray[refArrayInd] = currNode;
			refArrayInd++;
		}
		else
		{
			cout<<"Out of Reference Array Bounds!!\n\n";
		}

		// Make pointer of heap point to minimum value node
		rNode = min;
	}
	
	root = rNode;	
	return rNode;	
}

// Pairwise Combine to improve amortized complexity
void minFibHeap :: pairwiseCombine(minFibNode* root1, minFibNode* root2)
{	//cout << "Entered Pairwise Combine" << endl;
	minFibNode* firstChild = NULL;
	minFibNode *retrieveNode, *temp;
	
	// Swap so that value of node 'root1' will always be smaller than value of 'root2'
	if(root1->getValue() > root2->getValue())
	{
		//swap root1 and root2
		temp = root1;
		root1 = root2;
		root2 = temp;
	}

	// Make root2 child of root1 if child of root1 is NULL
	if(root1->getChild() == NULL)
	{	//cout << "Set Degree for root1 pairwise combine root1->getIndex(): " << root1->getIndex() << "root2->getIndex(): " << root2->getIndex() << endl;
	//cout << "root1->getDegree(): " << root1->getDegree() << endl;
		root1->setDegree(1);
		root1->setChild(root2);
		root1->setChildCut(false);
		root2->setLeft(root2);
		root2->setRight(root2);			
		root2->setParent(root1);		
	}

	// if child of root1 is not NULL make root2 as one of the siblings of child of root1 by meld operation
	else if(root1->getChild() != NULL)
	{	
		firstChild = root1->getChild();

		//meld
		firstChild->getRight()->setLeft(root2);
		root2->setLeft(firstChild);
		root2->setRight(firstChild->getRight());
		firstChild->setRight(root2);

		// Increase degree of parent by 1
		//cout << "Set Degree for root1 added 1 pairwise combine root1->getDegree(): " << root1->getDegree() << " root1->getIndex():" << root1->getIndex() << endl;
		root1->setDegree(root1->getDegree() + 1);
		root2->setParent(root1);
		root1->setChildCut(false);

	}

	// if degree table slot is NULL, add merged heap in that empty slot
	int d;
	d = root1->getDegree();
	
	if(degT[d] == NULL)
	{
		degT[d] = root1;
	}

	// if degree table slot is not NULL, merge heap in the slot and current heap 

	else if(degT[d]!=NULL)
	{
		//cout << "pairwiseCombine for root1->getIndex(): " << root1->getIndex() << endl;	
		retrieveNode = degT[d];
		//if(retrieveNode != NULL)	
		//cout << "pairwiseCombine for retrieveNode->getIndex(): " << retrieveNode->getIndex() << endl;
		degT[d]=NULL;
		pairwiseCombine(root1,retrieveNode);
	}

}


int minFibHeap :: removeMin(minFibNode* rootNode)
{

	minFibNode *minptr,*minptrChild, *temp, *currentNode, *retrieveNode, *minptrNext,*ptr;
	int minIndex;
	
	//cout << "Entered Remove Min" << endl;
	ptr = new minFibNode;

	minptr = rootNode;
	// Temporarily store the node so that it can be freed later
	temp = rootNode;
	minIndex = minptr->getIndex();
	//cout << "rootNode Index: " << rootNode -> getIndex() << endl;
	// Delete the element from the reference array
	refArray[minIndex] = NULL;

	int firstRetrieveFromDegreeTable = 1;
	
	// Allocate memory for the degree table
	degT = new minFibNode*[500];

	//initialize the Degree Table to NULL
	for(int j = 0 ; j < 500 ; j++)
	{
		degT[j] = NULL;
	}	

	// Degree Table 
	if(minptr->getChild() != NULL)
	{
		// Update the degree table and make parents of children of min node NULL		
		int cDeg, cDeg1;
		cDeg = minptr->getChild()->getDegree();
		degT[cDeg] = minptr->getChild();
		minptr->getChild()->setParent(NULL);
		minptrChild = minptr->getChild()->getRight();
		while(minptrChild != minptr->getChild())
		{
			minFibNode* tempRight = minptrChild -> getRight();
			//cout << "Do setting for minpointer child" << endl;
			minptrChild->setParent(NULL);
			cDeg1 = minptrChild->getDegree();
			if(degT[cDeg1] == NULL)
			{
				degT[cDeg1] = minptrChild;
			}			
			else
			{
				retrieveNode = degT[cDeg1];
				degT[cDeg1] = NULL;
				pairwiseCombine(minptrChild,retrieveNode);
			}
			minptrChild = tempRight;		
		}	

	}		

	//for the first node
	//Except for a merge remaining 
	minptr = minptr->getRight();
	minptrNext = minptr->getRight();

	// Do pairwise combine and update degree table to improve amortized cost
	while(minptr != rootNode)
	{
		minptrNext = minptr->getRight();
		int deg1;
		deg1 = minptr->getDegree();
		if(degT[deg1] != NULL)
		{		
			retrieveNode = degT[deg1];
			degT[deg1] = NULL;
			pairwiseCombine(minptr,retrieveNode);
		}
		
		else if(degT[deg1] == NULL)
		{
			degT[deg1] = minptr;
		}

		minptr = minptrNext;
	}		

	// Collect min heaps from the degree table by traversing every slot of the degree table to form min fibonacci heap and update min pointer
	for(int degreeTableIndex = 0; degreeTableIndex < 500; degreeTableIndex++)
	{ 
		if(degT[degreeTableIndex] != NULL)
		{
			currentNode = degT[degreeTableIndex];

			if(firstRetrieveFromDegreeTable == 1)
			{
				ptr = currentNode; // ptr or rootNode?
				ptr->setLeft(ptr);
				ptr->setRight(ptr);
				min = ptr;
				firstRetrieveFromDegreeTable++;
			}
			else
			{
				//if(currentNode->val < min->val)
				if(currentNode->getValue() < min->getValue())
				{
					min = currentNode;
				}
				ptr->getRight()->setLeft(currentNode); // ptr or rootNode?
				currentNode->setLeft(ptr);
				currentNode->setRight(ptr->getRight());
				ptr->setRight(currentNode);
				ptr = min;
			}
		}
	}
	
	rootNode = ptr;
	delete(temp);
	root = rootNode;
	return minIndex;
}


// Remove arbitrary element from the min Fibonacci heap
int minFibHeap :: arbitraryremove(minFibNode* rNode,int ind)
{
	minFibNode *childptr, *nextChild, *temp, *currNode, *minptr, *nodeToArbRem, *cutSubTree;	
	int minIndex, degArbRem;	
	minptr = NULL;
	temp = NULL;
	
	// Temporarily store the node so that it can be freed later
	if(refArray[ind] != NULL)
	{						
		minptr = refArray[ind];
		temp = minptr;
	}

	if(minptr == NULL)
	{
		cout<<"\nNo element exist at the referred index!!\n";
	}

	else
	{
		nodeToArbRem = minptr;
		// Store the index of the node to be arbitrarily removed
		minIndex = nodeToArbRem->getIndex();	
	
		if(nodeToArbRem->getParent()!=NULL)
		{			
			// Decrease the degree of the parent of the deleted node by 1
			int degTemp;		
			degTemp = nodeToArbRem->getParent()->getDegree();
			nodeToArbRem->getParent()->setDegree(degTemp - 1);

			// if child pointer of the parent of node to be deleted points to node to be deleted, reset child pointer of parent
			if(nodeToArbRem->getParent()->getChild() == nodeToArbRem)
			{
				if(nodeToArbRem->getRight() == nodeToArbRem)
				{
					nodeToArbRem->getParent()->setChild(NULL);
				}
				else
				{
					nodeToArbRem->getParent()->setChild(nodeToArbRem->getRight());
				}
			}
		}
		
		// Remove node to be arbitrarily removed from the doubly linked list
		nodeToArbRem->getRight()->setLeft(nodeToArbRem->getLeft());
		nodeToArbRem->getLeft()->setRight(nodeToArbRem->getRight());

		// Store degree of node to be deleted
		degArbRem = nodeToArbRem->getDegree();
	
		if(degArbRem!=0)
		{
			//move to the child of node to be deleted
			childptr = nodeToArbRem->getChild();

			// Set parent of child of node to be deleted to NULL
			childptr->setParent(NULL);
			nextChild = childptr->getRight();

			// Meld all childern of the node to be deleted with all nodes at the top level
			for(int i = 0; i < degArbRem; i++)
			{
				nextChild = childptr->getRight();
				if(min->getValue() > childptr->getValue())
				{
					min = childptr;
				}

				rNode->getRight()->setLeft(childptr);
				childptr->setLeft(rNode);
				childptr->setRight(rNode->getRight());
				rNode->setRight(childptr);

				rNode = min;
				childptr = nextChild;
				childptr->setParent(NULL);
			}
		}

		// Do the extra work of setting childcut fot getting good amortized complexity
		if(nodeToArbRem->getParent()!=NULL)
		{					
			currNode = nodeToArbRem;

			// Keep going up till we either hit node with childCut set False or root node
			while(currNode->getParent()->getChildCut() != false)
			{
				// Cut the node if its childcut is True for melding subtree rooted at it with the top level nodes
				cutSubTree = currNode->getParent();
				//cout << "Degree cutSubTree->getIndex(): " << cutSubTree->getIndex() << endl;
				

				int degPar;
				//cout << "Degree cutSubTree->getParent()->getIndex(): " << cutSubTree->getParent()->getIndex() << endl;
				degPar = cutSubTree->getParent()->getDegree(); 
				// Go on reducing degree of parents up the heap
				
				cutSubTree->getParent()->setDegree(degPar - 1);
				//cout << "Degree cutSubTree->getParent()->getDegree(): " << cutSubTree->getParent()->getDegree() << endl;

				// Seperate subtree rooted at cutSubtree
				if(cutSubTree->getParent()->getChild() == cutSubTree)
				{
					// if only one node
					if(cutSubTree->getRight() == cutSubTree)
					{
						cutSubTree->getParent()->setChild(NULL);
					}
					// if multiple nodes
					else
					{
						cutSubTree->getParent()->setChild(cutSubTree->getRight());
					}
				}
				
				cutSubTree->getRight()->setLeft(cutSubTree->getLeft());
				cutSubTree->getLeft()->setRight(cutSubTree->getRight());


				// Reset min pointer if root of the removed subtree has value less than that pointed by min pointer
				if(cutSubTree->getValue() < min->getValue())
				{
					min = cutSubTree;
				}
	
				// Add separated subtree next to the min node at the top level
				rNode->getRight()->setLeft(cutSubTree);
				cutSubTree->setLeft(rNode);
				cutSubTree->setRight(rNode->getRight());
				rNode->setRight(cutSubTree);
				rNode = min;

				// Move upwards to check the childcut value
				currNode = currNode->getParent();
			}

			// as soon as we encounter node with child cut value False, set the the childcut value of the node to true 
			currNode->getParent()->setChildCut(true);
		}

	}

	// Remove the element from the reference array
	int tempInd = temp->getIndex();
	if(temp!=NULL){
		refArray[tempInd] = NULL;
	}

	
	//Physically delete the node
	delete(temp);	
	root = rNode;
	return temp->getIndex();
}

// Decrease key of arbitrary node by the input amount
int minFibHeap :: decreasekey(minFibNode* rNode, int ind, int amtToDec)
{
	
	minFibNode *currNode, *minptr, *curr, *temp , *cutSubTree, *nodeToDecKey;
	minptr = NULL;
	bool tempChildCut;
	
	// Temporarily store the node whose value is to be decreased
	//cout << "Entered Dec Key" << endl;
	//cout << " min->getValue(): " << min->getValue() << endl;
	//cout << "Index to dec key: " << ind << endl;
		//cout << "Dec Key for min->getIndex(): " << min->getIndex() << " min->getValue(): " << min->getValue()  << "Dec Key for min->getRight()->getIndex(): " << min->getRight()->getIndex() << endl;
		
		//	if(min->getRight()->getChild()!=NULL)
	//cout << "Dec Key for min->getRight()->getChild()->getIndex(): " << min->getRight()->getChild()->getIndex() << endl;
	//if(min->getChild()!=NULL)
	//cout << "Dec Key for min->getChild()->getIndex(): " << min->getChild()->getIndex() << endl;
	if(refArray[ind] != NULL)
	{
		minptr = refArray[ind];
	}
	
	
	if(minptr == NULL)
	{
		//cout << "MinPtr is Null in Decrease Key" << endl;
		return -1;
	}
	else
	{
		//cout << "Dec Key for minptr->getIndex(): " << minptr->getIndex() << endl;
		//cout << "Dec Key for minptr->getDegree(): " << minptr->getDegree() << endl;
		nodeToDecKey = minptr;
		// Decrease key
		nodeToDecKey->setValue(nodeToDecKey->getValue() - amtToDec);
		
		
		// If found Node is not root of some min heap
		if(nodeToDecKey->getParent() != NULL)
		{		
			// If value of key is smaller than the parent after decreasing key	
			if(nodeToDecKey->getValue() < nodeToDecKey->getParent()->getValue())
			{
				// Reduce degree of parent as we are cutting subtree rooted at node whose key is decreased
				//"Reduce Degree of Parent " << ind << endl;
				int degPar;		
				degPar = nodeToDecKey->getParent()->getDegree();
				//cout << "Node to Dec Key nodeToDecKey->getIndex(): " << nodeToDecKey->getIndex() << endl;
				//cout << "Node to Dec Key parent nodeToDecKey->getParent()->getIndex(): " << nodeToDecKey->getParent()->getIndex() << endl;
				nodeToDecKey->getParent()->setDegree(degPar - 1);
				//cout << "Node to Dec Key parent nodeToDecKey->getParent()->getDegree(): " << nodeToDecKey->getParent()->getDegree() << endl;

				// Reset appropriate child pointer of the parent and separate the subtree
				if(nodeToDecKey->getParent()->getChild() == nodeToDecKey)
				{
					if(nodeToDecKey->getRight() == nodeToDecKey)
					{
						nodeToDecKey->getParent()->setChild(NULL);
					}
					else
					{
						nodeToDecKey->getParent()->setChild(nodeToDecKey->getRight());
						tempChildCut = nodeToDecKey->getParent()->getChildCut();
					}
				}
				nodeToDecKey->getRight()->setLeft(nodeToDecKey->getLeft());
				nodeToDecKey->getLeft()->setRight(nodeToDecKey->getRight());

				// Reset min pointer if required before melding separated subtree at the top level 
				currNode = nodeToDecKey;
				if(currNode->getValue() < min->getValue())
				{
					min = currNode;
				}
	
				// Meld separated subtree at the top level by adding it next to node pointed by min Fibonacci Heap Pointer and update min Fibonacci  					// Heap Pointer poining to min element 
				rNode->getRight()->setLeft(currNode);
				currNode->setLeft(rNode);
				currNode->setRight(rNode->getRight());
				rNode->setRight(currNode);
				rNode = min;

				curr = currNode;

				// Keep going up till we either hit node with childCut set False or root node
				if(tempChildCut){
				while(curr->getParent()!=NULL && curr->getParent()->getChildCut() != false)
				{
					
					minFibNode *tempPar = curr->getParent();
					// Cut the node if its childcut is True for melding subtree rooted at it with the top level nodes
					cutSubTree = curr->getParent();
					
					
					// Go on reducing degree of parents up the heap
					int parDeg;
					if(cutSubTree->getParent() != NULL){
						parDeg = cutSubTree->getParent()->getDegree();
						//cout << "Loop2" << endl;
				//cout << "Node to Dec Key cutSubTree->getIndex(): " << cutSubTree->getIndex() << endl;
				//cout << "Node to Dec Key parent cutSubTree->getParent()->getIndex(): " << cutSubTree->getParent()->getIndex() << endl;
				
						cutSubTree->getParent()->setDegree(parDeg - 1);
					//	cout << "Node to Dec Key parent cutSubTree->getParent()->getDegree(): " << cutSubTree->getParent()->getDegree() << endl;
				
						// Seperate subtree rooted at cutSubtree
						if(cutSubTree->getParent()->getChild() == cutSubTree)
						{
							// if only one node
							if(cutSubTree->getRight() == cutSubTree)
							{
								cutSubTree->getParent()->setChild(NULL);
							}
							// if multiple nodes
							else
							{
								cutSubTree->getParent()->setChild(cutSubTree->getRight());
							}
						}
					}
					cutSubTree->getRight()->setLeft(cutSubTree->getLeft());
					cutSubTree->getLeft()->setRight(cutSubTree->getRight());
					
					// Reset min pointer if required before melding separated subtree at the top level 
					if(cutSubTree->getValue() < min->getValue())
					{
						min = cutSubTree;
					}

				// Meld separated subtree at the top level by adding it next to node pointed by min Fibonacci Heap Pointer and update min Fibonacci  					// Heap Pointer poining to min element
				    rNode->getRight()->setLeft(cutSubTree);
					cutSubTree->setLeft(rNode);
					cutSubTree->setRight(rNode->getRight());
					rNode->setRight(cutSubTree);
					rNode = min;
					
					curr->setParent(NULL);
					curr = tempPar;
				}
				// as soon as we encounter node with child cut value False, set the the childcut value of the node to true 
					curr->setChildCut(true);

			}
			nodeToDecKey->setParent(NULL);
		}
		
		}
		//If node to decrease key is root of some min heap update min pointer if required
		else
		{
			if(nodeToDecKey->getValue() < min->getValue())
			{
				min = nodeToDecKey;
				rNode = min;
			}

		}
	}
	
	//cout << "Dec Key for min->getIndex(): " << min->getIndex() << " min->getValue(): " << min->getValue()  << endl;
	//cout << "Dec Key for rNode->getIndex(): " << rNode->getIndex() << " rNode->getValue(): " << rNode->getValue()  << endl;
	
	temp = rNode;
	temp->setParent(NULL);

	root = rNode;	
	//cout << "End of Dec key" << endl;
	return 1;
}


int fibMST(CGraph G, int *frmVerInd){
	//cout << "Entered fibMST" << endl;
	minFibHeap fHeap (G.numVer);
	int tempInd;
	int edgeValSum = 0;
	hash_map<int,int>::iterator it;
	minFibNode *tempChild;
	int count;
	for(int i = 0; i < G.numVer; i++){
		fHeap.insert(fHeap.getMinPtr(),G.Vertices[i].vVal);
		//cout << "Index Value of " << i << ": " << fHeap.getRefArrayElement(i)->getValue() << endl;
	}
	for(int verIndRem = 0; verIndRem < G.numVer; verIndRem++){
		if(fHeap.getMinPtr() != NULL){
			tempInd = fHeap.removeMin(fHeap.getMinFibHeap());
			//cout << "Removed Index: " << tempInd << endl;
			//cout << "Verindem: " << verIndRem << endl;
		}
	/*	for(int j = 0; j < G.numVer; j++){
			if(fHeap.getRefArrayElement(j) != NULL){
			if(fHeap.getRefArrayElement(j)->getRight() != NULL)
			cout << "Right Value of " << j << ": " << fHeap.getRefArrayElement(j)->getRight()->getIndex() << endl;
			if(fHeap.getRefArrayElement(j)->getLeft() != NULL)
			cout << "Left Value of " << j << ": " << fHeap.getRefArrayElement(j)->getLeft()->getIndex() << endl;
			if(fHeap.getRefArrayElement(j)->getParent() != NULL)
			cout << "Parent Value of " << j << ": " << fHeap.getRefArrayElement(j)->getParent()->getIndex() << endl;
			if(fHeap.getRefArrayElement(j)->getChild() != NULL)
			cout << "Child Value of " << j << ": " << fHeap.getRefArrayElement(j)->getChild()->getIndex() << endl;
			cout << "getDegree Value of " << j << ": " << fHeap.getRefArrayElement(j)->getDegree() << endl;
			cout << "getValue Value of " << j << ": " << fHeap.getRefArrayElement(j)->getValue() << endl;
			cout << endl;
		}
	}*/
		if(verIndRem != 0){
            tempList.push_back(tempInd);
		    tempList.push_back(frmVerInd[tempInd]);
		    tempList.push_back(G.Vertices[tempInd].vVal);			
		    edgeValSum = edgeValSum + G.Vertices[tempInd].vVal;
	}
		if(verIndRem == G.numVer - 1){
			break;
			//cout << "Ended Here" << endl;
		}
		it = G.edgeMap[tempInd].begin();
		while(it != G.edgeMap[tempInd].end()){
			if(it->second < G.Vertices[it->first].vVal){
				int decVal = G.Vertices[it->first].vVal - (it->second);
				G.Vertices[it->first].vVal = it->second;
				if(fHeap.getMinFibHeap())
					fHeap.decreasekey(fHeap.getMinPtr(),it->first,decVal);
				frmVerInd[it->first] = tempInd;	
			}
			++it;
		}
		count = verIndRem;
	/*	for(int j = 0; j < G.numVer; j++){
			if(fHeap.getRefArrayElement(j) != NULL){
			if(fHeap.getRefArrayElement(j)->getRight() != NULL)
			cout << "Right Value of " << j << ": " << fHeap.getRefArrayElement(j)->getRight()->getIndex() << endl;
			if(fHeap.getRefArrayElement(j)->getLeft() != NULL)
			cout << "Left Value of " << j << ": " << fHeap.getRefArrayElement(j)->getLeft()->getIndex() << endl;
			if(fHeap.getRefArrayElement(j)->getParent() != NULL)
			cout << "Parent Value of " << j << ": " << fHeap.getRefArrayElement(j)->getParent()->getIndex() << endl;
			if(fHeap.getRefArrayElement(j)->getChild() != NULL)
			cout << "Child Value of " << j << ": " << fHeap.getRefArrayElement(j)->getChild()->getIndex() << endl;
			cout << "getDegree Value of " << j << ": " << fHeap.getRefArrayElement(j)->getDegree() << endl;
			cout << "getValue Value of " << j << ": " << fHeap.getRefArrayElement(j)->getValue() << endl;
			cout << endl;
		}
	}*/	
	}
	return edgeValSum;
}


//Main method
int main(int argc, const char* argv[])
{
	//input from command line
	//clock_t Start, Time;
	CGraph G;
	if (argc == 3 || argc == 4)
	{
		if (tolower(*(argv[1] + 1)) == 'r')// || tolower(*(argv[1] + 1)) == 'f')
		{
			//cout << "Entered here" << endl;
			//random(100);
			int n = atoi(argv[2]);
			int d = atoi(argv[3]);
			int maxEdgeCount = (n*(n-1))/2;
			if(maxEdgeCount*d/100 < n-1){
				cout << "Not enough edges to form a MST for the given vertices" << endl;
				return -1;
			}
			G = createGraph(n, maxEdgeCount, d);
			//cout << "Hi" << endl;
			//Start = clock();
			primMst(G,'r');
			//primMst(G,'f');
			//	Time = clock() - Start;	
		//cout << endl;
		//cout << Time << endl;
		}
		else if (*(argv[1] + 1) == 's' || *(argv[1] + 1) == 'f'){
			ifstream fs (argv[2]);
			if(fs.is_open()){
				int t1, t2, t3, uN, uM, iter = 0;					
				fs >> uN >> uM;
			//	cout << uN << " " << uM << endl;
				int uEdge[uM];
				int uVer1Ind[uM], uVer2Ind[uM];
				while(fs >> t1 >> t2 >> t3){
					//cout << t3 << " " << t2 << " " << t1 << endl;
					uEdge[iter] = t3;
					uVer1Ind[iter] = t1;
					uVer2Ind[iter] = t2;
					//cout << t3 << " " << t2 << " " << t1 << endl;
					iter++;
				}	
			
			fs.close();
			G = createUserGraph(uN, uM, uVer1Ind, uVer2Ind, uEdge);
			if  (*(argv[1] + 1) == 'f')
				primMst(G, 'f');
			else if  (*(argv[1] + 1) == 's')
				primMst(G, 's');				
			}
			else{
				cout << "Could not open the file" << endl;
				return -1;
			}
		}
		else
			cout << "Wrong Argument Sequence" << endl;
			return -1;
	}
	else{
		cout << "Wrong Argument Sequence" << endl;
		return -1;	
	}		
	return 0;				
}

