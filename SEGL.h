#pragma once
#ifndef SEGL_H
#define SEGL_H

// SEGL is for Simple Efficient Graph Library
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
using namespace std;

template <class T> 
struct g_node {
	T value;
	set<pair<int, int>> adjList; // first is the node index, and second is the weight
	g_node() : value() {}
	g_node(const T& data) : value(data) {}
};

template <class T>
class SEGL {
private:
	int numNodes;
	int numEdges;
	vector<g_node<T>> nodes;
	void dfs(int node, vector<bool>& visited);
public:
	SEGL();
	SEGL(int size);

	void addNode();
	void addNode(const T& value);
	
	void addEdge(int node_1, int node_2);
	void addDirEdge(int from, int to);
	void addWeightedEdge(int node_1, int node_2, int weight);
	void addDirWeightedEdge(int from, int to, int weight);
	
	void removeEdge(int node_1, int node_2);
	void removeEdge(int node_1, int node_2, int weight);
	void removeDirEdge(int from, int to);
	void removeDirEdge(int from, int to, int weight);
	
	void setNodeValue(int node, const T& value);
	void setNodeValue(const T& old, const T& New); 
	T getNodeValue(int node);

	set<pair<int, int>> adjacent(int node);

	int size();

	bool existsNode(const T& value); 

	bool existsEdge(int from, int to); // if exists directed edge, then exists undirected edge
	bool existsEdge(int from, int to, int weight);

	int edgeWeight(int node_1, int node_2); // two nodes, return weight

	void transpose();

	int numConnectedComponents(); // TO-DO: work on direct access instead of implementing dfs each time

	bool existsCycle();

	bool isTree();

	bool existsPath(int from, int to);

	int shortestDistance(int from, int to); 

	pair<int, vector<int>> shortestPath(); // TO-DO: checking if there exists a path or no
	// To Do: understand the bellman-ford algorithm
	
	// traversal 
	
	bool inSameConnectedComponent();
	bool existsHamiltonian(); // *
	bool existsEulerian(); // *
	
	void condenseSCCs(); // *
	// condense strongly connected components (Kosaraju's Algorithm)
	
	vector<int> topoSort(); // *
	
	void clear();
	~SEGL();
};

template <class T>
SEGL<T>::SEGL() 
	: numNodes(0), numEdges(0)
{}

template <class T>
SEGL<T>::SEGL(int size) : numNodes(size), numEdges(0) {
	nodes.resize(size);
}

template <class T>
void SEGL<T>::addNode() {
	numNodes++;
	nodes.push_back(g_node<T>());
}

template <class T>
void SEGL<T>::addNode(const T& value) {
	numNodes++;
	nodes.push_back(g_node<T>(value));
}

template <class T>
void SEGL<T>::addDirEdge(int from, int to) {
	nodes[from].insert({ to, 1 });
	numEdges++;
}

template <class T> 
void SEGL<T>::addEdge(int node_1, int node_2) {
	addDirEdge(node_1, node_2);
	addDirEdge(node_2, node_1);
}

template <class T>
void SEGL<T>::addDirWeightedEdge(int from, int to, int weight) {
	nodes[from].insert({ to, weight });
	numEdges++;
}

template <class T> 
void SEGL<T>::addWeightedEdge(int node_1, int node_2, int weight) {
	addDirWeightedEdge(node_1, node_2, weight);
	addDirWeightedEdge(node_2, node_1, weight);
}

template <class T>
void SEGL<T>::removeDirEdge(int from, int to) {
	auto edge = lower_bound(nodes[from].adjList.begin(), nodes[from].adjList.end(), make_pair(to, INT_MIN));
	if (edge->first != to) return;
	auto it = edge++;
	if (edge == nodes[from].adjList.end()) 
		nodes[from].adjList.erase(it);
	else {
		while (edge != nodes[from].adjList.end() && it->first == to) {
			nodes[from].adjList.erase(it);
			it = edge++;
		}
		if (it != nodes[from].adjList.end() && it->first == to) nodes[from].adjList.erase(it);
	}
}

template <class T> 
void SEGL<T>::removeEdge(int node_1, int node_2) { 
	// without weight given, all multiple edges from node 1 to node 2 will be removed	
	removeDirEdge(node_1, node_2);
	removeDirEdge(node_2, node_1);
}

template <class T>
void SEGL<T>::removeDirEdge(int from, int to, int weight) { 
	// differs if there are multiple edges of the same weights
	auto edge = lower_bound(nodes[from].adjList.begin(), nodes[from].adjList.end(), make_pair(to, weight));
	if (edge->first != to || edge->second != weight) return;
	auto it = edge++;
	if (edge == nodes[from].adjList.end())
		nodes[from].adjList.erase(it);
	else {
		while (edge != nodes[from].adjList.end() && it->first == to && it->second == weight) {
			nodes[from].adjList.erase(it);
			it = edge++;
		}
		if (it != nodes[from].adjList.end() && it->first == to && it->second == weight) 
			nodes[from].adjList.erase(it);
	}
}

template <class T>
void SEGL<T>::removeEdge(int node_1, int node_2, int weight) {
	removeDirEdge(node_1, node_2, weight);
	removeDirEdge(node_2, node_1, weight);
}

template <class T>
bool SEGL<T>::existsEdge(int from, int to) {
	auto edge = lower_bound(nodes[from].adjList.begin(), nodes[from].adjList.end(), make_pair(to, INT_MIN));
	return (edge != nodes[from].adjList.end() && edge->first == to);
}

template <class T>
bool SEGL<T>::existsEdge(int from, int to, int weight) {
	return nodes[from].adjList.count({ to, weight });
}

template <class T>
bool SEGL<T>::existsNode(const T& value) {
	bool flag = false;
	for (int i = 0; i < nodes.size(); i++)
		if (nodes[i].value == value) {
			flag = true;
			break;
		}
	return flag;
}

template <class T>
int SEGL<T>::edgeWeight(int from, int to) { // returns minimum edge weight if exists, otherwise infinity
	auto edge = lower_bound(nodes[from].adjList.begin(), nodes[from].adjList.end(), make_pair(to, INT_MIN));
	if (edge == nodes[from].adjList.end() || edge->first != to) return INT_MAX;
}
template <class T>
void SEGL<T>::setNodeValue(int node, const T& value) {
	nodes[node].value = value;
}

template <class T>
void SEGL<T>::setNodeValue(const T& old, const T& New) {
	for (int i = 0; i < nodes.size(); i++)
		if (nodes[i].value == old) 
			nodes[i].value = New;
}

template <class T>
T SEGL<T>::getNodeValue(int node) {
	retrun nodes[node].value;
}

template <class T>
set<pair<int, int>> SEGL<T>::adjacent(int node) {
	return nodes[node].adjList;
}

template <class T>
int SEGL<T>::size() {
	return numNodes;
}

template <class T>
void SEGL<T>::transpose() {
	vector<set<pair<int, int>>> adjListnew(numNodes);
	// create a new adjacency list with the transpose
	for (int i = 0; i < numNodes; i++) 
		for (auto u : nodes[i].adjList) 
			adjListnew[u.first].insert(make_pair(i, u.second));
	for (int i = 0; i < numNodes; i++)
		nodes[i].adjList.clear();
	for (int i = 0; i < numNodes; i++)
		nodes[i].adjList = adjListnew[i];
}

template <class T>
int numConnectedComponents() { // DFS
	vector<bool> visited(numNodes);
	int count = 0;
	for (int i = 0; i < numNodes; i++)
		if (!visited[i]) {
			dfs(i, visited);
			count++;
		}
	return count;
}

template <class T>
void SEGL<T>::dfs(int node, vector<bool>& visited) {
	visited[node] = 1;
	for (auto u : nodes[node].adjList)
		if (!visited[u.f]) dfs(u.f);
}

template <class T>
bool SEGL<T>::existsCycle() {
	return (numEdges > numNodes - numConnectedComponents());
}

template <class T>
bool SEGL<T>::isTree() {
	return (numConnectedComponents() == 1 && numEdges == numNodes - 1);
}



#endif // !SEGL_H