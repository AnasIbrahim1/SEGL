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
	g_node(T data) : value(data) {}
};

template <class T>
class SEGL {
private:
	int numNodes;
	vector<g_node<T>> nodes;
public:
	SEGL();
	SEGL(int size);
	void addNode(T value);
	void addNode();
	void addEdge(int node_1, int node_2);
	void addDirEdge(int from, int to);
	void addWeightedEdge(int node_1, int node_2, int weight);
	void addDirWeightedEdge(int from, int to, int weight);
	void removeEdge(int node_1, int node_2);
	void removeEdge(int node_1, int node_2, int weight);
	void removeDirEdge(int from, int to);
	void removeDirEdge(int from, int to, int weight);
	vector<int> shortestPath(); // TO-DO: checking if there exists a path or not
	double shortestDistance(); // two input
	// To Do: understand the bellman-ford algorithm
	bool existsEdge();
	bool findNode(); //  
	int edgeWeight(); // two nodes, return weight
	// traversal 
	bool existsCycle();
	int numConnectedComponents(); // TO-DO: work on direct access instead of implementing dfs each time
	vector<int> findCycle();
	bool inSameConnectedComponent();
	bool existsHamiltonian();
	bool existsEulerian();
	int size();
	void condenseSCCs(); // condense strongly connected components (Kosaraju's Algorithm)
	void transpose();
	vector<int> topoSort();
	void clear();
	~SEGL();
};

template <class T>
SEGL<T>::SEGL() 
	: numNodes(0)
{}

template <class T>
SEGL<T>::SEGL(int size) : numNodes(size) {
	nodes.resize(size);
}

template <class T>
void SEGL<T>::addNode(T value) {
	numNodes++;
	nodes.push_back(g_node<T>(value));
}

template <class T>
void SEGL<T>::addDirEdge(int from, int to) {
	nodes[from].insert({ to, 1 });
}

template <class T> 
void SEGL<T>::addEdge(int node_1, int node_2) {
	addDirEdge(node_1, node_2);
	addDirEdge(node_2, node_1);
}

template <class T>
void SEGL<T>::addDirWeightedEdge(int from, int to, int weight) {
	nodes[from].insert({ to, weight });
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

#endif // !SEGL_H