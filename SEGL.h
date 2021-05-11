#pragma once
#ifndef SEGL_H
#define SEGL_H

// SEGL is for Simple Efficient Graph Library
#include <iostream>
#include <unordered_set>
#include <set>
#include <queue>
#include <vector>
#include <algorithm>
using namespace std;

template <class T>
struct g_node {
	T value;
	set<pair<int, int>> adjList; // first is the node index, and second is the weight
	g_node() : value() {}
	g_node(const T& data) : value(data) {}
	g_node(const g_node& old_node) {
		this->value = old_node.value;
		this->adjList.clear(); this->adjList = old_node.adjList;
	}
	void operator = (const g_node& old_node) {
		this->clear();
		this->value = old_node.value;
		this->adjList.clear(); this->adjList = old_node.adjList;
	}
	void clear() {
		adjList.clear();
	}
	~g_node() {
		this->clear();
	}
};

template <class T>
class SEGL {
private:
	int numNodes;
	int numEdges;
	vector<g_node<T>> nodes;
	vector<g_node<T>> nodes_transpose;
	void dfs(int node, vector<bool>& visited);
	bool existsPath(int from, int to, vector<bool>& visited);
	bool containNegativeWeights;
	unordered_set<int> hasOutdegree; // store the indices of the nodes that have such that outdegree != 0
	unordered_set<int> hasOutdegreeTranspose; // same but for transpose
	vector<pair<pair<int, int>, int>> listEdges();
	bool topSort(vector<int>& sorted, SEGL<T>& new_graph, vector<bool>& finished, int& finishedCount);

public:
	SEGL();
	SEGL(int size);

	SEGL(const SEGL& old_graph);
	void operator = (const SEGL& new_graph);

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
	void setNodeValue(const T& old, const T& New); // O(V) 
	T getNodeValue(int node);

	int indegree(int node);
	int outdegree(int node);

	set<pair<int, int>> adjacent(int node);

	vector<pair<pair<int, int>, int>> getEdgeList();

	int numOfNodes();
	int numOfDirEdges();
	int numOfEdges(); // assuming the graph is undirected

	bool existsNode(const T& value);

	bool existsEdge(int from, int to); // if exists directed edge, then exists undirected edge
	bool existsEdge(int from, int to, int weight);

	int edgeWeight(int node_1, int node_2);

	void transpose();

	int numConnectedComponents(); // assuming undirected
	// TO-DO: work on direct access instead of implementing dfs each time

	bool existsCycle(); // assuming undirected

	bool isTree(); // assuming undirected

	bool existsPath(int from, int to); // Normal DFS, to do is optimize

	int shortestDistance(int from, int to);

	vector<int> shortestDistanceAll(int from);

	pair<int, vector<int>> shortestPath(int from, int to);

	pair<bool, vector<int>> topSort(); // first is -1 if there doesn't exist a topological Sorting

	bool isDAG();

	//int numSCCs(); // Kosaraju's

	// traversal 
	//bool existsHamiltonianPath(); // *
	bool existsEulerianPath();

	//void condenseSCCs(); // *
	// condense strongly connected components (Kosaraju's Algorithm)

	void clear();
	~SEGL();
};

template <class T>
SEGL<T>::SEGL()
	: numNodes(0), numEdges(0)
{
	containNegativeWeights = false;
}

template <class T>
SEGL<T>::SEGL(int size) : numNodes(size), numEdges(0) {
	nodes.resize(size);
	nodes_transpose.resize(size);
	containNegativeWeights = false;
}

template <class T>
SEGL<T>::SEGL(const SEGL& old_graph) {
	this->clear();
	this->numNodes = old_graph.numNodes;
	this->numEdges = old_graph.numEdges;
	this->containNegativeWeights = old_graph.containNegativeWeights;
	this->nodes = old_graph.nodes;
	this->nodes_transpose = old_graph.nodes_transpose;
}

template <class T>
void SEGL<T>::operator = (const SEGL& old_graph) {
	this->clear();
	this->numNodes = old_graph.numNodes;
	this->numEdges = old_graph.numEdges;
	this->containNegativeWeights = old_graph.containNegativeWeights;
	this->nodes = old_graph.nodes;
	this->nodes_transpose = old_graph.nodes_transpose;
}

template <class T>
void SEGL<T>::addNode() {
	numNodes++;
	nodes.push_back(g_node<T>());
	nodes_transpose.push_back(g_node<T>());
}

template <class T>
void SEGL<T>::addNode(const T& value) {
	numNodes++;
	nodes.push_back(g_node<T>(value));
	nodes_transpose.push_back(g_node<T>(value));
}

template <class T>
void SEGL<T>::addDirEdge(int from, int to) {
	nodes[from].adjList.insert({ to, 1 });
	nodes_transpose[to].adjList.insert({ from, 1 });
	hasOutdegree.insert(from);
	hasOutdegreeTranspose.insert(to);
	numEdges++;
}

template <class T>
void SEGL<T>::addEdge(int node_1, int node_2) {
	addDirEdge(node_1, node_2);
	addDirEdge(node_2, node_1);
}

template <class T>
void SEGL<T>::addDirWeightedEdge(int from, int to, int weight) {
	if (weight < 0) containNegativeWeights = true;
	nodes[from].insert({ to, weight });
	nodes_transpose[to].insert({ from, weight });
	hasOutdegree.insert(from);
	hasOutdegreeTranspose.insert(to);
	numEdges++;
}

template <class T>
void SEGL<T>::addWeightedEdge(int node_1, int node_2, int weight) {
	addDirWeightedEdge(node_1, node_2, weight);
	addDirWeightedEdge(node_2, node_1, weight);
}

template <class T>
void SEGL<T>::removeDirEdge(int from, int to) {
	// erasing from normal graph
	auto edge = lower_bound(nodes[from].adjList.begin(), nodes[from].adjList.end(), make_pair(to, INT_MIN));
	if (edge == nodes[from].adjList.end() || edge->first != to) return;
	auto it = edge++;
	if (edge == nodes[from].adjList.end()) {
		nodes[from].adjList.erase(it);
	}
	else {
		while (edge != nodes[from].adjList.end() && it->first == to) {
			nodes[from].adjList.erase(it);
			it = edge++;
		}
		if (it != nodes[from].adjList.end() && it->first == to) nodes[from].adjList.erase(it);
	}
	if (nodes[from].adjList.size() == 0) hasOutdegree.erase(from);

	// erasing from transpose
	auto edge2 = lower_bound(nodes_transpose[to].adjList.begin(), nodes_transpose[to].adjList.end(), make_pair(from, INT_MIN));
	if (edge2 == nodes_transpose[to].adjList.end() || edge2->first != from) return;
	auto it2 = edge2++;
	if (edge2 == nodes_transpose[to].adjList.end())
		nodes_transpose[to].adjList.erase(it2);
	else {
		while (edge2 != nodes_transpose[to].adjList.end() && it->first == from) {
			nodes_transpose[to].adjList.erase(it);
			it2 = edge2++;
		}
		if (it2 != nodes_transpose[to].adjList.end() && it2->first == from) nodes_transpose[to].adjList.erase(it2);
	}
	if (nodes_transpose[to].adjList.size() == 0) hasOutdegreeTranspose.erase(to);
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
	if (nodes[from].adjList.size() == 0) hasOutdegree.erase(from);

	// transpose
	auto edge2 = lower_bound(nodes_transpose[to].adjList.begin(), nodes_transpose[to].adjList.end(), make_pair(from, weight));
	if (edge2->first != from || edge2->second != weight) return;
	auto it2 = edge2++;
	if (edge2 == nodes_transpose[to].adjList.end())
		nodes_transpose[to].adjList.erase(it2);
	else {
		while (edge2 != nodes_transpose[to].adjList.end() && it2->first == from && it2->second == weight) {
			nodes_transpose[to].adjList.erase(it2);
			it2 = edge2++;
		}
		if (it2 != nodes_transpose[to].adjList.end() && it2->first == from && it2->second == weight)
			nodes_transpose[to].adjList.erase(it2);
	}
	if (nodes_transpose[to].adjList.size() == 0) hasOutdegreeTranspose.erase(to);
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
int SEGL<T>::outdegree(int node) {
	return nodes[node].adjList.size();
}

template <class T>
int SEGL<T>::indegree(int node) {
	return nodes_transpose[node].adjList.size();
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
	return nodes[node].value;
}

template <class T>
set<pair<int, int>> SEGL<T>::adjacent(int node) {
	return nodes[node].adjList;
}

template <class T>
int SEGL<T>::numOfNodes() {
	return numNodes;
}

template <class T>
int SEGL<T>::numOfDirEdges() {
	return numEdges;
}

template <class T>
int SEGL<T>::numOfEdges() {
	return (numEdges / 2);
}

template <class T>
vector<pair<pair<int, int>, int>> SEGL<T>::listEdges() { // O(E)
	vector<pair<pair<int, int>, int>> edgeList;
	for (auto it = hasOutdegree.begin(); it != hasOutdegree.end(); it++)
		for (auto u : nodes[*it].adjList)
			edgeList.push_back({ {*it, u.first}, u.second });
	return edgeList;
}

template <class T>
void SEGL<T>::transpose() {
	swap(nodes, nodes_transpose);
	swap(hasOutdegree, hasOutdegreeTranspose);
}

template <class T>
int SEGL<T>::numConnectedComponents() { // DFS
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
bool SEGL<T>::existsCycle() { // assumes that the graph is undirected
	return (numEdges > numNodes - numConnectedComponents());
}

template <class T>
bool SEGL<T>::isTree() { // function assumes that the graph is undirected
	return (numConnectedComponents() == 1 && numEdges == numNodes - 1);
}

template <class T>
bool SEGL<T>::existsPath(int from, int to) {
	vector<bool> visited(numNodes);
	return existsPath(from, to, visited);
}

template <class T>
bool SEGL<T>::existsPath(int from, int to, vector<bool>& visited) {
	visited[from] = 1;
	if (from == to) return true;
	bool flag = false;
	for (auto u : nodes[from].adjList)
		if (!visited[u.f])
			flag = flag || existsPath(u.f, to, visited);
	return flag;
}

template <class T>
void SEGL<T>::clear() {
	numNodes = 0;
	numEdges = 0;
	nodes.clear();
	containNegativeWeights = false;
}

template <class T>
SEGL<T>::~SEGL() {
	this->clear();
}

template <class T>
vector<int> SEGL<T>::shortestDistanceAll(int from) {
	if (containNegativeWeights) { // Bellman-Ford O(V*E)
		vector<int> dist(numNodes, INT_MAX);
		dist[from] = 0;
		auto edgeList = listEdges();
		for (int i = 0; i < numNodes - 1; i++)
			for (auto edge : edgeList)
				dist[edge.first.second] = min(dist[edge.first.second], dist[edge.first.first] + edge.second);
		return dist;
	}
	else { // Dijkstra's, O(E + VLOG(V))
		vector<bool> visited(numNodes);
		vector<int> dist(numNodes, INT_MAX);
		dist[from] = 0;
		priority_queue<pair<int, int>> q;
		q.push({ 0, from });
		while (!q.empty()) {
			int a = q.top().second; q.pop();
			visited[a] = true;
			for (auto u : nodes[a].adjList)
				if (!visited[u.f]) {
					dist[u.f] = min(dist[u.f], dist[a] + u.s);
					q.push({ -1 * dist[u.f], u.f });
				}
		}
		return dist;
	}
}

template <class T>
int SEGL<T>::shortestDistance(int from, int to) {
	vector<int> dist = shortestDistanceAll(from);
	return dist[to];
}

template <class T>
pair<int, vector<int>> SEGL<T>::shortestPath(int from, int to) {
	vector<pair<int, int>> path(numNodes, { INT_MAX, -1 }); // first is distance, second is previous
	if (containNegativeWeights) { // Bellman-Ford O(V*E)
		path[from].first = 0;
		auto edgeList = listEdges();
		for (int i = 0; i < numNodes - 1; i++)
			for (auto edge : edgeList)
				if (path[edge.first.second].first > path[edge.first.first].first + edge.second)
					path[edge.first.second] = { path[edge.first.first].first + edge.second, edge.first.first };
	}
	else { // Dijkstra's, O(E + VLOG(V))
		vector<bool> visited(numNodes);
		priority_queue<pair<int, int>> q;
		path[from].first = 0;
		q.push({ 0, from });
		while (!q.empty()) {
			int a = q.top().second; q.pop();
			visited[a] = 1;
			for (auto u : nodes[a].adjList) {
				if (path[a].first + u.second < path[u.first].first) {
					path[u.first].first = path[a].first + u.second;
					path[u.first].second = a;
					q.push({ -1 * path[u.first].first, u.first });
				}
			}
		}
	}
	vector<int> patho;
	if (path[to].second == -1) return { -1, patho };
	patho.push_back(to);
	int previous = path[to].second;
	while (previous != from) {
		patho.push_back(previous);
		previous = path[previous].second;
	}
	patho.push_back(from);
	reverse(patho.begin(), patho.end());
	return { path[to].first, patho };
}

template <class T>
vector<pair<pair<int, int>, int>> SEGL<T>::getEdgeList() {
	return listEdges();
}

template <class T>
bool SEGL<T>::topSort(vector<int>& sorted, SEGL<T>& new_graph, vector<bool>& finished, int& finishedCount) {
	int temp = -1;
	for (int i = 0; i < numNodes; i++)
		if (nodes_transpose[i].adjList.size() == 0 && !finished[i]) {
			temp = i;
			break;
		}
	if (temp == -1) if (finishedCount < numNodes) return false;
	if (finishedCount == numNodes) return true;
	sorted.push_back(temp);
	finished[temp] = 1;
	finishedCount++;
	auto lis = new_graph.nodes[temp].adjList;
	for (auto u : list)
		new_graph.removeDirEdge(temp, u.f);
	return topSort(sorted, new_graph, finished, finishedCount);
}

template <class T>
pair<bool, vector<int>> SEGL<T>::topSort() {
	vector<bool> finished(numNodes);
	vector<int> topsort;
	SEGL<T> new_graph = *this;
	int finishedCount = 0;
	bool b = topSort(topsort, new_graph, finished, finishedCount);
	return { b, topsort };
}

template <class T>
bool SEGL<T>::isDAG() {
	vector<bool> finished(numNodes);
	vector<int> topsort;
	SEGL<T> new_graph = *this;
	int finishedCount = 0;
	return topSort(topsort, new_graph, finished, finishedCount);
}

//template <class T>
//bool SEGL<T>::existsHamiltonianPath() { // Woodwall's theorem
//	bool flag = false;
//	for (int u = 0; u < numNodes - 1; u++) // checking every two nodes
//		for (int v = u + 1; j < numNodes; v++) {
//			flag = flag || nodes[u].adjList.count(v);
//			flag = flag || (this->outdegree(u) + this->indegree(v) > numNodes);
//		}
//	if (flag) return true;
//	// INCOMPLETE
//}

template <class T>
bool SEGL<T>::existsEulerianPath() {
	int count1 = 0, count2 = 0;
	bool flag = false;
	for (int i = 0; i < numNodes; i++) {
		if (indegree(i) != outdegree(i)) {
			if (indegree(i) - outdegree(i) == 1) count1++;
			else if (outdegree(i) - indegree(i) == 1) count2++;
			else {
				flag = true;
				break;
			}
		}
	}
	if (flag) return false;
	if (!(count1 == 0 && count2 == 0) && !(count1 == 1 && count2 == 1)) return false;
	SEGL<T> new_graph = *this;
	// make the graph undirected
	for (int u : new_graph.hasOutdegree)
		for (auto v : new_graph.nodes[u].adjList)
			new_graph.addDirWeightedEdge(v.f, u, v.s);
	vector<bool> visited(numNodes);
	int count = 0;
	for (int i = 0; i < numNodes; i++)
		if (!visited[i] && new_graph.outdegree(i) != 0) {
			dfs(i, visited);
			count++;
		}
	if (count > 1) return false;
	else return true;
}

#endif 