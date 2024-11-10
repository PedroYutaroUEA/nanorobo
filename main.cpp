#include <iostream>
#include <limits>
#include <list>

typedef int Vertex;
typedef float Weight;

const Weight INF = std::numeric_limits<Weight>::infinity();

class Node
{
private:
  Vertex vertex;
  Weight weight;
public:
  Node(Vertex vertex, Weight weight) : vertex(vertex), weight(weight) {}
  Vertex getVertex() const { return vertex; }
  Weight getWeight() const { return weight; }
  void setVertex(Vertex vertex) { this->vertex = vertex; }
  void setWeight(Weight weight) { this->weight = weight; }
};

struct Edge
{
  Vertex u, v;
  Weight weight;
};

template <typename T>
class UnionFind {
private:
  T *vertices, *parent;
public:
  UnionFind(int nv) : vertices(new T[nv + 1]), parent(new T[nv + 1]) {}
  ~UnionFind() { delete[] vertices; delete[] parent; }
  void makeSet(T vertex) { parent[vertex] = vertex; }
  T findSet(T vertex)
  {
    if (parent[vertex] != vertex)
      parent[vertex] = findSet(parent[vertex]);
    return parent[vertex];
  }
  void unionSet(T vertex1, T vertex2)
  {
    T parent1 = findSet(vertex1);
    T parent2 = findSet(vertex2);
    if (parent1 != parent2)
      parent[parent1] = parent2;
  }
};

template <typename T>
class Graph {
private:
  int num_vertices, num_edges;
  std::list<T> *adj;
  bool *healthStatus;
public:
  Graph(int nv) : num_vertices(nv), adj(new std::list<T>[nv + 1]), healthStatus(new bool[nv + 1])
  {
    for (int i = 1; i <= nv; i++)
      healthStatus[i] = true;
  }
  void destructor() { delete[] adj; delete[] healthStatus; }
  void addEdge(Vertex u, Vertex v, Weight w)
  {
    adj[u].emplace_back(v, w);
    adj[v].emplace_back(u, w);
    num_edges++;
  }
  std::list<T> &getAdj(Vertex vertex)
  {
    return adj[vertex];
  }
  int getNumVertices() const { return num_vertices; }
  int getNumEdges() const { return num_edges; }
  bool getHealthStatus(Vertex vertex) const { return healthStatus[vertex]; }
  void setHealthStatus(Vertex vertex, bool status)
  {
    healthStatus[vertex] = status;
  }
  void inputNodes(int numEdges)
  {
    for (int i = 1; i <= numEdges; i++)
    {
      Vertex u, v;
      Weight w;
      std::cin >> u >> v >> w;
      addEdge(u, v, w);
    }
  }
};

template <typename T>
class PriorityQueue {
private:
  T *heap;
  int *pos, size;
  Vertex parent(int i) { return i / 2; }

  Vertex leftChild(int i) { return 2 * i; }
  Vertex rightChild(int i) { return 2 * i + 1; }
  void heapifyUp(int i)
  {
    while (i > 1 && heap[i].second < heap[parent(i)].second)
    {
      swap(heap[i], heap[parent(i)]);
      pos[heap[i].first] = i;
      pos[heap[parent(i)].first] = parent(i);
      i = parent(i);
    }
  }
  void heapifyDown(int i)
  {
    while (i <= size / 2)
    {
      int smallest = i;

      if (leftChild(i) < size &&
          heap[leftChild(i)].second < heap[smallest].second)
        smallest = leftChild(i);
      if (rightChild(i) < size &&
          heap[rightChild(i)].second < heap[smallest].second)
        smallest = rightChild(i);
      if (smallest != i) {
        swap(heap[i], heap[smallest]);
        pos[heap[i].first] = i;
        pos[heap[smallest].first] = smallest;
        i = smallest;
      } else
        break;
    }
  }

public:
  PriorityQueue(int num_vertices)
      : heap(new T[num_vertices + 1]),
        pos(new int[num_vertices + 1]), size(0)
  {
    for (int i = 1; i <= num_vertices + 1; i++)
      pos[i] = 0;
  }
  ~PriorityQueue()
  {
    delete[] heap;
    delete[] pos;
  }
  bool isEmpty() const { return size == 0; }
  void insert(Vertex u, Weight w)
  {
    size += 1;
    heap[size] = {u, w};
    pos[u] = size;
    heapifyUp(size);
  }
  T extractMin()
  {
    T min = heap[1];
    heap[1] = heap[size];
    pos[heap[1].first] = 1;
    size -= 1;
    heapifyDown(1);
    return min;
  }
  void decreaseKey(Vertex v, Weight w)
  {
    int i = pos[v];
    if (i != 0 && heap[i].second > w)
    {
      heap[i].second = w;
      heapifyUp(i);
    }
  }
};

template <typename T, typename N>
class Brain {
private:
  int num_vertices, num_edges;
  std::list<N> *adj;
  T **graph;
  bool *healthStatus;
public:
  Brain(int nv) : num_vertices(nv), num_edges(0)
  {
    adj = new std::list<N>[nv + 1];
    graph = new T *[nv + 1]();
    healthStatus = new bool[nv + 1];

    for (int i = 1; i <= nv; i++)
      healthStatus[i] = false;
  }

  ~Brain()
  {
    for (int i = 1; i <= num_vertices; i++)
    {
      delete graph[i];
    }
    delete[] graph;
    delete[] adj;
    delete[] healthStatus;
  }

  void createGraph(int i, int size)
  {
    graph[i] = new T(size);
  }
  void addEdge(Vertex u, Vertex v, Weight w)
  {
    adj[u].emplace_back(v, w);
    adj[v].emplace_back(u, w);
    num_edges++;
  }
  int getNumBlocks() const { return num_vertices; }
  T &getGraph(Vertex v) const { return *graph[v]; }
  int getNumEdges() const { return num_edges; }
  int getNumVertices() const { return num_vertices; }
  bool getHealthStatus(Vertex v) const { return healthStatus[v]; }
  void setHealthStatus(Vertex v, bool status) { healthStatus[v] = status; }
  std::list<N> &getAdj(Vertex v) { return adj[v]; }
  void inputBrain(int numEdges)
  {
    for (int i = 1; i <= numEdges; i++)
    {
      Vertex u, v;
      Weight w;
      std::cin >> u >> v >> w;
      addEdge(u, v, w);
    }
  }
  void inputGraphs()
  {
    for (int i = 1; i <= getNumBlocks(); i++)
    {
      int num_vertices, num_edges, num_sickness;
      std::cin >> num_vertices >> num_edges;
      createGraph(i, num_vertices);

      std::cin >> num_sickness;
      if (num_sickness > 0)
      {
        setHealthStatus(i, true);
        Vertex sick_neuron;
        for (int j = 1; j <= num_sickness; j++)
        {
          std::cin >> sick_neuron;
          getGraph(i).setHealthStatus(sick_neuron, true);
        }
      }

      getGraph(i).inputNodes(num_edges);
    }
  }
};

void initialize(const int numBlocks, Weight *distances, Vertex *previous, Vertex root)
{
  for (int i = 1; i <= numBlocks; i++)
  {
    distances[i] = INF;
    previous[i] = -1;
  }
  distances[root] = 0;
}
template<typename T>
void relax(PriorityQueue<T> &pq, Vertex u, Vertex v, Weight w, Weight *dist, Vertex *prev)
{
  Weight oldDist = 0;
  if (dist[u] + w < dist[v])
  {
    oldDist = dist[v];
    dist[v] = dist[u] + w;
    prev[v] = u;
    if (oldDist != INF)
      pq.decreaseKey(v, dist[v]);
    else
      pq.insert(v, dist[v]);
  }
}
  
template <typename T, typename N>
void dijkstra(Brain<T, N> &b, Vertex root, Weight *distances, Vertex *previous)
{
  PriorityQueue<std::pair<Vertex, Weight>> pq(b.getNumBlocks());
  initialize(b.getNumBlocks(), distances, previous, root);
  pq.insert(root, 0);

  while (!pq.isEmpty())
  {
    Vertex u = pq.extractMin().first;
    for (const auto &node : b.getAdj(u))
    {
      Vertex v = node.getVertex();
      Weight weight = node.getWeight();
      relax<std::pair<Vertex, Weight>>(pq, u, v, weight, distances, previous);
    }
  }
}
template <typename T>
std::list<Edge> getOrderedEdges(Graph<T> &g)
{
  std::list<Edge> edges;
  for (Vertex u = 1; u <= g.getNumVertices(); ++u)
  {
    std::list<T> &adj_list = g.getAdj(u);
    for (auto &node : adj_list)
      edges.push_back({u, node.getVertex(), node.getWeight()});
  }
  edges.sort(
      [](const Edge &e1, const Edge &e2) { return e1.weight < e2.weight; });
  return edges;
}
template <typename T>
Weight kruskal(Graph<T> &g)
{
  int num_vertices = g.getNumVertices();
  UnionFind<Vertex> uf(num_vertices);
  Weight total_weight = 0.0;
  // make sets for each vertex
  for (Vertex v = 1; v <= num_vertices; ++v)
    uf.makeSet(v);
  // order edges
  std::list<Edge> edges = getOrderedEdges(g);
  for (auto &edge : edges)
  {
    Vertex u = edge.u;
    Vertex v = edge.v;
    Weight weight = edge.weight;
    if (uf.findSet(u)!= uf.findSet(v))
    {
      uf.unionSet(u, v);
      total_weight += weight;
    }
  }
  return total_weight;
}

template <typename T, typename N>
Brain<T, N> createBrain()
{
  int num_vertices, num_edges;
  std::cin >> num_vertices >> num_edges;
  Brain<T, N> brain(num_vertices);
  brain.inputBrain(num_edges);
  return brain;
}
template <typename T, typename N>
Weight executeProcedure(Brain<T, N> &b, Vertex exit, Vertex *previous)
{
  Weight sumOfWeights = 0;
  Vertex current_v = exit;
  Vertex end = -1;
  while (current_v != end)
  {
    if (b.getHealthStatus(current_v) == true)
      sumOfWeights += kruskal(b.getGraph(current_v));
    current_v = previous[current_v];
  }
  
  delete[] previous;
  return sumOfWeights;
}
template <typename T, typename N>
Weight createTrajetory(Brain<T, N> &b, Vertex root, Vertex exit)
{
  Weight *distances = new Weight[b.getNumBlocks() + 1];
  Vertex *previous = new Vertex[b.getNumBlocks() + 1];
  dijkstra(b, root, distances, previous);

  delete[] distances;
  return executeProcedure(b, exit, previous);
}
int main()
{
  Brain<Graph<Node>, Node> brain = createBrain<Graph<Node>, Node>();
  
  Vertex root, exit;
  std::cin >> root >> exit;

  brain.inputGraphs();
  
  Weight result = createTrajetory(brain, root, exit);

  std::cout << result;

  for (int i = 1; i <= brain.getNumBlocks(); i++)
    brain.getGraph(i).destructor();
  return 0;
}