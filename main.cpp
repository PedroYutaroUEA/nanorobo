#include <iostream>
#include <list>
#include <limits>
#include <vector>
#define INF numeric_limits<Weight>::max()
using namespace std;

typedef unsigned int Vertex;
typedef float Weight, Dist;

struct Edge
{
  Vertex u, v;
  Weight weight;
};

class Graph
{
private:
  int num_vertices;
  list<pair<Vertex, Weight>> *adj;
  vector<Vertex> sick_neuroms;

public:
  Graph() : num_vertices(0), adj(nullptr) {}

  Graph(int num_vertices) : num_vertices(num_vertices)
  {
    adj = new list<pair<Vertex, Weight>>[num_vertices + 1];
  }

  ~Graph() { delete[] adj; }

  void addEdge(Vertex u, Vertex v, Weight w)
  {
    adj[u].push_back(make_pair(v, w));
    adj[v].push_back(make_pair(u, w));
  }

  const list<pair<Vertex, Weight>> &getAdj(Vertex vertex) const
  {
    return adj[vertex];
  }

  void setSickNeuroms(const vector<Vertex> &neuroms)
  {
    for (auto n : neuroms)
      sick_neuroms.push_back(n);
  }

  int getNumVertices() const
  {
    return num_vertices;
  }
};

template <typename T>
class PriorityQueue
{
private:
  vector<T> heap;

  void heapifyUp(int index)
  {
    while (index > 0 && heap[(index - 1) / 2].first > heap[index].first)
    {
      swap(heap[index], heap[(index - 1) / 2]);
      index = (index - 1) / 2;
    }
  }

  void heapifyDown(int index)
  {
    int size = heap.size();
    int smallest = index;
    int left = 2 * index + 1;
    int right = 2 * index + 2;

    if (left < size && heap[left].first < heap[smallest].first)
      smallest = left;
    if (right < size && heap[right].first < heap[smallest].first)
      smallest = right;
    if (smallest != index)
    {
      swap(heap[index], heap[smallest]);
      heapifyDown(smallest);
    }
  }

public:
  bool isEmpty() const
  {
    return heap.empty();
  }

  void insert(T node)
  {
    heap.push_back(node);
    heapifyUp(heap.size() - 1);
  }

  T extractMin()
  {
    if (heap.empty())
      throw runtime_error("Heap is empty");
    T minNode = heap[0];
    heap[0] = heap.back();
    heap.pop_back();
    if (!heap.empty())
      heapifyDown(0);
    return minNode;
  }

  T extractMax()
  {
    if (heap.empty())
      throw runtime_error("Heap is empty");

    // Encontrar o nó com o valor máximo
    int maxIndex = 0;
    for (int i = 1; i < heap.size(); ++i)
    {
      if (heap[i].first > heap[maxIndex].first)
      {
        maxIndex = i;
      }
    }
    T maxNode = heap[maxIndex];
    heap[maxIndex] = heap.back();
    heap.pop_back();
    if (!heap.empty())
      heapifyDown(maxIndex);
    return maxNode;
  }

  T getMin() const
  {
    if (heap.empty())
      throw runtime_error("Heap is empty");
    return heap[0];
  }

  T getMax() const
  {
    if (heap.empty())
      throw runtime_error("Heap is empty");

    // Encontrar o nó com o valor máximo
    int maxIndex = 0;
    for (int i = 1; i < heap.size(); ++i)
    {
      if (heap[i].first > heap[maxIndex].first)
      {
        maxIndex = i;
      }
    }
    return heap[maxIndex];
  }

  void decreaseKey(Vertex v, Dist newDist)
  {
    for (size_t i = 0; i < heap.size(); i++)
    {
      if (heap[i].second == v)
      {
        heap[i].first = newDist;
        heapifyUp(i);
        return;
      }
    }
  }
};

class UnionFind
{
private:
  vector<Vertex> vertices, parent;

public:
  UnionFind(int num_vertices) : vertices(num_vertices), parent(num_vertices) {}

  void makeSet(Vertex vertex)
  {
    parent[vertex] = vertex;
  }

  Vertex findSet(Vertex vertex)
  {
    if (parent[vertex] != vertex)
      parent[vertex] = findSet(parent[vertex]);
    return parent[vertex];
  }

  void unionSet(Vertex vertex1, Vertex vertex2)
  {
    Vertex parent1 = findSet(vertex1);
    Vertex parent2 = findSet(vertex2);
    if (parent1 != parent2)
      parent[parent1] = parent2;
  }
};

template <typename T>
Weight MST_Kurskal(Graph &g)
{
  int num_vertices = g.getNumVertices();
  UnionFind uf(num_vertices);
  Weight total_weight = 0.0;
  // make sets for each vertex
  for (Vertex v = 0; v < num_vertices; ++v)
    uf.makeSet(v);
  // order edges
  list<Edge> edges;
  for (Vertex u = 0; u < num_vertices; ++u)
  {
    const list<pair<Vertex, Weight>> &adj_list = g.getAdj(u);
    for (auto &edge : adj_list)
      edges.push_back({u, edge.first, edge.second});
  }
  edges.sort(
      [](const Edge &e1, const Edge &e2)
      { return e1.weight < e2.weight; });
  for (auto &edge : edges)
  {
    Vertex u = edge.u;
    Vertex v = edge.v;
    Weight weight = edge.weight;
    if (uf.findSet(u) != uf.findSet(v))
    {
      uf.unionSet(u, v);
      total_weight += weight;
    }
  }
  return total_weight;
}

class Brain
{
private:
  int num_vertices;
  vector<Graph> blocks;
  vector<vector<pair<Vertex, Weight>>> conections;

public:
  Brain(int num_vertices) : num_vertices(num_vertices)
  {
    blocks.resize(num_vertices + 1);
    conections.resize(num_vertices + 1);
  }
  void addConection(Vertex u, Vertex v, Weight w)
  {
    conections[u].emplace_back(v, w);
    conections[v].emplace_back(u, w);
  }
  void addBlock(Vertex id, Graph &g)
  {
    blocks[id] = g;
  }

  vector<pair<Vertex, Weight>> getAdj(Vertex vertex) const
  {
    vector<pair<Vertex, Weight>> adj;
    for (const auto& connection : conections[vertex]) {
      adj.push_back(connection);
    }
    return adj;
  }

  int getNumBlocks() const
  {
    return num_vertices;
  }
};
void dijkstra(const Brain &brain, Vertex root, Vertex exit, vector<Dist> &distances, vector<Vertex> &previous)
{
  int num_vertices = brain.getNumBlocks();
  distances.assign(num_vertices, INF);
  previous.assign(num_vertices, INF);

  PriorityQueue<pair<Dist, Vertex>> pq;
  distances[root] = 0;
  pq.insert(make_pair(0, root));

  while (!pq.isEmpty())
  {
    Vertex u = pq.extractMin().second;
    if (u == exit)
      return;
    for (const auto& connection : brain.getAdj(u))
    {
      Vertex v = connection.first;
      Weight weight = connection.second;
      if (distances[v] > distances[u] + weight)
      {
        distances[v] = distances[u] + weight;
        previous[v] = u;
        pq.insert(make_pair(distances[v], v));
      }
    }
  }
}
int main()
{
  int n_conections, n_blocks;
  cin >> n_blocks >> n_conections;
  Brain brain(n_blocks);
  Vertex u, v;
  Weight w;
  for (int i = 1; i <= n_conections; i++) {
    cin >> u >> v >> w;
    brain.addConection(u, v, w);
  }
  int root_block, exit_block;
  cin >> root_block >> exit_block;

  int n_vertices, n_edges, num_sickness;
  for (int i = 1; i <= brain.getNumBlocks(); i++) {
    cin >> n_vertices >> n_edges;
    Graph g(n_vertices);
    cin >> num_sickness;
    if (num_sickness > 0)
    {
      vector<Vertex> sick_neuroms(num_sickness);
      for (int j = 0; j < num_sickness; j++)
        cin >> sick_neuroms[j];
      g.setSickNeuroms(sick_neuroms);
    }
    for (int k = 1; k <= n_edges; k++) {
      cin >> u >> v >> w;
      g.addEdge(u, v, w);
    }
    brain.addBlock(i, g);
  }
  vector<Dist> distances;
  vector<Vertex> previous;
  dijkstra(brain, root_block, exit_block, distances, previous);

  for (Vertex i = 0; i < previous.size(); i++)
    cout << previous[i] << endl;
}