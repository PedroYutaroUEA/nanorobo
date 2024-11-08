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

class PriorityQueue
{
private:
  pair<Vertex, Weight> *heap;
  int *pos, size;
  Vertex parent(int i)
  {
    return i / 2;
  }

  Vertex leftChild(int i)
  {
    return 2 * i;
  }
  Vertex rightChild(int i)
  {
    return 2 * i + 1;
  }
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

      if (leftChild(i) < size && heap[leftChild(i)].second < heap[smallest].second)
        smallest = leftChild(i);
      if (rightChild(i) < size && heap[rightChild(i)].second < heap[smallest].second)
        smallest = rightChild(i);
      if (smallest != i)
      {
        swap(heap[i], heap[smallest]);
        pos[heap[i].first] = i;
        pos[heap[smallest].first] = smallest;
        i = smallest;
      }
      else
        break;
    }
  }

public:
  PriorityQueue(int num_vertices) : heap(new pair<Vertex, Weight>[num_vertices + 1]), pos(new int[num_vertices + 1]), size(0)
  {
    for (int i = 1; i <= num_vertices + 1; i++)
      pos[i] = 0;
  }
  // ~PriorityQueue()
  // {
  //   delete[] heap;
  //   delete[] pos;
  // }
  bool isEmpty() const
  {
    return size == 0;
  }
  void insert(Vertex u, Weight w)
  {
    size += 1;
    heap[size] = {u, w};
    pos[u] = size;
    heapifyUp(size);
  }
  pair<Vertex, Weight> extractMin()
  {
    pair<Vertex, Weight> min = heap[1];
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
    for (const auto &connection : conections[vertex])
    {
      adj.push_back(connection);
    }
    return adj;
  }

  int getNumBlocks() const
  {
    return num_vertices;
  }

  void dijkstra(Vertex root, Vertex exit)
  {
    PriorityQueue pq(num_vertices);
    Weight *distances = new Weight[num_vertices + 1];
    Vertex *previous = new Vertex[num_vertices + 1];
    for (int i = 1; i <= num_vertices; i++)
    {
      distances[i] = INF;
      previous[i] = -1;
    }
    distances[root] = 0;
    pq.insert(root, 0);

    while (!pq.isEmpty())
    {
      Vertex u = pq.extractMin().first;

      for (const auto &connection : getAdj(u))
      {
        Vertex v = connection.first;
        Weight weight = connection.second;
        Weight oldDist = 0;

        if (distances[u] + weight < distances[v])
        {
          oldDist = distances[v];
          distances[v] = distances[u] + weight;
          previous[v] = u;
          if (oldDist != INF)
            pq.decreaseKey(v, distances[v]);
          else
            pq.insert(v, distances[v]);
        }
      }
    }

    cout << "Previous vertices in the shortest path:" << endl;
    for (int i = 1; i <= num_vertices + 1; i++)
    {
      cout << "Previous vertex for " << i << ": " << previous[i] << endl;
    }
  }
};
// void dijkstra(const Brain &brain, Vertex root, Vertex exit, vector<Dist> &distances, vector<Vertex> &previous)
// {
//   int num_vertices = brain.getNumBlocks();
//   distances.assign(num_vertices + 1, INF);
//   previous.assign(num_vertices + 1, -1);

//   PriorityQueue<pair<Dist, Vertex>> pq;
//   distances[root] = 0;
//   pq.insert(make_pair(0, root));

//   while (!pq.isEmpty())
//   {
//     Vertex u = pq.extractMin().second;

//     for (const auto &connection : brain.getAdj(u))
//     {
//       Vertex v = connection.first;
//       Weight weight = connection.second;
//       Weight oldDist = 0;
//       if (distances[v] > distances[u] + weight)
//       {
//         oldDist = distances[v];
//         distances[v] = distances[u] + weight;
//         previous[v] = u;
//         if (oldDist != INF)
//           pq.decreaseKey(v, distances[v]);
//         else
//           pq.insert(make_pair(distances[v], v));
//       }
//     }
//   }
// }
int main()
{
  int n_conections, n_blocks;
  cin >> n_blocks >> n_conections;
  Brain brain(n_blocks);
  Vertex u, v;
  Weight w;
  for (int i = 1; i <= n_conections; i++)
  {
    cin >> u >> v >> w;
    brain.addConection(u, v, w);
  }
  int root_block, exit_block;
  cin >> root_block >> exit_block;

  int n_vertices, n_edges, num_sickness;
  for (int i = 1; i <= brain.getNumBlocks(); i++)
  {
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
    for (int k = 1; k <= n_edges; k++)
    {
      cin >> u >> v >> w;
      g.addEdge(u, v, w);
    }
    brain.addBlock(i, g);
  }
  // vector<Dist> distances;
  // vector<Vertex> previous;
  // dijkstra(brain, root_block, exit_block, distances, previous);
  brain.dijkstra(root_block, exit_block);
}