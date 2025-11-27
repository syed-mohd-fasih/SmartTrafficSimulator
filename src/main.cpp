#include "raylib.h"

#include <cmath>
#include <cstdio>
#include <ctime>

// -----------------------------
// Utility Vec2f
// -----------------------------
/**
 * @brief Simple 2D vector structure for positions.
 */
struct Vec2f
{
    float x, y;
    Vec2f() : x(0), y(0) {}
    Vec2f(float _x, float _y) : x(_x), y(_y) {}
};

/**
 * @brief Calculates the Euclidean distance between two Vec2f points.
 */
static float distf(const Vec2f &a, const Vec2f &b)
{
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    return sqrtf(dx * dx + dy * dy);
}

/**
 * @brief Linear interpolation (LERP) between two Vec2f points.
 */
static Vec2f lerp(const Vec2f &a, const Vec2f &b, float t)
{
    return Vec2f(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
}

// -----------------------------
// DynArr: minimal dynamic array (with deep copy)
// -----------------------------
/**
 * @brief Minimal custom dynamic array (vector) with deep copy semantics.
 */
template <typename T>
struct DynArr
{
    T *data;      // Pointer to the array data
    int count;    // Current number of elements
    int capacity; // Maximum number of elements before realloc

    DynArr() : data(nullptr), count(0), capacity(0) {}

    // Deep copy constructor (creates a distinct copy of all elements)
    DynArr(const DynArr &other) : data(nullptr), count(0), capacity(0)
    {
        if (other.count > 0)
        {
            reserve(other.count);
            count = other.count;
            for (int i = 0; i < count; ++i)
                data[i] = other.data[i];
        }
    }

    // Deep copy assignment operator
    DynArr &operator=(const DynArr &other)
    {
        if (this != &other)
        {
            delete[] data;
            data = nullptr;
            count = 0;
            capacity = 0;

            if (other.count > 0)
            {
                reserve(other.count);
                count = other.count;
                for (int i = 0; i < count; ++i)
                    data[i] = other.data[i];
            }
        }
        return *this;
    }

    ~DynArr() { delete[] data; }

    void clear() { count = 0; }
    bool empty() const { return count == 0; }
    int size() const { return count; }

    /**
     * @brief Ensures capacity is at least 'newcap', reallocating if necessary.
     */
    void reserve(int newcap)
    {
        if (newcap <= capacity)
            return;
        int c = capacity == 0 ? 4 : capacity;
        while (c < newcap)
            c <<= 1;
        T *nd = new T[c];
        for (int i = 0; i < count; i++)
            nd[i] = data[i];
        delete[] data;
        data = nd;
        capacity = c;
    }

    /**
     * @brief Adds an element to the end of the array.
     */
    void push(const T &v)
    {
        if (count + 1 > capacity)
            reserve(count + 1);
        data[count++] = v;
    }

    void popBack()
    {
        if (count > 0)
            --count;
    }

    T &operator[](int i) { return data[i]; }
    const T &operator[](int i) const { return data[i]; }
};

// -----------------------------
// SimpleQueue: circular queue using raw array
// -----------------------------
/**
 * @brief Simple circular queue for integer IDs (used for vehicle queues).
 */
struct SimpleQueueInt
{
    int *buf;       // Array buffer
    int cap;        // Capacity
    int head, tail; // Indices for dequeue (head) and enqueue (tail)

    SimpleQueueInt() : buf(nullptr), cap(0), head(0), tail(0) {}
    ~SimpleQueueInt() { delete[] buf; }

    bool empty() const { return head == tail; }

    /**
     * @brief Reallocates the buffer to a larger capacity.
     */
    void reserve(int n)
    {
        if (n <= cap)
            return;

        int newcap = cap == 0 ? 8 : cap;
        while (newcap < n)
            newcap <<= 1;

        int *nb = new int[newcap];

        // Copy existing elements in order (from head to tail)
        int sz = 0;
        if (buf != nullptr && cap > 0 && !empty())
        {
            int idx = head;
            while (idx != tail)
            {
                nb[sz++] = buf[idx];
                idx = (idx + 1) % cap;
            }
        }

        delete[] buf;
        buf = nb;
        cap = newcap;
        head = 0;  // Reset head to start of new buffer
        tail = sz; // Set tail after last copied element
    }

    /**
     * @brief Adds an element to the back of the queue. Expands if full.
     */
    void push(int v)
    {
        if (cap == 0)
            reserve(8);

        int next = (tail + 1) % cap;
        if (next == head) // Check if full
        {
            reserve(cap * 2);
            next = (tail + 1) % cap;
        }
        buf[tail] = v;
        tail = next;
    }

    /**
     * @brief Removes and returns the element from the front of the queue.
     */
    int pop()
    {
        int v = buf[head];
        head = (head + 1) % cap;
        return v;
    }

    // Returns the element at the front without removing it.
    int front() const { return buf[head]; }
};

// -----------------------------
// MinHeap for (id,priority)
// -----------------------------
/**
 * @brief Structure for an element in the MinHeap.
 */
struct HeapNode
{
    int id;   // Node ID
    float pr; // Priority (distance/cost)
};

/**
 * @brief Minimal Binary MinHeap implementation for Dijkstra/A* priority queue.
 */
struct MinHeap
{
    HeapNode *buf;
    int cap;
    int len; // Current number of elements

    MinHeap() : buf(nullptr), cap(0), len(0) {}
    ~MinHeap() { delete[] buf; }

    void reserve(int n)
    {
        // ... (resize logic - ensures capacity is a power of 2)
        if (n <= cap)
            return;
        int c = cap == 0 ? 8 : cap;
        while (c < n)
            c <<= 1;
        HeapNode *nb = new HeapNode[c];
        for (int i = 0; i < len; i++)
            nb[i] = buf[i];
        delete[] buf;
        buf = nb;
        cap = c;
    }

    /**
     * @brief Inserts a node and maintains the min-heap property (bubble up).
     */
    void push(int id, float pr)
    {
        if (len + 1 > cap)
            reserve(len + 1);
        buf[len].id = id;
        buf[len].pr = pr;
        int i = len++;
        while (i > 0)
        {
            int p = (i - 1) >> 1;
            if (buf[p].pr <= buf[i].pr)
                break;
            HeapNode tmp = buf[p];
            buf[p] = buf[i];
            buf[i] = tmp;
            i = p;
        }
    }

    bool empty() const { return len == 0; }

    /**
     * @brief Removes and returns the element with the minimum priority (min-heap root).
     */
    HeapNode pop()
    {
        HeapNode out = buf[0];
        buf[0] = buf[len - 1];
        len--;
        int i = 0;
        // Heapify-down operation
        for (;;)
        {
            int l = 2 * i + 1, r = 2 * i + 2, smallest = i;
            if (l < len && buf[l].pr < buf[smallest].pr)
                smallest = l;
            if (r < len && buf[r].pr < buf[smallest].pr)
                smallest = r;
            if (smallest == i)
                break;
            HeapNode tmp = buf[i];
            buf[i] = buf[smallest];
            buf[smallest] = tmp;
            i = smallest;
        }
        return out;
    }

    void clear() { len = 0; }
};

// -----------------------------
// Graph structures
// -----------------------------
/**
 * @brief Represents a directed edge in the graph.
 */
struct Edge
{
    int to;  // Destination node ID
    float w; // Static weight/distance
};

/**
 * @brief Represents a node (intersection) in the graph.
 */
struct Node
{
    Vec2f pos;            // Position in 2D space
    int id;               // Unique node ID
    float lightTimer;     // Timer for the traffic light cycle
    int lightState;       // 0 = Green, 1 = Red
    SimpleQueueInt queue; // Queue of vehicles waiting at this intersection
    Node() : pos(), id(-1), lightTimer(0), lightState(0) {}
};

/**
 * @brief The main graph structure (nodes and adjacency list).
 */
struct Graph
{
    DynArr<Node> nodes;       // List of all nodes/intersections
    DynArr<DynArr<Edge>> adj; // Adjacency list: adj[i] is a list of edges originating from node i

    Graph() {}

    // Initializes graph structure capacity
    void init(int n)
    {
        nodes.clear();
        adj.clear();
        for (int i = 0; i < n; i++)
        {
            Node nd;
            nd.id = nodes.size();
            nodes.push(nd);
            DynArr<Edge> earr;
            earr.clear();
            adj.push(earr);
        }
    }

    // Adds a new node at position 'p'
    int addNode(const Vec2f &p)
    {
        Node nd;
        nd.id = nodes.size();
        nd.pos = p;
        nodes.push(nd);
        DynArr<Edge> darr;
        darr.clear();
        adj.push(darr);
        return nd.id;
    }

    // Adds a directed edge from node 'a' to node 'b' with weight 'w'
    void addEdge(int a, int b, float w)
    {
        Edge e;
        e.to = b;
        e.w = w;
        adj[a].push(e);
    }

    int size() const { return nodes.size(); }
};

// -----------------------------
// Vehicles
// -----------------------------
enum VehicleType
{
    NORMAL = 0,
    EMERGENCY = 1 // Emergency vehicles ignore queues/lights
};

/**
 * @brief Represents a user-controlled vehicle with a defined path.
 */
struct Vehicle
{
    int id;
    VehicleType type;
    int src, dest;      // Start and End nodes
    DynArr<int> path;   // Sequence of node IDs to follow
    int currentSegment; // Index of the current starting node in 'path'
    float posAlong;     // Normalized position (0.0 to 1.0) along the current segment
    float speed;
    bool finished;

    // Queue-related state
    bool waiting;    // True if currently stopped at a red light or behind another car
    int waitingNode; // The ID of the node (intersection) where the car is waiting

    Vehicle()
        : id(-1), type(NORMAL), src(-1), dest(-1),
          currentSegment(0), posAlong(0), speed(80.0f),
          finished(false), waiting(false), waitingNode(-1) {}
};

/**
 * @brief Represents a simple non-player controlled vehicle (NPC).
 */
struct NPCVehicle
{
    int id;
    int currentNode; // The node the NPC just left (start of segment)
    int nextNode;    // The node the NPC is traveling towards (end of segment)
    float posAlong;  // Normalized position (0.0 to 1.0) along the segment
    float speed;
    bool waiting; // True if currently stopped at an intersection
    int waitingNode;
    NPCVehicle() : id(-1), currentNode(0), nextNode(0), posAlong(0.0f), speed(40.0f), waiting(false), waitingNode(-1) {}
};

// -----------------------------
// Globals
// -----------------------------
static Graph G;
static DynArr<Vehicle> vehicles;
static MinHeap emergencyHeap;
static int nextVehicleId = 1;
static int selectedStart = -1, selectedEnd = -1;
static bool showIds = true;
static bool runningSim = false;
static bool pathFound = false;
static DynArr<int> lastPath;    // Stores the last calculated path
static bool useAstar = false;   // True if A* was used, false for Dijkstra
static float lightCycle = 6.0f; // Total duration of a traffic light cycle (e.g., 3s green, 3s red)
static int screenW = 1280, screenH = 720;

// NPC globals
static DynArr<NPCVehicle> npcVehicles;
static int nextNPCId = 100000;

// -----------------------------
// Build simple grid graph
// -----------------------------
/**
 * @brief Initializes the graph as a 2D grid of interconnected nodes.
 */
void buildGridGraph(int cols, int rows, float marginX, float marginY, float gx, float gy)
{
    G.init(0);
    // Create nodes based on grid dimensions
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            Vec2f p(marginX + c * gx, marginY + r * gy);
            Node nd;
            nd.id = G.nodes.size();
            nd.pos = p;
            G.nodes.push(nd);
            DynArr<Edge> arr;
            arr.clear();
            G.adj.push(arr);
        }
    }

    // Lambda to get node index from grid coordinates
    auto idx = [&](int c, int r)
    { return r * cols + c; };

    // Create edges (roads) between adjacent nodes
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            int a = idx(c, r);
            if (c + 1 < cols)
            {
                // Horizontal edge (and back)
                int b = idx(c + 1, r);
                float w = distf(G.nodes[a].pos, G.nodes[b].pos);
                G.addEdge(a, b, w);
                G.addEdge(b, a, w);
            }
            if (r + 1 < rows)
            {
                // Vertical edge (and back)
                int b = idx(c, r + 1);
                float w = distf(G.nodes[a].pos, G.nodes[b].pos);
                G.addEdge(a, b, w);
                G.addEdge(b, a, w);
            }
        }
    }

    // Initialize traffic light state and queues for all nodes
    for (int i = 0; i < G.size(); ++i)
    {
        G.nodes[i].lightTimer = (float)(i % 3); // Stagger initial timers
        G.nodes[i].lightState = 0;              // Start state (Green)
        G.nodes[i].queue.reserve(8);
    }
}

// -----------------------------
// Pathfinding (Dijkstra + A*)
// -----------------------------
static DynArr<float> g_scores; // Distance from start node (g(u))
static DynArr<int> cameFrom;   // Parent node in the shortest path tree

/**
 * @brief Ensures g_scores and cameFrom arrays are the correct size and reset.
 */
void ensureArrSize(int n)
{
    // ... (logic to resize and reset arrays to INFINITY / -1)
    if (g_scores.size() < n)
    {
        int cur = g_scores.size();
        for (int i = cur; i < n; i++)
            g_scores.push(INFINITY);
    }
    else
    {
        for (int i = 0; i < n; i++)
            g_scores[i] = INFINITY;
    }
    if (cameFrom.size() < n)
    {
        int cur = cameFrom.size();
        for (int i = cur; i < n; i++)
            cameFrom.push(-1);
    }
    else
    {
        for (int i = 0; i < n; i++)
            cameFrom[i] = -1;
    }
}

/**
 * @brief Traces back the path from target 't' to source 's' using the cameFrom array.
 */
bool reconstructPath(int s, int t, DynArr<int> &out)
{
    out.clear();
    if (cameFrom[t] == -1 && s != t)
        return false; // No path found
    int cur = t;
    while (cur != -1)
    {
        out.push(cur);
        if (cur == s)
            break;
        cur = cameFrom[cur];
    }
    // Reverse the path (currently t -> s) to be (s -> t)
    for (int i = 0, j = out.size() - 1; i < j; ++i, --j)
    {
        int tmp = out[i];
        out[i] = out[j];
        out[j] = tmp;
    }
    return true;
}

/**
 * @brief Calculates the cost of traversing an edge, factoring in traffic queue delays.
 * This function provides the dynamic weight used by A* and Dijkstra for user vehicles.
 */
float getDynamicEdgeWeight(int fromNode, int toNode, float staticDistance)
{
    // The cost is calculated based on delays at the 'toNode' intersection.
    Node &targetNode = G.nodes[toNode];

    // Base estimated travel time without delays (Time = Distance / Speed)
    float baseTime = staticDistance / 80.0f;
    float penaltyFactor = 1.0f; // Penalty: 1 second per car in the queue

    float trafficDelay = 0.0f;

    // Apply delay ONLY if the light is RED (1)
    if (targetNode.lightState == 1)
    {
        // Calculate current queue size at the target intersection
        int queueSize = (targetNode.queue.tail - targetNode.queue.head + targetNode.queue.cap) % targetNode.queue.cap;

        // Apply a penalty based on queue length
        trafficDelay = queueSize * penaltyFactor;
    }

    // The final weight is the estimated time = base travel time + traffic delay.
    return baseTime + trafficDelay;
}

/**
 * @brief Runs Dijkstra's algorithm from source 's' to target 't' using dynamic weights.
 */
bool runDijkstra(int s, int t)
{
    int n = G.size();
    ensureArrSize(n);
    MinHeap heap;
    heap.reserve(n + 4);

    g_scores[s] = 0.0f;
    heap.push(s, 0.0f); // Push start node

    while (!heap.empty())
    {
        HeapNode hn = heap.pop();
        int u = hn.id;
        float d = hn.pr;

        if (d != g_scores[u])
            continue; // Skip stale entries

        if (u == t)
        {
            // Path found, reconstruct and exit
            bool ok = reconstructPath(s, t, lastPath);
            pathFound = ok;
            return ok;
        }

        // Iterate over neighbors
        for (int ei = 0; ei < G.adj[u].size(); ++ei)
        {
            Edge &e = G.adj[u][ei];
            int v = e.to;

            // Calculate dynamic, traffic-aware edge cost (weight)
            float currentDistance = distf(G.nodes[u].pos, G.nodes[v].pos);
            float edgeCost = getDynamicEdgeWeight(u, v, currentDistance);

            float nd = d + edgeCost; // New accumulated cost (time)

            if (nd < g_scores[v])
            {
                g_scores[v] = nd;
                cameFrom[v] = u;
                heap.push(v, nd); // Push/update neighbor in the heap
            }
        }
    }
    pathFound = false;
    return false;
}

/**
 * @brief Runs A* algorithm from source 's' to target 't' using dynamic weights.
 */
bool runAstar(int s, int t)
{
    int n = G.size();
    ensureArrSize(n);
    MinHeap heap;
    heap.reserve(n + 4);
    DynArr<float> fscore; // f(u) = g(u) + h(u)
    fscore.clear();
    for (int i = 0; i < n; i++)
        fscore.push(INFINITY);

    g_scores[s] = 0.0f;
    fscore[s] = distf(G.nodes[s].pos, G.nodes[t].pos); // Initial F-score (heuristic only)
    heap.push(s, fscore[s]);

    while (!heap.empty())
    {
        HeapNode hn = heap.pop();
        int u = hn.id;

        if (u == t)
        {
            // Path found
            bool ok = reconstructPath(s, t, lastPath);
            pathFound = ok;
            return ok;
        }

        // Iterate over neighbors
        for (int ei = 0; ei < G.adj[u].size(); ++ei)
        {
            Edge &e = G.adj[u][ei];
            int v = e.to;

            // Calculate dynamic, traffic-aware edge cost (weight)
            float currentDistance = distf(G.nodes[u].pos, G.nodes[v].pos);
            float edgeCost = getDynamicEdgeWeight(u, v, currentDistance);

            float tentative = g_scores[u] + edgeCost; // New G-score (cost from start)

            if (tentative < g_scores[v])
            {
                cameFrom[v] = u;
                g_scores[v] = tentative;

                // F-score: new G-score + static distance heuristic (ensures admissibility)
                float est = tentative + distf(G.nodes[v].pos, G.nodes[t].pos);
                heap.push(v, est);
            }
        }
    }
    pathFound = false;
    return false;
}

// -----------------------------
// Vehicles management
// -----------------------------
/**
 * @brief Creates a new vehicle along the currently calculated path.
 */
Vehicle *spawnVehicle(int src, int dest, VehicleType type)
{
    // ... (checks and vehicle setup)
    if (src < 0 || dest < 0 || !pathFound)
        return nullptr;

    Vehicle v;
    v.id = nextVehicleId++;
    v.type = type;
    v.src = src;
    v.dest = dest;
    v.currentSegment = 0;
    v.posAlong = 0.0f;
    v.finished = false;
    v.waiting = false;
    v.waitingNode = -1;
    // Copy the calculated path
    v.path.clear();
    for (int i = 0; i < lastPath.size(); ++i)
        v.path.push(lastPath[i]);
    v.speed = (type == EMERGENCY) ? 140.0f : 80.0f;

    vehicles.push(v);
    Vehicle *pv = &vehicles[vehicles.size() - 1];

    if (type == EMERGENCY)
        emergencyHeap.push(pv->id, 0.0f);

    return pv;
}

/**
 * @brief Resets all vehicles and simulation state.
 */
void resetVehicles()
{
    vehicles.clear();
    npcVehicles.clear();
    emergencyHeap.clear();
    lastPath.clear();
    nextVehicleId = 1;
    selectedStart = -1, selectedEnd = -1;
    runningSim = false;
    pathFound = false;
    useAstar = false;
    nextNPCId = 100000;

    // Clear all traffic queues
    for (int i = 0; i < G.size(); i++)
    {
        Node &n = G.nodes[i];
        while (!n.queue.empty())
            n.queue.pop();
    }
}

// -----------------------------
// NPC control
// -----------------------------
/**
 * @brief Spawns a number of non-player controlled vehicles at random locations.
 */
void spawnNPCs(int count = 1)
{
    if (G.size() == 0)
        return;
    for (int k = 0; k < count; ++k)
    {
        NPCVehicle n;
        n.id = nextNPCId++;
        // Choose random start node
        n.currentNode = GetRandomValue(0, G.size() - 1);
        // Choose random next node from neighbors
        if (G.adj[n.currentNode].size() > 0)
        {
            int idx = GetRandomValue(0, G.adj[n.currentNode].size() - 1);
            n.nextNode = G.adj[n.currentNode][idx].to;
        }
        else
            n.nextNode = n.currentNode;
        n.posAlong = 0.0f;
        n.speed = (float)GetRandomValue(30, 60); // Slower than user cars
        n.waiting = false;
        n.waitingNode = -1;
        npcVehicles.push(n);
    }
}

/**
 * @brief Removes NPCs from the simulation.
 */
void removeNPCs(int count = 1)
{
    for (int i = 0; i < count; ++i)
    {
        if (!npcVehicles.empty())
            npcVehicles.popBack();
    }
}

// -----------------------------
// Simulation updates
// -----------------------------
/**
 * @brief Updates the timer and state for all traffic lights.
 */
void updateTrafficLights(float dt)
{
    for (int i = 0; i < G.size(); ++i)
    {
        Node &n = G.nodes[i];

        // Update light timer
        n.lightTimer += dt;
        if (n.lightTimer >= lightCycle)
            n.lightTimer -= lightCycle;

        // State: 0 = Green (first half), 1 = Red (second half)
        n.lightState = (n.lightTimer < (lightCycle * 0.5f)) ? 0 : 1;
    }
}

/**
 * @brief Updates the position and traffic logic for user-controlled vehicles.
 */
void updateVehicles(float dt)
{
    for (int i = 0; i < vehicles.size(); ++i)
    {
        Vehicle &v = vehicles[i];
        if (v.finished || v.path.size() < 2 || v.currentSegment >= v.path.size() - 1)
        {
            v.finished = true;
            continue;
        }

        int curNode = v.path[v.currentSegment];
        int nextNode = v.path[v.currentSegment + 1];
        Node &n = G.nodes[curNode];
        bool green = (n.lightState == 0);

        // --- Intersection Arrival/Waiting Logic ---
        if (v.posAlong <= 0.001f) // Just arrived at the intersection (curNode)
        {
            if (v.type == EMERGENCY)
            {
                // Emergency vehicles ignore lights and queues
                v.waiting = false;
                v.waitingNode = -1;
            }
            else // NORMAL vehicle logic
            {
                // Condition to proceed: Green Light AND (Queue is empty OR Vehicle is at the front)
                bool canProceed = green && (n.queue.empty() || n.queue.front() == v.id);

                if (canProceed)
                {
                    // If moving from a waiting state, pop self from queue
                    if (v.waiting && !n.queue.empty() && n.queue.front() == v.id)
                    {
                        n.queue.pop();
                    }
                    v.waiting = false;
                    v.waitingNode = -1;
                }
                else // Must wait (Red light OR stuck behind another car)
                {
                    if (!v.waiting)
                    {
                        // Join the queue for the first time
                        n.queue.push(v.id);
                        v.waiting = true;
                        v.waitingNode = curNode;
                    }
                    continue; // Stop movement and wait
                }
            }
        }

        // --- Move along segment ---
        Vec2f pcur = G.nodes[curNode].pos;
        Vec2f pnext = G.nodes[nextNode].pos;
        float segLen = distf(pcur, pnext);

        if (segLen < 1e-4f) // Segment length is zero (shouldn't happen in grid)
        {
            v.currentSegment++;
            v.posAlong = 0.0f;
            if (v.currentSegment >= v.path.size() - 1)
                v.finished = true;
            continue;
        }

        // Calculate normalized distance change
        float delta = (v.speed * dt) / segLen;
        v.posAlong += delta;

        if (v.posAlong >= 1.0f) // Segment finished
        {
            v.currentSegment++;
            v.posAlong = 0.0f;
            if (v.currentSegment >= v.path.size() - 1)
            {
                v.finished = true; // Destination reached
            }
        }
    }
}

/**
 * @brief Updates the position and traffic logic for NPC vehicles.
 */
void updateNPCs(float dt)
{
    for (int i = 0; i < npcVehicles.size(); ++i)
    {
        NPCVehicle &n = npcVehicles[i];
        if (n.currentNode < 0 || n.currentNode >= G.size() || n.nextNode < 0 || n.nextNode >= G.size())
            continue;

        Node &curNodeObj = G.nodes[n.currentNode];

        // --- Intersection Arrival/Waiting Logic ---
        if (n.posAlong <= 0.001f) // Just arrived at the intersection
        {
            bool green = (curNodeObj.lightState == 0);

            // NPC proceed condition: Green AND (Queue is empty OR NPC is at the front)
            bool canProceed = green && (curNodeObj.queue.empty() || curNodeObj.queue.front() == n.id);

            if (canProceed)
            {
                // If it was waiting, pop self from queue
                if (n.waiting && !curNodeObj.queue.empty() && curNodeObj.queue.front() == n.id)
                {
                    curNodeObj.queue.pop();
                }
                n.waiting = false; // Reset waiting state
                n.waitingNode = -1;
                // Fall through to movement logic
            }
            else // Must wait: Red OR not at the front of the queue
            {
                if (!n.waiting) // ONLY PUSH IF NOT ALREADY WAITING
                {
                    if (curNodeObj.queue.cap > 0)
                    {
                        curNodeObj.queue.push(n.id);
                        n.waiting = true;
                        n.waitingNode = n.currentNode;
                    }
                }
                continue; // Stop movement and wait
            }
        }

        // --- Move along segment ---
        Vec2f a = G.nodes[n.currentNode].pos;
        Vec2f b = G.nodes[n.nextNode].pos;
        float segLen = distf(a, b);

        if (segLen < 1e-4f) // Segment length is zero
        {
            n.currentNode = n.nextNode;
            n.posAlong = 0.0f;
            // Choose a new random destination node
            DynArr<Edge> &nb = G.adj[n.currentNode];
            n.nextNode = (nb.size() > 0) ? nb[GetRandomValue(0, nb.size() - 1)].to : n.currentNode;
            continue;
        }

        // Calculate normalized distance change
        float delta = (n.speed * dt) / segLen;
        n.posAlong += delta;

        if (n.posAlong >= 1.0f) // Segment finished
        {
            n.currentNode = n.nextNode;
            n.posAlong = 0.0f;
            // Choose a new random destination node
            DynArr<Edge> &nb = G.adj[n.currentNode];
            n.nextNode = (nb.size() > 0) ? nb[GetRandomValue(0, nb.size() - 1)].to : n.currentNode;
        }
    }
}

// -----------------------------
// Drawing helpers
// -----------------------------
/**
 * @brief Draws a single graph node (intersection) with traffic light status.
 */
void drawNode(int i)
{
    Node &n = G.nodes[i];
    Color col = LIGHTGRAY;
    // Highlight based on light state
    if (n.lightState == 0)
        col = (Color){200, 230, 200, 255}; // Greenish
    else
        col = (Color){230, 230, 200, 255}; // Reddish

    DrawCircle((int)n.pos.x, (int)n.pos.y, 8, col);
    DrawCircleLines((int)n.pos.x, (int)n.pos.y, 8, DARKGRAY);
    // Draw traffic light indicator
    Color lcol = (n.lightState == 0) ? DARKGREEN : RED;
    DrawRectangle((int)n.pos.x + 10, (int)n.pos.y - 6, 8, 12, lcol);
    if (showIds)
    {
        char buf[32];
        sprintf(buf, "%d", n.id);
        DrawText(buf, (int)n.pos.x + 18, (int)n.pos.y - 8, 10, RAYWHITE);
    }
}

/**
 * @brief Draws a single edge (road segment).
 */
void drawEdge(int a, int ei)
{
    Edge &e = G.adj[a][ei];
    Vec2f pa = G.nodes[a].pos;
    Vec2f pb = G.nodes[e.to].pos;
    DrawLineEx((Vector2){pa.x, pa.y}, (Vector2){pb.x, pb.y}, 2.0f, LIGHTGRAY);
}

/**
 * @brief Draws the calculated path (lastPath).
 */
void drawPath(const DynArr<int> &path)
{
    if (path.size() < 2)
        return;
    // Draw thick blue line for the path segments
    for (int i = 0; i < path.size() - 1; ++i)
    {
        Vec2f a = G.nodes[path[i]].pos;
        Vec2f b = G.nodes[path[i + 1]].pos;
        DrawLineEx((Vector2){a.x, a.y}, (Vector2){b.x, b.y}, 4.0f, (Color){100, 180, 255, 180});
    }
    // Highlight nodes along the path
    for (int i = 0; i < path.size(); ++i)
    {
        Vec2f p = G.nodes[path[i]].pos;
        DrawCircle((int)p.x, (int)p.y, 5, BLUE);
    }
}

/**
 * @brief Draws user-controlled vehicles.
 */
void drawVehicles()
{
    for (int i = 0; i < vehicles.size(); ++i)
    {
        Vehicle &v = vehicles[i];
        if (v.finished || v.path.size() < 2 || v.currentSegment >= v.path.size() - 1)
            continue;

        // Calculate vehicle position using LERP
        int s = v.currentSegment;
        Vec2f a = G.nodes[v.path[s]].pos;
        Vec2f b = G.nodes[v.path[s + 1]].pos;
        float t = v.posAlong;
        float x = a.x + (b.x - a.x) * t;
        float y = a.y + (b.y - a.y) * t;

        // Draw normal (Dark Blue) or emergency (Orange) vehicle
        Color c = (v.type == EMERGENCY) ? ORANGE : DARKBLUE;
        DrawCircle((int)x, (int)y, 6, c);
        DrawCircleLines((int)x, (int)y, 6, BLACK);

        // Draw vehicle ID
        char buf[32];
        sprintf(buf, "%d", v.id);
        DrawText(buf, (int)x + 8, (int)y - 8, 10, RAYWHITE);
    }
}

/**
 * @brief Draws NPC vehicles.
 */
void drawNPCs()
{
    for (int i = 0; i < npcVehicles.size(); ++i)
    {
        NPCVehicle &n = npcVehicles[i];
        if (n.currentNode < 0 || n.currentNode >= G.size() || n.nextNode < 0 || n.nextNode >= G.size())
            continue;

        // Calculate NPC position using LERP
        Vec2f a = G.nodes[n.currentNode].pos;
        Vec2f b = G.nodes[n.nextNode].pos;
        Vec2f p = lerp(a, b, n.posAlong);

        // Draw NPC (Green)
        DrawCircle((int)p.x, (int)p.y, 5, GREEN);
        DrawCircleLines((int)p.x, (int)p.y, 6, BLACK);
    }
}

/**
 * @brief Draws the status bar and controls information at the bottom of the screen.
 */
void drawUI()
{
    // ... (UI drawing logic using DrawText and DrawRectangle)
    DrawRectangle(0, screenH - 140, screenW, 140, (Color){20, 20, 20, 200});
    DrawText("SmartTraffic - Keyboard Shortcuts:", 12, screenH - 135, 26, SKYBLUE);
    DrawText("Left-click nodes to set Start then End", 12, screenH - 100, 18, GOLD);

    DrawText("D = ", 12, screenH - 74, 18, PINK);
    DrawText("Dijkstra", 45, screenH - 74, 18, RAYWHITE);
    DrawText("    A = ", 12 + 100, screenH - 74, 18, PURPLE);
    DrawText("A*", 12 + 155, screenH - 74, 18, RAYWHITE);
    DrawText("    V = ", 12 + 185, screenH - 74, 18, BLUE);
    DrawText("spawn vehicle", 12 + 245, screenH - 74, 18, RAYWHITE);
    DrawText("    E = ", 12 + 360, screenH - 74, 18, ORANGE);
    DrawText("spawn emergency", 12 + 410, screenH - 74, 18, RAYWHITE);
    DrawText("UP/DOWN = spawn/remove NPCs", 12 + 580, screenH - 74, 18, LIGHTGRAY);

    DrawText("S = ", 12, screenH - 55, 18, GREEN);
    DrawText("Start/Stop sim", 45, screenH - 55, 18, RAYWHITE);
    DrawText("    R = ", 12 + 165, screenH - 55, 18, RED);
    DrawText("Reset vehicles", 12 + 220, screenH - 55, 18, RAYWHITE);
    DrawText("    G = ", 12 + 345, screenH - 55, 18, YELLOW);
    DrawText("toggle node ids", 12 + 400, screenH - 55, 18, RAYWHITE);

    char buf[256];
    sprintf(buf, "Start: %d   End: %d   Path found: %s   Algorithm: %s   Vehicles: %d   NPCs: %d",
            selectedStart, selectedEnd, pathFound ? "YES" : "NO", useAstar ? "A*" : "Dijkstra", vehicles.size(), npcVehicles.size());
    DrawText(buf, 12, screenH - 25, 18, BEIGE);

    DrawText("Legend:", screenW - 350, screenH - 130, 26, SKYBLUE);
    DrawRectangle(screenW - 270, screenH - 130, 260, 120, (Color){40, 40, 40, 180});
    DrawText("Node:", screenW - 260, screenH - 125, 14, LIGHTGRAY);
    DrawCircle(screenW - 105, screenH - 117, 6, LIGHTGRAY);
    DrawText("Edge:", screenW - 260, screenH - 105, 14, LIGHTGRAY);
    DrawLine(screenW - 120, screenH - 97, screenW - 90, screenH - 97, LIGHTGRAY);
    DrawText("Path:", screenW - 260, screenH - 85, 14, LIGHTGRAY);
    DrawRectangle(screenW - 120, screenH - 80, 30, 6, (Color){100, 180, 255, 180});
    DrawText("Normal Vehicle:", screenW - 260, screenH - 65, 14, LIGHTGRAY);
    DrawCircle(screenW - 105, screenH - 60, 6, DARKBLUE);
    DrawText("Emergency Vehicle:", screenW - 260, screenH - 45, 14, LIGHTGRAY);
    DrawCircle(screenW - 105, screenH - 40, 6, ORANGE);
    DrawText("NPC Vehicle:", screenW - 260, screenH - 25, 14, LIGHTGRAY);
    DrawCircle(screenW - 105, screenH - 20, 6, GREEN);

    if (selectedStart != -1)
    {
        Vec2f p = G.nodes[selectedStart].pos;
        DrawCircleLines((int)p.x, (int)p.y, 14, GREEN);
        DrawText("Start", (int)p.x - 18, (int)p.y - 28, 10, GREEN);
    }
    if (selectedEnd != -1)
    {
        Vec2f p = G.nodes[selectedEnd].pos;
        DrawCircleLines((int)p.x, (int)p.y, 14, RED);
        DrawText("End", (int)p.x - 12, (int)p.y - 28, 10, RED);
    }
}

// -----------------------------
// Main
// -----------------------------
int main()
{
    // ... (Initialization)
    Image icon = LoadImage("./assets/icon.png");
    srand((unsigned)time(NULL));
    InitWindow(screenW, screenH, "SmartTraffic - Data Structures Project Demo");
    SetWindowIcon(icon);
    UnloadImage(icon);
    SetTargetFPS(60);

    // Setup grid parameters and build graph
    int cols = 12, rows = 8;
    float marginX = 80, marginY = 40;
    float gx = (screenW - marginX * 2) / (float)(cols - 1);
    float gy = (screenH - 180 - marginY * 2) / (float)(rows - 1);
    buildGridGraph(cols, rows, marginX, marginY, gx, gy);

    // Initialize pathfinding arrays
    g_scores.clear();
    cameFrom.clear();
    ensureArrSize(G.size());

    vehicles.clear();
    npcVehicles.clear();

    // Spawn initial NPCs
    spawnNPCs(6);

    // Main game loop
    double prev = GetTime();
    while (!WindowShouldClose())
    {
        double cur = GetTime();
        float dt = (float)(cur - prev);
        prev = cur;

        // --- Input Handling ---
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
        {
            // Select start and end nodes based on mouse click
            Vector2 mp = GetMousePosition();
            int nearest = -1;
            float ndist = 1e9f;
            for (int i = 0; i < G.size(); ++i)
            {
                float d = distf(G.nodes[i].pos, Vec2f(mp.x, mp.y));
                if (d < ndist)
                {
                    ndist = d;
                    nearest = i;
                }
            }
            if (nearest != -1 && ndist < 16.0f)
            {
                if (selectedStart == -1)
                    selectedStart = nearest;
                else if (selectedEnd == -1 && nearest != selectedStart)
                    selectedEnd = nearest;
                else
                {
                    selectedStart = nearest;
                    selectedEnd = -1;
                }
            }
        }

        // Run Dijkstra (D)
        if (IsKeyPressed(KEY_D))
        {
            if (selectedStart != -1 && selectedEnd != -1)
            {
                useAstar = false;
                runDijkstra(selectedStart, selectedEnd);
            }
        }

        // Run A* (A)
        if (IsKeyPressed(KEY_A))
        {
            if (selectedStart != -1 && selectedEnd != -1)
            {
                useAstar = true;
                runAstar(selectedStart, selectedEnd);
            }
        }

        // Spawn Normal Vehicle (V)
        if (IsKeyPressed(KEY_V))
        {
            if (selectedStart != -1 && selectedEnd != -1)
                spawnVehicle(selectedStart, selectedEnd, NORMAL);
        }

        // Spawn Emergency Vehicle (E)
        if (IsKeyPressed(KEY_E))
        {
            if (selectedStart != -1 && selectedEnd != -1)
                spawnVehicle(selectedStart, selectedEnd, EMERGENCY);
        }

        // Toggle Simulation (S)
        if (IsKeyPressed(KEY_S))
            runningSim = !runningSim;

        // Reset (R)
        if (IsKeyPressed(KEY_R))
        {
            resetVehicles();
            selectedStart = selectedEnd = -1;
            pathFound = false;
            lastPath.clear();
        }

        // Toggle Node IDs (G)
        if (IsKeyPressed(KEY_G))
            showIds = !showIds;

        // NPC spawn/remove (UP/DOWN)
        if (IsKeyPressed(KEY_UP))
            spawnNPCs(1);
        if (IsKeyPressed(KEY_DOWN))
            removeNPCs(1);

        // --- Simulation Update ---
        if (runningSim)
        {
            updateTrafficLights(dt);
            updateVehicles(dt);
            updateNPCs(dt);
        }

        // --- Drawing ---
        BeginDrawing();
        ClearBackground((Color){30, 30, 40, 255});

        // Draw Edges
        for (int i = 0; i < G.size(); ++i)
        {
            for (int j = 0; j < G.adj[i].size(); ++j)
                drawEdge(i, j);
        }

        // Draw Path
        if (pathFound)
            drawPath(lastPath);

        // Draw Nodes
        for (int i = 0; i < G.size(); ++i)
            drawNode(i);

        // Draw NPCs and Vehicles
        drawNPCs(); // Draw NPCs first (under user vehicles)
        drawVehicles();

        // Draw UI and HUD
        drawUI();

        EndDrawing();
    }

    // Cleanup
    CloseWindow();
    return 0;
}