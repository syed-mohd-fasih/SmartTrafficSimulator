#include "raylib.h"

#include <cmath>
#include <cstdio>
#include <ctime>

// -----------------------------
// Utility Vec2f
// -----------------------------
struct Vec2f
{
    float x, y;
    Vec2f() : x(0), y(0) {}
    Vec2f(float _x, float _y) : x(_x), y(_y) {}
};

static float distf(const Vec2f &a, const Vec2f &b)
{
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    return sqrtf(dx * dx + dy * dy);
}

// -----------------------------
// DynArr: minimal dynamic array (with deep copy)
// -----------------------------
template <typename T>
struct DynArr
{
    T *data;
    int count;
    int capacity;

    DynArr() : data(nullptr), count(0), capacity(0) {}

    // Deep copy constructor
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

    // Deep copy assignment
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
    int size() const { return count; }

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

    void push(const T &v)
    {
        if (count + 1 > capacity)
            reserve(count + 1);
        data[count++] = v;
    }

    T &operator[](int i) { return data[i]; }
    const T &operator[](int i) const { return data[i]; }
};

// -----------------------------
// SimpleQueue: circular queue using raw array
// -----------------------------
struct SimpleQueueInt
{
    int *buf;
    int cap;
    int head, tail;

    SimpleQueueInt() : buf(nullptr), cap(0), head(0), tail(0) {}
    ~SimpleQueueInt() { delete[] buf; }

    bool empty() const { return head == tail; }

    void reserve(int n)
    {
        if (n <= cap)
            return;

        int newcap = cap == 0 ? 8 : cap;
        while (newcap < n)
            newcap <<= 1;

        int *nb = new int[newcap];

        // copy existing elements in order (head -> tail)
        int sz = 0;
        if (buf != nullptr && cap > 0 && !empty())
        {
            int idx = head;
            while (idx != tail)
            {
                nb[sz++] = buf[idx];
                idx = (idx + 1) % cap; // ✅ FIXED: step by 1, not +cap
            }
        }

        delete[] buf;
        buf = nb;
        cap = newcap;
        head = 0;
        tail = sz;
    }

    void push(int v)
    {
        if (cap == 0)
            reserve(8);

        int next = (tail + 1) % cap;
        if (next == head)
        {
            // expand, but keep existing elements
            reserve(cap * 2);
            next = (tail + 1) % cap;
        }
        buf[tail] = v;
        tail = next;
    }

    int pop()
    {
        int v = buf[head];
        head = (head + 1) % cap;
        return v;
    }

    // peek front (undefined if empty)
    int front() const { return buf[head]; }
};

// -----------------------------
// MinHeap for (id,priority)
// -----------------------------
struct HeapNode
{
    int id;
    float pr;
};
struct MinHeap
{
    HeapNode *buf;
    int cap;
    int len;
    MinHeap() : buf(nullptr), cap(0), len(0) {}
    ~MinHeap() { delete[] buf; }
    void reserve(int n)
    {
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
    HeapNode pop()
    {
        HeapNode out = buf[0];
        buf[0] = buf[len - 1];
        len--;
        int i = 0;
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
struct Edge
{
    int to;
    float w;
};
struct Node
{
    Vec2f pos;
    int id;
    float lightTimer;
    int lightState;
    SimpleQueueInt queue;
    Node() : pos(), id(-1), lightTimer(0), lightState(0) {}
};

struct Graph
{
    DynArr<Node> nodes;
    DynArr<DynArr<Edge>> adj;
    Graph() {}
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
    EMERGENCY = 1
};

struct Vehicle
{
    int id;
    VehicleType type;
    int src, dest;
    DynArr<int> path;
    int currentSegment;
    float posAlong;
    float speed;
    bool finished;

    // queue-related state (to avoid multiple pushes)
    bool waiting;
    int waitingNode;

    Vehicle()
        : id(-1), type(NORMAL), src(-1), dest(-1),
          currentSegment(0), posAlong(0), speed(80.0f),
          finished(false), waiting(false), waitingNode(-1) {}
};

// -----------------------------
// Globals
// -----------------------------
static Graph G;
static DynArr<Vehicle> vehicles;
static MinHeap emergencyHeap; // store emergency vehicle ids (not heavily used)
static int nextVehicleId = 1;
static int selectedStart = -1, selectedEnd = -1;
static bool showIds = true;
static bool runningSim = false;
static bool pathFound = false;
static DynArr<int> lastPath;
static bool useAstar = false;
static float lightCycle = 6.0f;
static int screenW = 1280, screenH = 720;

// -----------------------------
// Build simple grid graph
// -----------------------------
void buildGridGraph(int cols, int rows, float marginX, float marginY, float gx, float gy)
{
    G.init(0);
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
    auto idx = [&](int c, int r)
    { return r * cols + c; };

    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            int a = idx(c, r);
            if (c + 1 < cols)
            {
                int b = idx(c + 1, r);
                float w = distf(G.nodes[a].pos, G.nodes[b].pos);
                G.addEdge(a, b, w);
                G.addEdge(b, a, w);
            }
            if (r + 1 < rows)
            {
                int b = idx(c, r + 1);
                float w = distf(G.nodes[a].pos, G.nodes[b].pos);
                G.addEdge(a, b, w);
                G.addEdge(b, a, w);
            }
        }
    }
    for (int i = 0; i < G.size(); ++i)
    {
        G.nodes[i].lightTimer = (float)(i % 3);
        G.nodes[i].lightState = 0;
        G.nodes[i].queue.reserve(8);
    }
}

// -----------------------------
// Pathfinding (Dijkstra + A*)
// -----------------------------
static DynArr<float> g_scores;
static DynArr<int> cameFrom;

void ensureArrSize(int n)
{
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

bool reconstructPath(int s, int t, DynArr<int> &out)
{
    out.clear();
    if (cameFrom[t] == -1 && s != t)
        return false;
    int cur = t;
    while (cur != -1)
    {
        out.push(cur);
        if (cur == s)
            break;
        cur = cameFrom[cur];
    }
    // reverse
    for (int i = 0, j = out.size() - 1; i < j; ++i, --j)
    {
        int tmp = out[i];
        out[i] = out[j];
        out[j] = tmp;
    }
    return true;
}

bool runDijkstra(int s, int t)
{
    int n = G.size();
    ensureArrSize(n);
    MinHeap heap;
    heap.reserve(n + 4);
    g_scores[s] = 0.0f;
    heap.push(s, 0.0f);
    while (!heap.empty())
    {
        HeapNode hn = heap.pop();
        int u = hn.id;
        float d = hn.pr;
        if (d != g_scores[u])
            continue;
        if (u == t)
        {
            bool ok = reconstructPath(s, t, lastPath);
            pathFound = ok;
            return ok;
        }
        for (int ei = 0; ei < G.adj[u].size(); ++ei)
        {
            Edge &e = G.adj[u][ei];
            float nd = d + e.w;
            if (nd < g_scores[e.to])
            {
                g_scores[e.to] = nd;
                cameFrom[e.to] = u;
                heap.push(e.to, nd);
            }
        }
    }
    pathFound = false;
    return false;
}

bool runAstar(int s, int t)
{
    int n = G.size();
    ensureArrSize(n);
    MinHeap heap;
    heap.reserve(n + 4);
    DynArr<float> fscore;
    fscore.clear();
    for (int i = 0; i < n; i++)
        fscore.push(INFINITY);
    g_scores[s] = 0.0f;
    fscore[s] = distf(G.nodes[s].pos, G.nodes[t].pos);
    heap.push(s, fscore[s]);
    while (!heap.empty())
    {
        HeapNode hn = heap.pop();
        int u = hn.id;
        if (u == t)
        {
            bool ok = reconstructPath(s, t, lastPath);
            pathFound = ok;
            return ok;
        }
        for (int ei = 0; ei < G.adj[u].size(); ++ei)
        {
            Edge &e = G.adj[u][ei];
            float tentative = g_scores[u] + e.w;
            if (tentative < g_scores[e.to])
            {
                cameFrom[e.to] = u;
                g_scores[e.to] = tentative;
                float est = tentative + distf(G.nodes[e.to].pos, G.nodes[t].pos);
                heap.push(e.to, est);
            }
        }
    }
    pathFound = false;
    return false;
}

// -----------------------------
// Vehicles management
// -----------------------------
Vehicle *spawnVehicle(int src, int dest, VehicleType type)
{
    if (src < 0 || dest < 0)
        return nullptr;
    // compute path first (Dijkstra)
    if (!pathFound)
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

void resetVehicles()
{
    vehicles.clear();
    emergencyHeap.clear();
    nextVehicleId = 1;
}

// -----------------------------
// Simulation updates
// -----------------------------
void updateTrafficLights(float dt)
{
    for (int i = 0; i < G.size(); ++i)
    {
        Node &n = G.nodes[i];
        n.lightTimer += dt;
        if (n.lightTimer >= lightCycle)
            n.lightTimer -= lightCycle;
        n.lightState = (n.lightTimer < (lightCycle * 0.5f)) ? 0 : 1;

        // emergency override: if any emergency in queue, force state 0 briefly
        if (n.queue.cap > 0 && !n.queue.empty())
        {
            int head = n.queue.head;
            int tail = n.queue.tail;
            int idx = head;
            bool found = false;

            while (idx != tail)
            {
                int vid = n.queue.buf[idx];

                // LINEAR SEARCH for vehicle
                Vehicle *pv = nullptr;
                for (int k = 0; k < vehicles.size(); ++k)
                {
                    if (vehicles[k].id == vid)
                    {
                        pv = &vehicles[k];
                        break;
                    }
                }

                if (pv && pv->type == EMERGENCY)
                {
                    found = true;
                    break;
                }
                idx = (idx + 1) % n.queue.cap;
            }
            if (found)
                n.lightState = 0;
        }
    }
}

void updateVehicles(float dt)
{
    for (int i = 0; i < vehicles.size(); ++i)
    {
        Vehicle &v = vehicles[i];
        if (v.finished)
            continue;
        if (v.path.size() < 2)
        {
            v.finished = true;
            continue;
        }

        // safety: bounds check on segment index
        if (v.currentSegment >= v.path.size() - 1)
        {
            v.finished = true;
            continue;
        }

        int curNode = v.path[v.currentSegment];
        int nextNode = v.path[v.currentSegment + 1];

        if (curNode < 0 || curNode >= G.size() || nextNode < 0 || nextNode >= G.size())
        {
            v.finished = true;
            continue;
        }

        if (v.posAlong <= 0.001f)
        {
            Node &n = G.nodes[curNode];
            bool green = (n.lightState == 0);

            // EMERGENCY BYPASS
            if (v.type == EMERGENCY)
            {
                // skip lights and queue
                v.waiting = false;
                v.waitingNode = -1;

                // clear queue at this node
                while (!n.queue.empty())
                    n.queue.pop();
            }
            else
            {
                // NORMAL logic
                if (!green)
                {
                    n.queue.push(v.id);
                    v.waiting = true;
                    v.waitingNode = curNode;
                    continue;
                }

                if (!n.queue.empty())
                {
                    int front = n.queue.front();
                    if (front != v.id)
                        continue;
                    n.queue.pop();
                    v.waiting = false;
                }
            }

            {
                bool allowed = true;
                if (!n.queue.empty())
                {
                    int front = n.queue.front();
                    if (front != v.id)
                        allowed = false;
                    else
                    {
                        n.queue.pop();
                        v.waiting = false;
                        v.waitingNode = -1;
                    }
                }
                if (!allowed)
                    continue;
            }
        }

        Vec2f pcur = G.nodes[curNode].pos;
        Vec2f pnext = G.nodes[nextNode].pos;
        float segLen = distf(pcur, pnext);
        if (segLen < 1e-4f)
        {
            v.currentSegment++;
            v.posAlong = 0.0f;
            if (v.currentSegment >= v.path.size() - 1)
            {
                v.finished = true;
                continue;
            }
            continue;
        }
        float delta = (v.speed * dt) / segLen;
        v.posAlong += delta;
        if (v.posAlong >= 1.0f)
        {
            v.currentSegment++;
            v.posAlong = 0.0f;
            if (v.currentSegment >= v.path.size() - 1)
            {
                v.finished = true;
                continue;
            }
        }
    }
}

// -----------------------------
// Drawing helpers
// -----------------------------
void drawNode(int i)
{
    Node &n = G.nodes[i];
    Color col = LIGHTGRAY;
    if (n.lightState == 0)
        col = (Color){200, 230, 200, 255};
    else
        col = (Color){230, 230, 200, 255};
    DrawCircle((int)n.pos.x, (int)n.pos.y, 8, col);
    DrawCircleLines((int)n.pos.x, (int)n.pos.y, 8, DARKGRAY);
    Color lcol = (n.lightState == 0) ? GREEN : RED;
    DrawRectangle((int)n.pos.x + 10, (int)n.pos.y - 6, 8, 12, lcol);
    if (showIds)
    {
        char buf[32];
        sprintf(buf, "%d", n.id);
        DrawText(buf, (int)n.pos.x + 18, (int)n.pos.y - 8, 10, RAYWHITE);
    }
}

void drawEdge(int a, int ei)
{
    Edge &e = G.adj[a][ei];
    Vec2f pa = G.nodes[a].pos;
    Vec2f pb = G.nodes[e.to].pos;
    DrawLineEx((Vector2){pa.x, pa.y}, (Vector2){pb.x, pb.y}, 2.0f, LIGHTGRAY);
}

void drawPath(const DynArr<int> &path)
{
    if (path.size() < 2)
        return;
    for (int i = 0; i < path.size() - 1; ++i)
    {
        Vec2f a = G.nodes[path[i]].pos;
        Vec2f b = G.nodes[path[i + 1]].pos;
        DrawLineEx((Vector2){a.x, a.y}, (Vector2){b.x, b.y}, 4.0f, (Color){100, 180, 255, 180});
    }
    for (int i = 0; i < path.size(); ++i)
    {
        Vec2f p = G.nodes[path[i]].pos;
        DrawCircle((int)p.x, (int)p.y, 5, BLUE);
    }
}

void drawVehicles()
{
    for (int i = 0; i < vehicles.size(); ++i)
    {
        Vehicle &v = vehicles[i];
        if (v.finished)
            continue;
        if (v.path.size() < 2)
            continue;
        if (v.currentSegment >= v.path.size() - 1)
            continue;

        int s = v.currentSegment;
        Vec2f a = G.nodes[v.path[s]].pos;
        Vec2f b = G.nodes[v.path[s + 1]].pos;
        float t = v.posAlong;
        float x = a.x + (b.x - a.x) * t;
        float y = a.y + (b.y - a.y) * t;
        Color c = (v.type == EMERGENCY) ? ORANGE : DARKBLUE;
        DrawCircle((int)x, (int)y, 6, c);
        DrawCircleLines((int)x, (int)y, 6, BLACK);
        char buf[32];
        sprintf(buf, "%d", v.id);
        DrawText(buf, (int)x + 8, (int)y - 8, 10, RAYWHITE);
    }
}

void drawUI()
{
    DrawRectangle(0, screenH - 120, screenW, 120, (Color){20, 20, 20, 200});
    DrawText("SmartTraffic - Keyboard Shortcuts:", 12, screenH - 110, 26, RAYWHITE);
    DrawText("Left-click nodes to set Start then End", 12, screenH - 83, 18, RAYWHITE);

    DrawText("D = ", 12, screenH - 64, 18, BLUE);
    DrawText("Dijkstra", 45, screenH - 64, 18, RAYWHITE);
    DrawText("    A = ", 12 + 100, screenH - 64, 18, GREEN);
    DrawText("A*", 12 + 155, screenH - 64, 18, RAYWHITE);
    DrawText("    V = ", 12 + 180, screenH - 64, 18, BLUE);
    DrawText("spawn vehicle", 12 + 245, screenH - 64, 18, RAYWHITE);
    DrawText("    E = ", 12 + 370, screenH - 64, 18, ORANGE);
    DrawText("spawn emergency", 12 + 435, screenH - 64, 18, RAYWHITE);

    DrawText("S = ", 12, screenH - 45, 18, GREEN);
    DrawText("Start/Stop sim", 45, screenH - 45, 18, RAYWHITE);
    DrawText("    R = ", 12 + 150, screenH - 45, 18, RED);
    DrawText("Reset vehicles", 12 + 205, screenH - 45, 18, RAYWHITE);
    DrawText("    G = ", 12 + 330, screenH - 45, 18, YELLOW);
    DrawText("toggle node ids", 12 + 385, screenH - 45, 18, RAYWHITE);

    char buf[256];
    sprintf(buf, "Start: %d   End: %d   Path found: %s   Algorithm: %s   Vehicles: %d",
            selectedStart, selectedEnd, pathFound ? "YES" : "NO", useAstar ? "A*" : "Dijkstra", vehicles.size());
    DrawText(buf, 12, screenH - 25, 18, RAYWHITE);

    DrawText("Legend:", screenW - 350, screenH - 110, 20, RAYWHITE);
    DrawRectangle(screenW - 270, screenH - 110, 260, 100, (Color){40, 40, 40, 180});
    DrawText("Node:", screenW - 260, screenH - 105, 14, RAYWHITE);
    DrawCircle(screenW - 105, screenH - 97, 6, LIGHTGRAY);
    DrawText("Edge:", screenW - 260, screenH - 85, 14, RAYWHITE);
    DrawLine(screenW - 120, screenH - 77, screenW - 90, screenH - 77, LIGHTGRAY);
    DrawText("Path:", screenW - 260, screenH - 65, 14, RAYWHITE);
    DrawRectangle(screenW - 120, screenH - 60, 30, 6, (Color){100, 180, 255, 180});
    DrawText("Normal Vehicle:", screenW - 260, screenH - 45, 14, RAYWHITE);
    DrawCircle(screenW - 105, screenH - 40, 6, DARKBLUE);
    DrawText("Emergency Vehicle:", screenW - 260, screenH - 25, 14, RAYWHITE);
    DrawCircle(screenW - 105, screenH - 20, 6, ORANGE);
}

// -----------------------------
// Main
// -----------------------------
int main()
{
    Image icon = LoadImage("./assets/icon.png");
    srand((unsigned)time(NULL));
    InitWindow(screenW, screenH, "SmartTraffic - Data Structures Project Demo");
    SetWindowIcon(icon);
    UnloadImage(icon);
    SetTargetFPS(60);

    int cols = 12, rows = 8;
    float marginX = 80, marginY = 40;
    float gx = (screenW - marginX * 2) / (float)(cols - 1);
    float gy = (screenH - 180 - marginY * 2) / (float)(rows - 1);
    buildGridGraph(cols, rows, marginX, marginY, gx, gy);

    g_scores.clear();
    cameFrom.clear();
    ensureArrSize(G.size());

    vehicles.clear();

    double prev = GetTime();
    while (!WindowShouldClose())
    {
        double cur = GetTime();
        float dt = (float)(cur - prev);
        prev = cur;

        // Input
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
        {
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
        if (IsKeyPressed(KEY_D))
        {
            if (selectedStart != -1 && selectedEnd != -1)
            {
                useAstar = false;
                runDijkstra(selectedStart, selectedEnd);
            }
        }
        if (IsKeyPressed(KEY_A))
        {
            if (selectedStart != -1 && selectedEnd != -1)
            {
                useAstar = true;
                runAstar(selectedStart, selectedEnd);
            }
        }
        if (IsKeyPressed(KEY_V))
        {
            if (selectedStart != -1 && selectedEnd != -1)
                spawnVehicle(selectedStart, selectedEnd, NORMAL);
        }
        if (IsKeyPressed(KEY_E))
        {
            if (selectedStart != -1 && selectedEnd != -1)
                spawnVehicle(selectedStart, selectedEnd, EMERGENCY);
        }
        if (IsKeyPressed(KEY_S))
            runningSim = !runningSim;
        if (IsKeyPressed(KEY_R))
        {
            resetVehicles();
            selectedStart = selectedEnd = -1;
            pathFound = false;
            lastPath.clear();
        }
        if (IsKeyPressed(KEY_G))
            showIds = !showIds;

        if (runningSim)
        {
            updateTrafficLights(dt);
            updateVehicles(dt);
        }

        BeginDrawing();
        ClearBackground((Color){30, 30, 40, 255});

        for (int i = 0; i < G.size(); ++i)
        {
            for (int j = 0; j < G.adj[i].size(); ++j)
                drawEdge(i, j);
        }

        if (pathFound)
            drawPath(lastPath);
        for (int i = 0; i < G.size(); ++i)
            drawNode(i);
        drawVehicles();

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

        drawUI();
        EndDrawing();
    }

    return 0;
}
