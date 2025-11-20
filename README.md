# 🚦 Simulated Traffic System

_A C++ + Raylib simulation demonstrating pathfinding and traffic management._

---

## 📌 Overview

This project is a **Smart Traffic System** built in **C++** using **custom data structures** and **Raylib** for visualization.
It simulates a small city map, traffic lights, and vehicle flow — and allows users to compute the **shortest path** between two points using:

-   **Dijkstra's Algorithm**
-   **A\* Search Algorithm**

This acts as a simplified “Google Maps–style” routing demo while showcasing core Data Structures.

---

## 🎯 Features

### 🛣️ **Pathfinding**

-   Shortest/fastest route between two locations
-   Visual node + edge animations
-   Supports Dijkstra and A\*

### 🚗 **Traffic Simulation**

-   Vehicles moving along roads
-   Traffic lights with cycles
-   Emergency vehicle priority system

### 🏗️ **Data Structures Used**

| Structure                  | Purpose                                            |
| -------------------------- | -------------------------------------------------- |
| **Graph (Adjacency List)** | Represents city roads & intersections              |
| **Priority Queue (Heap)**  | Dijkstra/A\* frontier + emergency vehicle priority |
| **Queue**                  | Normal traffic flow simulation                     |
| **Hash Map**               | Vehicle lookup and state management                |
| **Array**                  | Traffic light cycles & vehicle states              |

All structures are **implemented manually** (no STL containers for these features).

---

## 🖥️ Screens & UI

-   Map view
-   Start & destination selection
-   Visualized traversal
-   Legends for signals and nodes
-   Clean UI using Raylib (no external UI frameworks)

---

## ⚙️ Installation

### 1. Install Raylib (Windows)

-   Install Raylib for windows from (Raylib.com)[raylib.com]

-   Make sure installation is in your `C:\raylib`

-   Clone this repo

### 2. Build & Run

-   Open folder in VSCode

-   Press: `F5` to Build + Run

---

## 📚 How It Works

-   The city is modeled as a graph of intersections (nodes) and roads (edges).
-   When the user selects two points, the system runs Dijkstra/A\*.
-   Real-time animation highlights the search frontier and chosen path.
-   Traffic lights dynamically control flow.
-   Emergency vehicles are prioritized using a **custom max-heap**.
