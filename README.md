# Greeduce
A heuristic solver for the Hitting Set and Dominating Set problems submitted to PACE Challenge 2025

Authors: Adam Polak and Jonas Schmidt

## Installation

The solver is in a single C++ file with no external dependencies. In order to compile it, simply run make.

## Running

The solver reads a problem instance from stdin, in either .gr or .hgr format, as described at https://pacechallenge.org/2025

## Solver description

On a high level, the solver combines the standard best-bang-for-the-buck greedy algorithm with reduction rules. It runs the greedy algorithm multiple times, with randomized tie-breaking, until the time limit is reached.

It is primarily a Hitting Set solver. When given an instance of the Dominating Set problem, it simply converts it to an equivalent Hitting Set instance.

### Reduction rules

The solver uses the Unit (Hyper-)Edge Rule, (Hyper-)Edge Domination Rule, and Vertex Domination Rule, as described in https://epubs.siam.org/doi/epdf/10.1137/1.9781611977042.17

The rules are applied exhaustively, in a lazy manner: the solver keeps a queue of objects (vertices and hyperedges) that can be still potentially reduced, which is first initialized to contain the whole instance, and whenever and object is reduced all its incident objects are added to the queue (if they are not there at the moment).

The rules are implemented to run in time O(d^3) per object, where d is the maximum degree, to benefit from the fact that the PACE Challenge input instances tend to be sparse. In particular, the exhaustive application of the rules does not require iterating over all pairs of vertices, which would be prohibitively slow.

### Greedy algorithm

The greedy algorithm iteratively picks a vertex with the highest number of (yet unhit) hyperedges and adds it to the solution. If there are multiple such vertices, a random one is picked.

### Greedy + Reduce = Greeduce

The algorithm that we dub Greeduce combines the previous two ideas. It applies the reduction rules exhaustively, and when no further reduction is possible it (1) greedily picks a vertex, (2) adds it to the solution, (3) removes from the instance the vertex and incident hyperedges, and (4) adds other vertices incident to these hyperedges to the queue with candidates for reductions.

### Evolutionary improvement

Greeduce is very fast on sparse instances, so within the 5-minute PACE time limit the solver can run it multiple times and output the best solution. Since the greedy vertex selection breaks ties randomly, multiple runs would indeed result in different solutions. However, the solver tries to make a better use of the multiple runs. The greedy algorithm is modified to accept a hint, which is just a subset of vertices supposed to be close to a good solution. These vertices are always prioritized when making a greedy choice (with ties broken first by the current degree, and then by randomness). Each run of Greeduce uses as the hint the best solution found so far with some vertices removed at random.

### Increasing neighborhood limit for reductions

On high-degree vertices, the reduction rules can get very slow: even a single run of Greeduce can run out of the 5-minute time limit. In order to circumvent this issue we introduce a neighborhood limit so that the two domination rules only check neighborhoods of vertices/hyperedges up to that limit. The limit is set to k in the k-th run of Greeduce.