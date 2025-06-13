#include <signal.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <chrono>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

const bool DBG = false;  // print debug info to stderr

const int TIME_LIMIT_MS = 5 * 60 * 1000;
const int REDUCTION_TIME_LIMIT_MS = 0.5 * TIME_LIMIT_MS;

typedef vector<bool> Solution;

minstd_rand prg;

int num_objects, num_vertices, num_hyperedges;
vector<vector<int>> incidence;
vector<int> vertex_to_original_id, always_in_solution;

Solution best;
int best_cost;

bool volatile sigterm_received = false;

void sigterm_handler(int _signal) {
  sigterm_received = true;
}

auto start_time = std::chrono::system_clock::now();

int elapsed_time_ms() {
  return chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-start_time).count();
}

void check_signal() {
  if (sigterm_received || elapsed_time_ms() > TIME_LIMIT_MS) {
    cout << best_cost + always_in_solution.size() << endl;
    for (int i = 0; i < num_vertices; ++i)
      if (best[i])
        cout << vertex_to_original_id[i] << endl;
    for (int v : always_in_solution)
      cout << v << endl;
    exit(0);
  }
}

/////////////////// input //////////////////////////////////////////////////////

void read_input() {
  string line;
  for (; line.empty() || line[0] != 'p'; getline(cin, line));
  istringstream iss(line);
  string ignore, problem_type;
  iss >> ignore >> problem_type;
  if (problem_type == "hs") {
    iss >> num_vertices >> num_hyperedges;
    num_objects = num_vertices + num_hyperedges;
    incidence.resize(num_objects);
    vertex_to_original_id.resize(num_vertices);
    for (int i = 0; i < num_vertices; ++i)
      vertex_to_original_id[i] = i + 1;
    for (int i = 0; i < num_hyperedges; ++i) {
      getline(cin, line);
      istringstream edge_iss(line);
      int v;
      while (edge_iss >> v) {
        --v;
        incidence[v].push_back(num_vertices + i);
        incidence[num_vertices + i].push_back(v);
      }
    }
  } else if (problem_type == "ds") {
    int num_edges;
    iss >> num_vertices >> num_edges;
    num_hyperedges = num_vertices;
    num_objects = 2 * num_vertices;
    incidence.resize(num_objects);
    vertex_to_original_id.resize(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
      incidence[i].push_back(num_vertices + i);
      incidence[num_vertices + i].push_back(i);
      vertex_to_original_id[i] = i + 1;
    }
    for (int i = 0; i < num_edges; ++i) {
      int a, b;
      cin >> a >> b;
      assert(1 <= a && a <= num_vertices && 1 <= b && b <= num_vertices);
      --a;
      --b;
      incidence[a].push_back(num_vertices + b);
      incidence[b].push_back(num_vertices + a);
      incidence[num_vertices + a].push_back(b);
      incidence[num_vertices + b].push_back(a);
    }
  } else {
    assert(false);
  }
  for (int i = 0; i < num_objects; ++i)
    sort(incidence[i].begin(), incidence[i].end());
}

/////////////////// utils //////////////////////////////////////////////////////

inline bool sorted_vector_contains(const vector<int> &v, int x) {
  return binary_search(v.begin(), v.end(), x);
}

void prune_inactive_objects(const vector<bool>& active) {
  assert(active.size() == num_objects);
  int new_num_objects = 0;
  vector<int> new_object_id(num_objects, -1);
  for (int i = 0; i < num_objects; ++i)
    if (active[i])
      new_object_id[i] = new_num_objects++;
  int last_new_vertex_id = -1;
  for (int i = 0; i < num_objects; ++i) {
    if (active[i]) {
      vector<int> new_incidence_list;
      for (int j : incidence[i])
        if (active[j])
          new_incidence_list.push_back(new_object_id[j]);
      incidence[new_object_id[i]] = new_incidence_list;
      if (i < num_vertices) {
        vertex_to_original_id[new_object_id[i]] = vertex_to_original_id[i];
        last_new_vertex_id = new_object_id[i];
      }
    }
  }
  assert(last_new_vertex_id != -1);
  num_objects = new_num_objects;
  num_vertices = last_new_vertex_id + 1;
  num_hyperedges = num_objects - num_vertices;
  vertex_to_original_id.resize(num_vertices);
  incidence.resize(num_objects);
  if (DBG) cerr << "nodes always in solution " << always_in_solution.size() << endl;
}

void remove_unnecessary_nodes(Solution& solution) {
  vector<int> hitcount(num_hyperedges);
  for (int i = 0; i < num_vertices; ++i)
    if (solution[i])
      for (int h : incidence[i])
        ++hitcount[h - num_vertices];
  vector<bool> unnecessary = solution;
  for (int h = 0; h < num_hyperedges; ++h) {
    assert(hitcount[h] > 0);
    if (hitcount[h] == 1)
      for (int v : incidence[h + num_vertices])
        unnecessary[v] = false;
  }
  for (int v = 0; v < num_vertices; ++v) {
    if (unnecessary[v]) {
      solution[v] = false;
      for (int h : incidence[v]) {
        --hitcount[h - num_vertices];
        assert(hitcount[h - num_vertices] > 0);
        if (hitcount[h - num_vertices] == 1)
          for (int u : incidence[h])
            unnecessary[u] = false;
      }
    }
  }
}

/////////////////// greedy + reduce ////////////////////////////////////////////

void greeduce(Solution hint, size_t max_candidates_for_reduction, bool reduce_only_and_save) {
  int greeduce_start_time = elapsed_time_ms();

  int num_active_hyperedges = num_hyperedges;
  vector<bool> active(num_objects, true);
  vector<int> degree(num_objects);
  for (int i = 0; i < num_objects; ++i)
    degree[i] = incidence[i].size();
  queue<int> Q;
  for (int i = 0; i < num_objects; ++i)
    Q.push(i);
  vector<bool> enqueued(num_objects, true);

  vector<pair<int, int>> PQ;
  if (reduce_only_and_save) {
    PQ.emplace_back(-1, -1);
  } else {
    for (int level = 0; level < 2; ++level) {
      int start = PQ.size();
      for (int i = 0; i < num_vertices; ++i)
        if (hint[i] == level)
          PQ.emplace_back(1, i);
      shuffle(PQ.begin() + start, PQ.end(), prg);
      for (size_t i = start; i < PQ.size(); ++i) {
        if (PQ[i].first >= degree[PQ[i].second]) continue;
        pair<int, int> temp(PQ[i].first + 1, PQ[i].second);
        PQ.push_back(temp);
      }
    }
    reverse(PQ.begin(), PQ.end());
  }

  Solution solution(num_vertices, false);
  for (auto iter : PQ) {
    check_signal();
    while (!Q.empty()) {
      if (reduce_only_and_save && elapsed_time_ms() > REDUCTION_TIME_LIMIT_MS) break;
      check_signal();
      int x = Q.front();
      Q.pop();
      assert(enqueued[x]);
      enqueued[x] = false;
      if (max_candidates_for_reduction == 0) continue;
      if (!active[x]) continue;
      // unit edge rule
      if (x >= num_vertices && degree[x] == 1) {
        int u = -1;
        for (int v : incidence[x]) {
          if (active[v]) {
            assert(u == -1);
            u = v;
          }
        }
        assert(u != -1);
        active[u] = false;
        if (reduce_only_and_save) {
          always_in_solution.push_back((vertex_to_original_id[u]));
        } else {
          solution[u] = true;
        }
        for (int h : incidence[u]) {
          if (!active[h])
            continue;
          active[h] = false;
          --num_active_hyperedges;
          for (int v : incidence[h]) {
            --degree[v];
            if (active[v] && !enqueued[v]) {
              Q.push(v);
              enqueued[v] = true;
            }
          }
        }
        continue;
      }
      // domination rules
      set<int> candidates;
      for (int y : incidence[x]) {
        if (candidates.size() + incidence[y].size() > max_candidates_for_reduction)
          continue;
        candidates.insert(incidence[y].begin(), incidence[y].end());
      }
      for (int y : candidates) {
        assert((x < num_vertices) == (y < num_vertices));
        if (x==y || !active[y]) continue;
        bool x_is_subset_of_y = true;
        if (degree[x] > degree[y]) continue;
        for (int z : incidence[x]) {
          if (active[z] && !sorted_vector_contains(incidence[y], z)) {
            x_is_subset_of_y = false;
            break;
          }
        }
        if (x_is_subset_of_y) {
          int r = x < num_vertices ? x : y;
          active[r] = false;
          if (r >= num_vertices)
            --num_active_hyperedges;
          for (int z : incidence[r]) {
            --degree[z];
            if (active[z] && !enqueued[z]) {
              Q.push(z);
              enqueued[z] = true;
            }
          }
          if (r==x) break;
        }
      }
    }
    // at this point no more reductions available
    if (reduce_only_and_save) {
      prune_inactive_objects(active);
      return;
    }
    // make a greedy choice
    int v = iter.second;
    if (!active[v])
      continue;
    if (iter.first != degree[v])
      continue;
    active[v] = false;
    for (int h : incidence[v]) {
      if (!active[h]) continue;
      solution[v] = true;
      active[h] = false;
      --num_active_hyperedges;
      for (int u : incidence[h]) {
        --degree[u];
        if (active[u] && !enqueued[u]) {
          Q.push(u);
          enqueued[u] = true;
        }
      }
    }
    if (num_active_hyperedges == 0) break;
  }

  assert(num_active_hyperedges == 0);

  remove_unnecessary_nodes(solution);

  int sum = 0;
  for (bool x : solution)
    sum += x;
  if (sum < best_cost) {
    best = solution;
    best_cost = sum;
    if (DBG) cerr << "best " << sum << ' ' << elapsed_time_ms() << endl;
  }
  if (DBG) cerr << "greeduce " << max_candidates_for_reduction << " " << elapsed_time_ms() - greeduce_start_time << endl;

  check_signal();
}

void fast_lazy_reductions() {
  Solution dummy;
  greeduce(dummy, -1, true);
}

////////////////////////////////////////////////////////////////////////////////

void print_statictics(const string& header) {
  if (!DBG)
    return;
  vector<int> max_deg(num_hyperedges), deg(num_vertices);
  int total_deg = 0;
  for (int i = 0; i < num_vertices; ++i) {
    for (int h : incidence[i])
      max_deg[h - num_vertices] = max(max_deg[h - num_vertices], (int)incidence[i].size());
    deg[i] = incidence[i].size();
    total_deg += incidence[i].size();
  }
  double efficiency_lower_bound = 0.0;
  for (int d : max_deg)
    efficiency_lower_bound += 1.0/(double)d;
  efficiency_lower_bound = ceil(efficiency_lower_bound);
  nth_element(deg.begin(), deg.begin() + num_vertices * 99 / 100, deg.end());
  cerr << header
       << " N " << num_vertices
       << " H " << num_hyperedges
       << " AvD " << total_deg / num_vertices
       << " 99D " << deg[num_vertices * 99 / 100]
       << " LB " << efficiency_lower_bound
       << " T " << elapsed_time_ms() << endl;
}

////////////////////////////////////////////////////////////////////////////////

int main() {
  ios_base::sync_with_stdio(false);

  signal(SIGTERM, sigterm_handler);

  read_input();
  print_statictics("before reductions");
  fast_lazy_reductions();
  print_statictics("after reductions");

  best.assign(num_vertices, true);
  best_cost = num_vertices;

  uniform_int_distribution random_vertex(0, num_vertices - 1);

  for (int iter = 0;; ++iter) {
    Solution hint = best;
    if (iter > 0) {
      int mutation_size = best_cost < 5000 ? 50 : 15000;
      for (int k = 0; k < mutation_size; ++k)
        hint[random_vertex(prg)] = false;
    }
    greeduce(hint, min(iter, (int)1e6), false);
  }
}
