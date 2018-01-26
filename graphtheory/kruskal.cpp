/*
Kruskal's Algorithm using Disjoint-set forests.

Solves Kruskal (MST): Really Special Subtree from HackerRank.
*/

#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

// Structures
struct node {
    int rank;
    struct node* p;
};

// Functions
void MakeSet(struct node* x) {
  x->rank = 0;
  x->p = x;
}

struct node* FindSet(struct node* x) {
  // cout << "FINDSET" << endl;
  if (x->p != x) {
    // cout << "FINDSET IF" << endl;
    x->p = FindSet(x->p);
  }
  return x->p;
}

void Link(struct node* x, struct node* y) {

  if (x->rank > y->rank) {
    y->p = x;
  }
  else {
    x->p = y;
    if (x->rank == y->rank) {
      y->rank += 1;
    }
  }
  // cout << "LINK" << endl;
}

void Union (struct node* x, struct node* y) {
  Link(FindSet(x),FindSet(y));
}

int Kruskal(vector< vector<int> >& edges, int n) {
  vector<int> a (0);
  vector<struct node> gv(n);

  for (int i=0; i<n; i++) {
    MakeSet(&gv[i]);
  }

  for (int e=0; e<edges.size(); e++) {
    // cout << e << endl;
    // cout << &(gv[edges[e][2]].p) << endl;
    // cout << edges[e][2] << " " << edges[e][3] << endl;

    if (FindSet(&gv[edges[e][2]-1]) != FindSet(&gv[edges[e][3]-1])) {
      // cout << edges[e][2] << " " << edges[e][3] << endl;
      a.push_back(e);
      Union(&gv[edges[e][2]-1],&gv[edges[e][3]-1]);
    }
  }

  int tot = 0;
  for (auto e: a) {
    tot += edges[e][0];
  }

  return tot;
}

// Main
int main() {

  int n = 4, m = 5;
  vector< vector<int> > edges (m, vector<int>(4));

  edges[0][2] = 1; edges[0][3] = 2; edges[0][0] = 1; edges[0][1] = edges[0][2] + edges[0][3];
  edges[1][2] = 3; edges[1][3] = 2; edges[1][0] = 150; edges[1][1] = edges[1][2] + edges[1][3];
  edges[2][2] = 4; edges[2][3] = 3; edges[2][0] = 99; edges[2][1] = edges[2][2] + edges[2][3];
  edges[3][2] = 1; edges[3][3] = 4; edges[3][0] = 100; edges[3][1] = edges[3][2] + edges[3][3];
  edges[4][2] = 3; edges[4][3] = 1; edges[4][0] = 200; edges[4][1] = edges[4][2] + edges[4][3];

  sort(edges.begin(), edges.end());
  /*for (auto e: edges) {
    cout << e[2] << " " << e[3] << " " << e[0] << endl;
  }*/

  cout << Kruskal(edges,n) << endl;
  return 0;
}
