#include "gurobi_c++.h"
#include "KGraph.h"
#include "LazyConstraints.h"
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

string itos(int i);
string statusNumtoString(int num);
// random generator function:
int myrandom(int i);

long smallestFeasibleLatency2Robust(KGraph &g);

// overall latency-s CDS function, which calls others depending on values of (r,s). Assumes hop-based distances
vector<long> solveLatencyCDS(KGraph &g, long r, long s); 

// latency-s CDS funtions for specific values of (r,s) and hop-based distances
vector<long> solveCDSMIP(KGraph &g, bool &subOpt);   // r=1, s=2
vector<long> solveROBCDSMIP(KGraph &g, bool &subOpt);// r=2, s=2
vector<long> solveReCDSMIP(KGraph &g, bool &subOpt); // r=2, s=2
vector<long> solveMCDS(KGraph &g, long s);			 // r=1, s=arbitrary
vector<long> solveROBMCDS(KGraph &g, long s);		 // r=2, s=arbitrary
vector<long> solveReMCDS(KGraph &g, long s);		 // r=2, s=arbitrary
vector<long> solveMCDSrel(KGraph &g, long s); //r=1, s>2, LP Relaxation using just initial constraints

// node-weighted distance function for latency-s CDS
vector<long> solveMCDSweighted(KGraph &g, long s, vector<long> W);

// polynomial-size formulations for latency-s CDS when distances are hop-based
vector<long> solveVBCDSMIP(KGraph &g); // r=1, s=2. Same as solveCDSMIP above.
vector<long> solveVBMCDS(KGraph &g, long s); //r=1, s>2, VB version as-is (i.e., no extra domination constraints)
vector<long> solveVBMCDS2(KGraph &g, long s); //r=1, s>2, VB version 2 with domination constraints

// functions for verifying a latency-s CDS.
bool IsLatencyConstrainedCDS(KGraph &g, vector<bool> &D, long s);
bool IsLatencyConstrainedRCDS(KGraph &g, vector<bool> D, long s);
long EccentricitySubroutine(KGraph &g, vector<bool> &D, long v);
vector<long> ComputeSSSPinGBv(KGraph &g, vector<bool> &B, long v);

// functions used to minimalize length-s vertex cuts (to strengthen the formulation)
vector<long> Minimalize(KGraph &g, vector<bool> &B, long s);	// B is the infeasible CDS "solution"
vector<long> Minimalize(KGraph &g, vector<long> &CUT, long s);  // CUT are the vertices not in the CDS solution
vector< vector<long> > EnumerateFarPairs(KGraph &g, vector<bool> &B, long s);
vector< vector<long> > EnumerateFarPairsRob(KGraph &g, vector<bool> &B, long s);
vector<long> ComplementVector(vector<bool> &B, long n);

// slower versions of minimalize functions
vector<long> MinimalizeBasic(KGraph &g, vector<bool> &B, long s);
vector<long> MinimalizeBasic(KGraph &g, vector<long> &CUT, long s);
vector<long> MinimalizeNone(KGraph &g, vector<bool> &b, long s);
vector<long> MinimalizeNone(KGraph &g, vector<long> &CUT, long s);

// heuristics that we use
vector<bool> HeuristicLCDS(KGraph &g, long s);
vector<bool> HeuristicLCDS(KGraph &g, vector<bool> &heuristicSoln, long s);
vector<bool> HeuristicLCDSBestIn(KGraph &g, long s);

// unused heuristics (i.e., ones that we tried, but ended up not working very well compard to others)
vector<long> BiasedShuffle(KGraph &g);
vector<bool> RandomHeuristicLCDS(KGraph &g, long s, long iterations);
vector<bool> HeuristicLCDS2(KGraph &g, long s);
vector<bool> HeuristicLCDSminimalWeighted(KGraph &g, vector<bool> &heuristicSoln, long s, vector<long> W);//minimalizing heuristic in weighted variant 
vector<long> SortVerticesByIncreasingDegree(KGraph &g);

// heuristics that we use for r-robust latency-s CDS (case r=2)
vector<bool> HeuristicRLCDS(KGraph &g, vector<bool> &heuristicSoln,long s);
vector<bool> MinimalizeRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s);

// unused heuristics for r-robust latency-s CDS (r=2)
vector<bool> HeuristicRLCDS2(KGraph &g, long s);

// Tarjan's linear-time algorithm for finding biconnected components
vector< vector<long> > FindBiconnectedComponents(KGraph &g, vector<bool> &AV);
void Bico_Sub(long v, long u, long &i, KGraph &g, vector<long> &number, vector<long> &lowopt, stack<long> &le, stack<long> &re, vector< vector<long> > &BC);

// all of the same functions as above, but for the node-weighted variant
bool IsLatencyConstrainedCDSweighted(KGraph &g, vector<bool> &D, long s, vector<long> W);
bool IsLatencyConstrainedCDSweighted(KGraph &g, vector<long> &D, long s, vector<long> W);
vector< vector<long> > EnumerateFarPairsWeighted(KGraph &g, vector<bool> &B, long s, vector<long> W);
vector<long> MinimalizeWeighted(KGraph &g, vector<long> &CUT, long s, vector<long> W);
vector<long> MinimalizeWeighted(KGraph &g, vector<bool> &B, long s, vector<long> W);
vector<long> ComputeSSSPinGBvWeighted(KGraph &g, vector<bool> &B, long origin, vector<long> W);

// unused weighted functions
vector<bool> HeuristicLCDSweighted(KGraph &g, long s, vector<long> W);
vector<bool> HeuristicLCDSweighted2(KGraph &g, vector<bool> &heuristicSoln, long s, vector<long> W);
vector<bool> HeuristicLCDSweighted2(KGraph &g, long s, vector<long> W);
