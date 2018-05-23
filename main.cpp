#include "GRBInterface.h"
#include "KGraph.h"
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;

int main(int argc, char *argv[])
{
	time_t start_time = clock();
	if (argc < 2)
	{
		cout << "ERROR: Not enough arguments.";
	}
	else if (strcmp(argv[1], "CDS") == 0) // Solves minimum CDS problem with no latency constraints (i.e., latency-s CDS with s=n-1).
	{
		if (argc < 4)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[3], argv[3], argv[2]);
		long s = g.n - 1;
		vector<long> D;
		long r = 1;
		time_t start = clock();
		D = solveMCDS(g, s);
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "LCDSK") == 0) // Solving latency-s CDS (where s = diameter + k ). Only for r=1. 
	{
		if (argc < 5)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[4], argv[4], argv[3]);
		long k = atol(argv[2]);
		bool subOpt;
		long s = g.DiameterUnweighted() + k;
		vector<long> D;
		long r = 1;
		time_t start = clock();
		if (s > 2)
			D = solveMCDS(g, s);
		else if (s == 2)
			D = solveCDSMIP(g, subOpt);
		else if (s < 2)
			cout << "No CDS nodes exists!" << endl;
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "LCDSKRel") == 0) // Solving the LP relaxation (initial constraints only) of latency-s CDS
	{
		if (argc < 5)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[4], argv[4], argv[3]);
		long k = atol(argv[2]);
		bool subOpt;
		long s = g.DiameterUnweighted() + k;
		vector<long> D;
		long r = 1;
		time_t start = clock();
		if (s > 2)
			D = solveMCDSrel(g, s);
		else if (s == 2)
			D = solveCDSMIP(g, subOpt);
		else if (s < 2)
			cout << "No CDS nodes exists!" << endl;
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "VBLCDSK") == 0) // Solving latency-s CDS problem based on VB formulation as-is (i.e., no domination constraints)
	{
		if (argc < 5)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[4], argv[4], argv[3]);
		long k = atol(argv[2]);
		long s = g.DiameterUnweighted() + k;
		vector<long> D;
		long r = 1;
		time_t start = clock();
		if (s > 2)
			D = solveVBMCDS(g, s);
		else if (s == 2)
			D = solveVBCDSMIP(g);
		else if (s < 2)
			cout << "No CDS nodes exists!" << endl;
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "VBLCDSK2") == 0) // Solving latency-s CDS problem based on VB formulation with extra domination constraints
	{
		if (argc < 5)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[4], argv[4], argv[3]);
		long k = atol(argv[2]);
		long s = g.DiameterUnweighted() + k;
		vector<long> D;
		time_t start = clock();
		if (s > 2)
			D = solveVBMCDS2(g, s);
		else if (s == 2)
			D = solveVBCDSMIP(g);
		else if (s < 2)
			cout << "No CDS nodes exists!" << endl;
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "RLCDS") == 0) //r-robust latency-s CDS (r=2 only)
	{
		if (argc < 4)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[3], argv[3], argv[2]);
		long s = smallestFeasibleLatency2Robust(g);
		vector<long> D;
		bool subOpt;
		time_t start = clock();
		if (s > 2)
			D = solveROBMCDS(g, s);
		else if (s == 2)
			D = solveROBCDSMIP(g, subOpt);
		else if (s < 2)
			cout << "NO cds exists!" << endl;
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "LCDSW") == 0)// latency-s CDS with node-weighted distances
	{
		if (argc < 4)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g(argv[3], argv[3], argv[2]);
		
		vector<long> W(g.n, (long)0);
		vector<long> distv;
		for (long i = 0; i < g.n; i++)
		{
			double w = 0;
			distv = g.ShortestPathsUnweighted(i);
			for (long j = 0; j < distv.size(); j++)
			{
				w += distv[j];
			}
			w = (double)1000 * (g.n - 1) / w;
			W[i] = (long)floor(w);
			//W[i] = 1;   /// use to debug. should give same answers as unweighted case
			cerr << "W [ " << i << " ] = " << W[i] << endl;
		}
		
		long s = g.DiameterWeighted(W);
		time_t start = clock();
		vector<long> D = solveMCDSweighted(g, s, W);
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "diameter") == 0) // Computing (hop-based) diameter of a graph
	{
		KGraph g(argv[3], argv[3], argv[2]);
		time_t start = clock();
		cerr << "diameter =  " << g.DiameterUnweighted() << endl;
		cerr << "time in secs = " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}
	else if (strcmp(argv[1], "graphviz") == 0) //Graphviz codes
	{
		if (argc<5)
		{
			cout << argv[0] << " graphviz <graph_type> <graph_file> <graphviz_file>\n<graph_type>\tdimacs, snap_d\n";
			return 0;
		}
		KGraph G("temp", argv[3], argv[2]);
		G.WriteGVizGraph(argv[4]);
		cout << "Time = " << (double)(clock() - start_time) / CLOCKS_PER_SEC << " secs." << endl;
		return 0;
	}
	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
