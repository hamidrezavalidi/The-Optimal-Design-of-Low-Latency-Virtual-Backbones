#include "GRBInterface.h"
#include "LazyConstraints.h"
#include "KGraph.h"
#include <sstream>
#include <vector>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }

string statusNumtoString(int num)
{
	if (num == 1) return "Model is loaded, but no solution information is available.";
	else if (num == 2) return "Model was solved to optimality (subject to tolerances), and an optimal solution is available.";
	else if (num == 3) return "Model was proven to be infeasible.";
	else if (num == 4) return "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.";
	else if (num == 5) return "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.";
	else if (num == 6) return "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.";
	else if (num == 7) return "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.";
	else if (num == 8) return "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.";
	else if (num == 9) return "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.";
	else if (num == 10) return "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.";
	else if (num == 11) return "Optimization was terminated by the user.";
	else if (num == 12) return "Optimization was terminated due to unrecoverable numerical difficulties.";
	else if (num == 13) return "Unable to satisfy optimality tolerances; a sub-optimal solution is available.";
	else if (num == 14) return "An asynchronous optimization call was made, but the associated optimization run is not yet complete.";
	else if (num >= 15 || num <= 0) return "No specific error could be recognized.";
}

long smallestFeasibleLatency2Robust(KGraph &g)
{
	long s = g.DiameterUnweighted();
	vector<bool> allNodes(g.n, true);
	vector<long> dist;

	for (long i = 0; i < g.n; i++)
	{
		allNodes[i] = false;

		// find the smallest s such that V-v is a latency-s CDS
		for (long k = 0; k < g.n; k++)
		{
			dist = ComputeSSSPinGBv(g, allNodes, k);
			for (long j = 0; j < g.n; j++) if(dist[j]>s) s = dist[j];
		}
		allNodes[i] = true;
	}
	return s;
}

vector<long> solveLatencyCDS(KGraph &g, long r, long s)
{
	bool subOpt;
	vector<long> D;
	if (r == 1 && s == 2) D = solveCDSMIP(g, subOpt);
	else if (r == 2 && s == 2) D = solveROBCDSMIP(g, subOpt);
	else if (r == 1 && s > 2) D = solveMCDS(g, s);
	else if (r == 2 && s > 2) D = solveROBMCDS(g, s);
	else cout << "ERROR: Your values of (r,s) are NOT supported." << endl;

	cerr << "Nodes in CDS are : ";
	for (long i = 0; i < D.size(); i++) cerr << D[i] << " ";
	cerr << endl;

	return D;
}

vector<long> solveROBCDSMIP(KGraph &g, bool &subOpt)
{
	long s = 2;
	vector<long> D;
	// check if the graph has any articulation vetices
	vector<bool> ArticulationVertices; // (g.n, false);
									   //vector< vector<long>> Z = FindBiconnectedComponents (g, ArticulationVertices);
	FindBiconnectedComponents(g, ArticulationVertices);
	long NumArticulationVertices = count(ArticulationVertices.begin(), ArticulationVertices.end(), true);
	if (NumArticulationVertices > 0)
	{
		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " ";
		cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		return D;
	}
	else
	{
		
		//subOpt = true;

		try
		{
			cerr << "Creating Gurobi model variable." << endl;
			GRBEnv env = GRBEnv();
			env.set(GRB_IntParam_OutputFlag, 0);
			env.set(GRB_DoubleParam_TimeLimit, 3600);
			GRBModel model = GRBModel(env);
			GRBVar *X = model.addVars(g.n, GRB_BINARY);
			model.update();

			cerr << "Adding objective function.\n";
			for (long w = 0; w < g.n; w++)
				X[w].set(GRB_DoubleAttr_Obj, 1);
			model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
			model.update();
			for (long i = 0; i < g.n; i++)
			{
				//Adding length-s vertex cut constraints
				vector<bool> neighbors(g.n, false);
				long v;
				for (long j = 0; j < g.degree[i]; j++)
				{
					v = g.adj[i][j];
					neighbors[v] = true;
				}

				for (long j = i + 1; j < g.n; j++)
				{
					if (neighbors[j]) continue;

					vector<long> p = g.CommonNeighborsList(i, j);

					GRBLinExpr expr = 0;
					for (long k = 0; k < p.size(); k++)
						expr += X[p[k]];

					string s = "Con_" + itos(i) + itos(j);
					model.addConstr(expr >= 2, s);
					//cerr << "adding constraint.\n";

				}
			}

			cerr << "Model is solving now." << endl;
			model.update();

			vector <bool> initialSolution(g.n, false);
			double TotaltimeinHeuristic = 0;
			time_t start1 = clock();

			vector<bool> soln2 = HeuristicLCDSBestIn(g, s);
			cerr << "Best-in heuristic soln = " << count(soln2.begin(), soln2.end(), true) << endl;
			soln2 = HeuristicRLCDS(g, soln2, s);
			cerr << "Best-in heuristic soln (for r = 2, after augmenting the solution) = " << count(soln2.begin(), soln2.end(), true) << endl;
			soln2 = MinimalizeRLCDS(g, soln2, s);
			cerr << "Best-in heuristic soln (for r=2, after minimalizaing the solution) = " << count(soln2.begin(), soln2.end(), true) << endl;
			initialSolution = soln2;
			
			TotaltimeinHeuristic = (double)(clock() - start1) / CLOCKS_PER_SEC;
			for (long i = 0; i < g.n; i++)
			{
				if (initialSolution[i])
				{
					X[i].set(GRB_DoubleAttr_Start, 1.0);
				}
				else
				{
					X[i].set(GRB_DoubleAttr_Start, 0.0);
				}
			}
			model.update();
			model.optimize();

			cerr << "Extracting solution" << endl;
			int status = model.get(GRB_IntAttr_Status);
			cerr << statusNumtoString(status) << endl;

			for (long l = 0; l < g.n; l++)
				if (X[l].get(GRB_DoubleAttr_X) > 0.5)
					D.push_back(l);

			long NumNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
			double lb = model.get(GRB_DoubleAttr_ObjBound);
			cerr << "MIP LB = " << lb << endl;

			cerr << "Number of lazy cuts : " << " 0" << endl;
			cerr << "Number of branch-and-bound nodes :  " << NumNodes << endl;


			cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << count(initialSolution.begin(), initialSolution.end(), true) << " " << D.size() << " " << lb << " " << NumNodes << " " << "N/A" << " " << TotaltimeinHeuristic << " " << "N/A" << " ";

		}
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during optimization" << endl;
		}

		return D;

	}

}

vector<long> solveReCDSMIP(KGraph &g, bool &subOpt)
{
	long s = 2;
	vector<long> D;
	// check if the graph has any articulation vetices
	vector<bool> ArticulationVertices; // (g.n, false);
									   //vector< vector<long>> Z = FindBiconnectedComponents (g, ArticulationVertices);
	FindBiconnectedComponents(g, ArticulationVertices);
	long NumArticulationVertices = count(ArticulationVertices.begin(), ArticulationVertices.end(), true);
	if (NumArticulationVertices > 0)
	{
		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " ";
		cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		return D;
	}
	else
	{

		try
		{
			cerr << "Creating Gurobi model variable." << endl;
			GRBEnv env = GRBEnv();
			env.set(GRB_IntParam_OutputFlag, 0);
			env.set(GRB_DoubleParam_TimeLimit, 3600);
			GRBModel model = GRBModel(env);
			GRBVar *X = model.addVars(g.n, GRB_BINARY);
			model.update();

			cerr << "Adding objective function.\n";
			for (long w = 0; w < g.n; w++)
				X[w].set(GRB_DoubleAttr_Obj, 1);
			model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
			model.update();
			vector<long> minimalVertexCut;
			long v;
			for (long i = 0; i < g.n; i++)
			{
				//Adding vertex cut constraints
				minimalVertexCut = Minimalize(g, g.adj[i], g.n - 1);
				GRBLinExpr expr1 = 0;
				for (long j = 0; j < minimalVertexCut.size(); j++)
				{
					v = minimalVertexCut[j];
					expr1 += X[v];
				}
				model.addConstr(expr1 >= 2);
				//Adding length-s vertex cut constraints
				vector<bool> neighbors(g.n, false);
				long v;
				for (long j = 0; j < g.degree[i]; j++)
				{
					v = g.adj[i][j];
					neighbors[v] = true;
				}

				for (long j = i + 1; j < g.n; j++)
				{
					if (neighbors[j]) continue;

					vector<long> p = g.CommonNeighborsList(i, j);

					GRBLinExpr expr = 0;
					for (long k = 0; k < p.size(); k++)
						expr += X[p[k]];

					string s = "Con_" + itos(i) + itos(j);
					model.addConstr(expr >= 1, s);
					//cerr << "adding constraint.\n";

				}
			}

			cerr << "Model is solving now." << endl;
			model.update();

			cerr << "Adding lazy constraints";

			LazyConstraints5 cb = LazyConstraints5(X, g, s);	// tell Gurobi which function generates the lazy cuts.
			model.setCallback(&cb);

			model.optimize();

			cerr << "Extracting solution" << endl;
			int status = model.get(GRB_IntAttr_Status);
			cerr << statusNumtoString(status) << endl;

			for (long l = 0; l < g.n; l++)
				if (X[l].get(GRB_DoubleAttr_X) > 0.5)
					D.push_back(l);

			long NumNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
			double gap = model.get(GRB_DoubleAttr_MIPGap);
			cerr << "MIP gap = " << gap << endl;

			cerr << "Number of lazy cuts : " << " 0" << endl;
			cerr << "Number of branch-and-bound nodes :  " << NumNodes << endl;


			cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << "NA" << " " << D.size() << " " << gap << " " << NumNodes << " " << LazyConstraints5::numLazyCuts + LazyConstraints5::numLazyCuts2 << " " << "NA" << " " << LazyConstraints5::TotalCallbackTime + LazyConstraints5::TotalCallbackTime2 << " ";

		}
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during optimization" << endl;
		}

		return D;

	}

}


vector<long> solveCDSMIP(KGraph &g, bool &subOpt)
{
	vector<long> D;
	try
	{
		cerr << "Creating Gurobi model variable." << endl;
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		model.update();

		cerr << "Adding objective function.\n";
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();
		for (long i = 0; i < g.n; i++)
		{
			vector<bool> neighbors(g.n, false);
			long v;
			for (long j = 0; j < g.degree[i]; j++)
			{
				v = g.adj[i][j];
				neighbors[v] = true;
			}

			for (long j = i + 1; j < g.n; j++)
			{
				if (neighbors[j]) continue;

				vector<long> p = g.CommonNeighborsList(i, j);
				GRBLinExpr expr = 0;
				for (long k = 0; k < p.size(); k++)
					expr += X[p[k]];
				string s = "Con_" + itos(i) + itos(j);
				model.addConstr(expr >= 1, s);
			}
		}

		cerr << "Model is solving now." << endl;
		model.update();
	
		vector<bool> ArticulationVertices(g.n, false);
		vector< vector<long>> Z = FindBiconnectedComponents(g, ArticulationVertices);

		for (int i = 0; i < g.n; i++)
		{
			if (ArticulationVertices[i])
			{
				X[i].set(GRB_DoubleAttr_LB, 1.0);
			}
		}


		model.update();

		model.optimize();

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;

		for (long l = 0; l<g.n; l++)
			if (X[l].get(GRB_DoubleAttr_X)>0.5)
				D.push_back(l);

		long NumNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;
		
		cerr << "Number of lazy cuts : " << " 0" << endl;
		cerr << "Number of branch-and-bound nodes :  " << NumNodes << endl;

		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << 2 << " " << "NA" << " " << D.size() << " " << lb << " " << NumNodes << " " << "N/A" << " " << "NA" << " " << "N/A" << " ";
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}

	return D;

}

vector<long> solveVBCDSMIP(KGraph &g)
{
	bool subopt;
	return solveCDSMIP(g,subopt);
}

vector<long> SortVerticesByIncreasingDegree(KGraph &g) //Sorting vertices by increasing degree
{
	vector<long> emptyBucket;
	vector< vector<long> > degreeBuckets(g.n,emptyBucket);

	long deg;
	for (long i = 0; i < g.n; i++)
	{
		deg = g.degree[i];
		degreeBuckets[deg].push_back(i);
	}

	vector<long> sortedList(g.n,-1);
	long count = 0;
	long v;
	for (long i = 0; i < degreeBuckets.size(); i++)
	{
		for (long j = 0; j < degreeBuckets[i].size(); j++)
		{
			v = degreeBuckets[i][j];
			sortedList[count] = v;
			count++;
		}
	}

	return sortedList;
}

bool IsLatencyConstrainedCDS(KGraph &g, vector<bool> &D, long s)
{
	long ev;
	for (long v = 0; v < g.n; v++)
	{
		ev = EccentricitySubroutine(g, D, v);
		if (ev > s) return false;
	}
	return true;
}
bool IsLatencyConstrainedCDSweighted(KGraph &g, vector<long> &D, long s, vector<long> W)
{
	long v;
	vector<bool> Dbool(g.n, false);
	for (long i = 0; i < D.size(); i++)
	{
		v = D[i];
		Dbool[v] = true;
	}
	return IsLatencyConstrainedCDSweighted(g, Dbool, s, W);
}
bool IsLatencyConstrainedCDSweighted(KGraph &g, vector<bool> &D, long s, vector<long> W)
{
	long ev;
	vector<long> SP;
	for (long v = 0; v < g.n; v++)
	{
		SP = ComputeSSSPinGBvWeighted(g, D, v, W);
		ev = *max_element(SP.begin(), SP.end());
		if (ev > s) return false;
	}
	return true;
}

vector<long> ComplementVector(vector<bool> &B, long n)
{
	vector<long> Q;					// complement of B (i.e., V\B), where V = {0, 1, ..., n-1 }
	for (long i = 0; i < n; i++)
	{
		if (!B[i]) Q.push_back(i);
	}
	return Q;
}
vector<long> Minimalize(KGraph &g, vector<long> &CUT, long s)
{
	vector<bool> B(g.n, true);
	long v;
	for (long i = 0; i < CUT.size(); i++)
	{
		v = CUT[i];
		B[v] = false;
	}
	return Minimalize(g, B, s);
}
vector<long> MinimalizeBasic(KGraph &g, vector<long> &CUT, long s)
{
	vector<bool> B(g.n, true);
	long v;
	for (long i = 0; i < CUT.size(); i++)
	{
		v = CUT[i];
		B[v] = false;
	}
	return MinimalizeBasic(g, B, s);
}
vector< vector<long> > EnumerateFarPairs(KGraph &g, vector<bool> &B, long s)
{
	vector< vector<long> > FarPairs;
	vector<long> dist;
	vector<long> Pair;

	for (long i = 0; i < g.n; i++)
	{
		dist = ComputeSSSPinGBv(g, B, i);
		for (long j = i + 1; j < g.n; j++)
		{
			if (dist[j] > s)
			{
				Pair.push_back(i);
				Pair.push_back(j);
				FarPairs.push_back(Pair);
				Pair.clear();
			}
		}
	}

	return FarPairs;
}

vector< vector<long> > EnumerateFarPairsRob(KGraph &g, vector<bool> &B, long s)
{
	vector< vector<long> > FarPairs;
	vector<long> dist;
	vector<long> Pair;

	for (long i = 0; i < g.n; i++)
	{
			dist = ComputeSSSPinGBv(g, B, i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist[j] > s)
				{
					Pair.push_back(i);
					Pair.push_back(j);
					FarPairs.push_back(Pair);
					Pair.clear();
				}
				else
				{
					for (long k = 0; k < B.size(); k++) 
					{
						if (!B[k]) continue;
						B[k] = false;
						dist = ComputeSSSPinGBv(g, B, i);
						if (dist[j] > s)
						{
							Pair.push_back(i);
							Pair.push_back(j);
							FarPairs.push_back(Pair);
							Pair.clear();
							break;
						}
					}
				}
			}
		
	}

	return FarPairs;
}

vector< vector<long> > EnumerateFarPairsWeighted(KGraph &g, vector<bool> &B, long s, vector<long> W)
{
	vector< vector<long> > FarPairs;
	vector<long> dist;
	vector<long> Pair;

	for (long i = 0; i < g.n; i++)
	{
		dist = ComputeSSSPinGBvWeighted(g, B, i, W);
		for (long j = 0; j < g.n; j++)
		{
			{
				if (dist[j] > s)
				{
					Pair.push_back(i);
					Pair.push_back(j);
					FarPairs.push_back(Pair);
					Pair.clear();
				}
			}
		}
	}

	return FarPairs;
}
vector<long> MinimalizeWeighted(KGraph &g, vector<long> &CUT, long s, vector<long> W)
{
	vector<bool> B(g.n, true);  // these are the CDS nodes (i.e., the complement of the cut
	long v;
	for (long i = 0; i < CUT.size(); i++)
	{
		v = CUT[i];
		B[v] = false;
	}
	return MinimalizeWeighted(g, B, s, W);
}

vector<long> MinimalizeWeighted(KGraph &g, vector<bool> &B, long s, vector<long> W)
{
	vector<long> Q = ComplementVector(B, g.n);					// complement of initial B (i.e., V\B)
	long q;
	long a, b;

	vector< vector<long> > FarPairs = EnumerateFarPairsWeighted(g, B, s, W);
	vector<long> dist;
	for (long i = 0; i < Q.size(); i++)
	{
		q = Q[i];
		dist = ComputeSSSPinGBvWeighted(g, B, q, W);
		bool StillTooFar = false;
		for (long j = 0; j < FarPairs.size() && !StillTooFar; j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];

			if (dist[a] + dist[b] - W[q] + W[a] > s)
			{
				StillTooFar = true;
			}
		}

		if (StillTooFar)
		{
			// move q to set B
			B[q] = true;

			// create new set of FarPairs
			vector< vector<long> > newFarPairs;
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];

				if (dist[a] + dist[b]- W[q] + W[a] > s)
				{
					newFarPairs.push_back(FarPairs[j]);
				}
			}
			FarPairs = newFarPairs;
		}
	}
	return ComplementVector(B, g.n);
}

vector<long> Minimalize(KGraph &g, vector<bool> &B, long s)
{
	vector<long> Q = ComplementVector(B, g.n);					// complement of initial B (i.e., V\B)
	long q;
	long a, b;

	vector< vector<long> > FarPairs = EnumerateFarPairs(g, B, s);
	vector<long> dist;

	// make B maximal but not a latency-s CDS
	for (long i = 0; i < Q.size(); i++)
	{
		q = Q[i];

		// compute SSSP in graph G^B_q
		dist = ComputeSSSPinGBv(g, B, q);

		// check if FarPairs will become empty if we add q
		bool StillTooFar = false;
		for (long j = 0; j < FarPairs.size() && !StillTooFar; j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];

			if (dist[a] + dist[b] > s)
			{
				StillTooFar = true;
			}
		}

		if (StillTooFar)
		{
			// move q to set B
			B[q] = true; 

			// create new set of FarPairs
			vector< vector<long> > newFarPairs;
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];

				if (dist[a] + dist[b] > s)
				{
					newFarPairs.push_back(FarPairs[j]);
				}
			}
			FarPairs = newFarPairs;
		}
	}

	// return complement of *final* B
	return ComplementVector(B, g.n);
}
vector<long> MinimalizeNone(KGraph &g, vector<bool> &B, long s)
{
	vector<long> Q = ComplementVector(B, g.n);
	return Q;
}

vector<long> MinimalizeNone(KGraph &g, vector<long> &CUT, long s)
{
	return CUT;
}

vector<long> MinimalizeBasic(KGraph &g, vector<bool> &B, long s)
{
	vector<long> Q = ComplementVector(B,g.n);		// complement of initial B (i.e., V\B)

	long q;
	// make B maximal but not a latency-s CDS
	for (long i = 0; i < Q.size(); i++)
	{
		q = Q[i];
		B[q] = true;
		if (IsLatencyConstrainedCDS(g, B, s))
		{
			B[q] = false;
		}
	}
	
	// return complement of *final* B
	return ComplementVector(B, g.n);
}

long EccentricitySubroutine(KGraph &g, vector<bool> &D, long v)
{
	vector<long> dist = ComputeSSSPinGBv(g, D, v);
	long ev = 0;
	for (long i = 0; i < g.n; i++)
	{
		if (dist[i] > ev) ev = dist[i];
	}
	return ev;

}
vector<long> ComputeSSSPinGBv(KGraph &g, vector<bool> &B, long v) //Computing single source shortest path
{
	vector<long> dist(g.n, g.n);
	vector<long> parents;
	vector<long> children;

	dist[v] = 0;
	long u, w;

	children.push_back(v);

	while (!children.empty())
	{
		parents = children;
		children.clear();
		for (long i = 0; i < parents.size(); i++)
		{
			w = parents[i];
			for (long j = 0; j < g.degree[w]; j++) {
				u = g.adj[w][j];
				if (dist[u] == g.n) {
					dist[u] = dist[w] + 1;
					if (B[u])
						children.push_back(u);
				}
			}
		}
	}
	return dist;
}

vector<long> ComputeSSSPinGBvWeighted(KGraph &g, vector<bool> &B, long origin, vector<long> W) //Computing weighted single source shortest path
{
	long infty = long(floor((double)0.4*LONG_MAX));
	vector<long> dist(g.n,infty); //shortest distance from origin node to each other node. dist[i] = \infty means i not reachable
	vector<bool> Q(g.n, true); // the set of vertices whose distance labels are not yet permanent
	//W[origin] = 0;   // for convenience, since we are not counting the weight of the first vertex in the path as part of its length
	dist[origin] = 0; //the origin node is distance 0 from itself
	long minDistance;
	long minVertex;

	// do node-weighted version of Dijkstra's shortest path algorithm.
	for (long i = 0; i < g.n; i++)
	{
		// find a vertex u from Q of minimum distance label dist[u]
		minDistance = LONG_MAX;
		minVertex = -1;
		for (long v = 0; v < g.n; v++)
		{
			if (!Q[v]) continue;
			if (dist[v] < minDistance)
			{
				minDistance = dist[v];
				minVertex = v;
			}
		}

		// remove minVertex from Q
		Q[minVertex] = false;

		if (B[minVertex] || minVertex == origin)  // only allow communication paths to emanate from the origin or from nodes in the CDS B
		{
			// update distance labels of neighbors
			for (long j = 0; j < g.degree[minVertex]; j++)
			{
				long v = g.adj[minVertex][j];
				if (Q[v] && dist[minVertex] + W[minVertex] < dist[v])
				{
					dist[v] = dist[minVertex] + W[minVertex];
				}
			}
		}
	}
	return dist;
}

vector<bool> HeuristicLCDS(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, true);
	return HeuristicLCDS(g, heuristicSoln, s);
}
vector<bool> HeuristicLCDS(KGraph &g, vector<bool> &heuristicSoln, long s)
{
	vector<long> sortedList = SortVerticesByIncreasingDegree(g);
	long v;

	for (int i = 0; i < g.n; i++) {
		v = sortedList[i];
		if (!heuristicSoln[v]) continue;
		heuristicSoln[v] = false;
		if (!IsLatencyConstrainedCDS(g, heuristicSoln, s))
		{
			heuristicSoln[v] = true;
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<bool> HeuristicLCDSminimalWeighted(KGraph &g, vector<bool> &heuristicSoln, long s, vector<long> W)
{
	vector<long> sortedList = SortVerticesByIncreasingDegree(g);
	long v;

	for (int i = 0; i < g.n; i++) {
		v = sortedList[i];
		if (!heuristicSoln[v]) continue;
		heuristicSoln[v] = false;
		if (!IsLatencyConstrainedCDSweighted(g, heuristicSoln, s, W))
		{
			heuristicSoln[v] = true;
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<bool> HeuristicLCDSweighted(KGraph &g, long s, vector<long> W)
{
	vector<bool> heuristicSoln(g.n, false);
	long wmax = -1;
	for (long i = 0; i < W.size(); i++) 
	{
		if (W[i] > wmax)
			wmax = W[i];
	}
	vector< vector<long> > FarPairs = EnumerateFarPairsWeighted(g, heuristicSoln, wmax, W);
	vector<long> dist;
	vector<long> score;
	vector<long> scoreZero(g.n, 0);

	long a, b;
	long max;
	long imax;
	while (!FarPairs.empty())
	{
		cerr << FarPairs.size() << " ";
		score = scoreZero;
		max = -1;
		imax = -1;
		for (long i = 0; i < g.n; i++)
		{
			if (heuristicSoln[i]) continue;
			dist = ComputeSSSPinGBvWeighted(g, heuristicSoln, i, W);
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] - W[i] + W[a] <= s)
				{
					score[i]++;
				}

			}

			if (score[i] > max) {
				max = score[i];
				imax = i;
			}
		}
		heuristicSoln[imax] = true;

		vector< vector<long> > newFarPairs;
		vector<long> newPair;
		dist = ComputeSSSPinGBvWeighted(g, heuristicSoln, imax, W);
		for (long j = 0; j < FarPairs.size(); j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];
			if (dist[a] + dist[b] - W[imax] + W[a] > s)
			{
				newPair.clear();
				newPair.push_back(a);
				newPair.push_back(b);
				newFarPairs.push_back(newPair);
			}
		}
		FarPairs = newFarPairs;
		if (FarPairs.size() == 1)
		{
			a = FarPairs[0][0];
			b = FarPairs[0][1];
			dist = ComputeSSSPinGBvWeighted(g, heuristicSoln, a, W);
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<bool> HeuristicLCDSweighted2(KGraph &g, long s, vector<long> W)
{
	vector<bool> heuristicSoln(g.n, true);
	return HeuristicLCDSweighted2(g, heuristicSoln, s, W);
}

vector<bool> HeuristicLCDSweighted2(KGraph &g, vector<bool> &heuristicSoln, long s, vector<long> W)
{
	vector<long> sortedList = SortVerticesByIncreasingDegree(g);
	long v;

	for (int i = 0; i < g.n; i++) {
		v = sortedList[i];
		if (!heuristicSoln[v]) continue;
		heuristicSoln[v] = false;
		if (!IsLatencyConstrainedCDSweighted(g, heuristicSoln, s, W))
		{
			heuristicSoln[v] = true;
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<bool> HeuristicLCDS2(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, true);
	vector<long> sortedList = SortVerticesByIncreasingDegree(g);
	long v;

	for (int i = 0; i < g.n; i++) 
	{
		v = sortedList[i];
		heuristicSoln[v] = false;
		if (!IsLatencyConstrainedCDS(g, heuristicSoln, s))
		{
			heuristicSoln[v] = true;
		}
		else
		{
			for (long j = 0; j < g.n && !heuristicSoln[v]; j++)
			{
				if (j != v && heuristicSoln[j])
				{
					heuristicSoln[j] = false;
					if (!IsLatencyConstrainedCDS(g, heuristicSoln, s)) heuristicSoln[v] = true;
					heuristicSoln[j] = true;
				}
			}
		}
	}
	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

bool IsLatencyConstrainedRCDS(KGraph &g, vector<bool> D, long s)
{
	for (long i = 0; i < g.n; i++)
	{
		if (!D[i]) continue;
		D[i] = false;
		if (!IsLatencyConstrainedCDS(g, D, s))
		{
			return false;
		}
		D[i] = true;
	}
	return true;
}

vector<bool> MinimalizeRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s)
{
	for (long i = 0; i < g.n; i++)
	{
		if (!heuristicSoln[i]) continue;
		// now node i belongs to the heuristic solution for r=2. Can we remove it and still be feasibe?
		heuristicSoln[i] = false;
		if (!IsLatencyConstrainedRCDS(g, heuristicSoln, s)) heuristicSoln[i] = true;
		// else we didn't need node i in the solution, and we can keep it set to false.
	}
	return heuristicSoln;
}

vector<bool> HeuristicRLCDS(KGraph &g, vector<bool> &heuristicSoln, long s)
{
	vector<bool> initialHeuristicSoln = heuristicSoln;
	vector< vector<long> > FarPairs;
	vector< vector<long> > newFarPairs;
	vector<long> dist;
	vector<long> newPair;
	vector<long> score;
	long a, b;
	long max;
	long imax;
	vector<long> scoreZero(g.n, 0);

	for (long v = 0; v < g.n; v++)
	{
		if (!initialHeuristicSoln[v]) continue;  //if node was not in initial heuristic soln, then removing it will not make current heuristic soln infeasibe
		heuristicSoln[v] = false;
		FarPairs = EnumerateFarPairs(g, heuristicSoln, s);
		//cerr << "v = " << v << endl;
		while (!FarPairs.empty())
		{
			score = scoreZero;
			score[v] = -1;
			max = -1;
			imax = -1;
			for (long i = 0; i < g.n; i++)
			{
				if (heuristicSoln[i] || i==v) continue;
				dist = ComputeSSSPinGBv(g, heuristicSoln, i);
				for (long j = 0; j < FarPairs.size(); j++)
				{
					a = FarPairs[j][0];
					b = FarPairs[j][1];
					if (dist[a] + dist[b] <= s)
					{
						score[i]++;
					}

				}

				if (score[i] > max) {
					max = score[i];
					imax = i;
				}
			}
			heuristicSoln[imax] = true;
			newFarPairs.clear();
			
			dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] > s)
				{
					newPair.clear();
					newPair.push_back(a);
					newPair.push_back(b);
					newFarPairs.push_back(newPair);
				}
			}

			FarPairs = newFarPairs;

		}
		heuristicSoln[v] = true;
	}

	return heuristicSoln;
}


vector<bool> HeuristicRLCDS2(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, false);
	vector< vector<long> > FarPairs = EnumerateFarPairs(g, heuristicSoln, 1);
	vector<long> dist;
	vector<long> score;
	vector<long> scoreZero(g.n, 0);
	vector< vector<long> > newFarPairs;
	vector<long> newPair;

	long a, b;
	long max;
	long imax;
	while (!FarPairs.empty())
	{
		score = scoreZero;
		max = -1;
		imax = -1;
		for (long i = 0; i < g.n; i++)
		{
			if (heuristicSoln[i]) continue;
			dist = ComputeSSSPinGBv(g, heuristicSoln, i);
			for (long j = 0; j < FarPairs.size(); j++)
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] <= s)
				{
					score[i]++;
				}

			}

			if (score[i] > max) {
				max = score[i];
				imax = i;
			}
		}
		heuristicSoln[imax] = true;

		dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
		for (long j = 0; j < FarPairs.size(); j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];
			if (dist[a] + dist[b] > s)
			{
				newPair.clear();
				newPair.push_back(a);
				newPair.push_back(b);
				newFarPairs.push_back(newPair);
			}
			else
			{
				for (long k = 0; k < heuristicSoln.size(); k++)
				{
					if (k == imax) continue;
					if (!heuristicSoln[k]) continue;
					heuristicSoln[k] = false;
					dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
					heuristicSoln[k] = true;
					if (dist[a] + dist[b] > s)
					{
						newPair.clear();
						newPair.push_back(a);
						newPair.push_back(b);
						newFarPairs.push_back(newPair);
						break;
					}
				}
			}
		}

		FarPairs = newFarPairs;

	}

	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<bool> HeuristicLCDSBestIn(KGraph &g, long s)
{
	vector<bool> heuristicSoln(g.n, false);
	vector< vector<long> > FarPairs = EnumerateFarPairs(g, heuristicSoln, 1);
	vector<long> dist;
	vector<long> score;
	vector<long> scoreZero(g.n, 0);

	long a, b;
	long max;
	long imax;
	while (!FarPairs.empty()) 
	{
		score = scoreZero;
		max = -1;
		imax = -1;
		for (long i = 0; i < g.n; i++) 
		{
			if (heuristicSoln[i]) continue;
			dist = ComputeSSSPinGBv(g, heuristicSoln, i);
			for (long j = 0; j < FarPairs.size(); j++) 
			{
				a = FarPairs[j][0];
				b = FarPairs[j][1];
				if (dist[a] + dist[b] <= s)
				{
					score[i]++;
				}
				
			}

			if (score[i] > max) {
				max = score[i];
				imax = i;
			}
		}
		heuristicSoln[imax] = true;

		vector< vector<long> > newFarPairs;
		vector<long> newPair;
		dist = ComputeSSSPinGBv(g, heuristicSoln, imax);
		for (long j = 0; j < FarPairs.size(); j++)
		{
			a = FarPairs[j][0];
			b = FarPairs[j][1];
			if (dist[a] + dist[b] > s)
			{
				newPair.clear();
				newPair.push_back(a);
				newPair.push_back(b);
				newFarPairs.push_back(newPair);
			}
		}

		FarPairs = newFarPairs;
		
	}

	cerr << "Heuristic solution found,  number of nodes = " << count(heuristicSoln.begin(), heuristicSoln.end(), true) << endl;
	return heuristicSoln;
}

vector<long> BiasedShuffle(KGraph &g)
{
	vector<long> order(g.n,-1);
	long rn;
	long remdeg = 2*g.m; // remaining degree
	long tempdeg;
	vector<bool> selected(g.n, false);  // whether or not the vertex is in order yet.
	bool numFound;

	for (long i = 0; i < g.n; i++)	// randomly pick a vertex that has not been selected, biased by degree
	{
		rn = rand() % remdeg;

		numFound = false;
		tempdeg = 0;
		for (long j = 0; j < g.n && !numFound; j++)
		{
			if (selected[j]) continue;
			tempdeg += g.degree[j];
			if (rn < tempdeg)
			{
				order[g.n-i-1] = j;
				selected[j] = true;
				remdeg -= g.degree[j];
				numFound = true;
			}
		}
	}
	return order;
}
vector<bool> RandomHeuristicLCDS(KGraph &g, long s, long iterations)
{
	vector<bool> incumbent;
	long incsize = g.n+1;
	long newsize;
	vector<long> sortedList;

	for(long iter = 0; iter < iterations; iter++)
	{
		vector<bool> heuristicSoln(g.n, true);
		sortedList = BiasedShuffle(g);
		long v;

		for (int i = 0; i < g.n; i++) 
		{
			v = sortedList[i];
			heuristicSoln[v] = false;
			if (!IsLatencyConstrainedCDS(g, heuristicSoln, s))
			{
				heuristicSoln[v] = true;
			}
		}
		newsize = count(heuristicSoln.begin(), heuristicSoln.end(), true);
		if (newsize < incsize)
		{
			incumbent = heuristicSoln;
			incsize = newsize;
		}
	}

	cerr << "Heuristic solution found,  number of nodes = " << count(incumbent.begin(), incumbent.end(), true) << endl;
	return incumbent;
}


vector< vector<long> > FindBiconnectedComponents(KGraph &g, vector<bool> &AV)
{
	/* I tried to use the naming conventions presented in Tarjan's 1972 paper.
	I assume that the graph is connected, so that only one call to the recursion is necessary. */

	// declare vars
	long u = -1, v = 0, i = 0;
	vector<long> number(g.n, (long)-1);
	vector<long> lowopt(g.n, g.n);
	vector< vector<long> > BC;		// biconnected components
	stack<long> le, re;				// used to store a stack of edges. le is "left edge" and re is "right edge". An edge is stored (le[i],re[i]). 

									// perform DFS-based algorithm
	Bico_Sub(v, u, i, g, number, lowopt, le, re, BC);

	vector<long> countComp(g.n, 0);
	for (long p = 0; p<BC.size(); p++) // count how many components each vertex belongs to
		for (long q = 0; q<BC[p].size(); q++)
			countComp[BC[p][q]]++;

	vector<bool> AV_temp(g.n, false);
	AV = AV_temp;
	for (long p = 0; p<g.n; p++) // if a vertex belongs to >1 component, then it is a cut vertex
		if (countComp[p]>1)
			AV[p] = true;

	return BC;
}

void Bico_Sub(long v, long u, long &i, KGraph &g, vector<long> &number, vector<long> &lowopt, stack<long> &le, stack<long> &re, vector< vector<long> > &BC)
{
	i++;
	number[v] = i;
	lowopt[v] = number[v];
	long w;
	for (long j = 0; j<g.degree[v]; j++)
	{
		w = g.adj[v][j];
		if (number[w] == -1)
		{
			le.push(v);
			re.push(w);
			Bico_Sub(w, v, i, g, number, lowopt, le, re, BC);
			lowopt[v] = (long)min(lowopt[v], lowopt[w]);
			if (lowopt[w] >= number[v])
			{
				vector<long> temp_BC;
				vector<bool> bBC(g.n, false);
				while (!le.empty() && !re.empty() && number[le.top()] >= number[w])
				{
					bBC[le.top()] = true;
					bBC[re.top()] = true;
					le.pop();
					re.pop();
				}
				if (!le.empty() && le.top() == v)
				{
					bBC[le.top()] = true;
					bBC[re.top()] = true;
					le.pop();
					re.pop();
				}
				else
				{
					cerr << "ERROR: edge (v,w) not on top of stack" << endl;
				}
				for (long p = 0; p<g.n; p++)
					if (bBC[p])
						temp_BC.push_back(p);
				BC.push_back(temp_BC);
			}
		}
		else if (number[w]<number[v] && w != u)
		{
			le.push(v);
			re.push(w);
			lowopt[v] = min(lowopt[v], number[w]);
		}
	}
}

int myrandom(int i) 
{ 
	return std::rand() % i; 
}

vector<long> solveMCDS(KGraph &g, long s) 
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}
	try{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);

		model.update();

		cerr << "Adding objective function" << endl;
		for (int i = 0; i<g.n; i++)
		{
			X[i].set(GRB_DoubleAttr_Obj, 1);
		}
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();

		cerr << "Adding domination constraints" << endl;
		long v;
		vector<long> minimalCut;
		string MinimalizeSetting = "Basic";
		for (int i = 0; i<g.n; i++) // for each vertex i, make sure you choose at least one vertex from N(i). This is a vertex cut!
		{
			if(MinimalizeSetting == "None")
				minimalCut = MinimalizeNone(g, g.adj[i], s);
			else if(MinimalizeSetting == "Basic")
				minimalCut = MinimalizeBasic(g, g.adj[i], s);
			else if(MinimalizeSetting == "Fast")
				minimalCut = Minimalize(g, g.adj[i], s);   // a minimal subset of N(i) that is a length-s cut.
			else cerr << "ERROR: not a supported value for MinimalizeSetting." << endl;
			GRBLinExpr expr = 0;
			for (long j = 0; j < minimalCut.size(); j++)
			{
				v = minimalCut[j];
				expr += X[v];
			}
			model.addConstr(expr >= 1);
		}
		model.update();

		//Initial Solution
		double TotaltimeinHeuristic = 0;
		time_t start1 = clock();
		vector<bool> soln2 = HeuristicLCDSBestIn(g, s);
		cerr << "Best-in heuristic soln = " << count(soln2.begin(), soln2.end(), true) << endl;
		soln2 = HeuristicLCDS(g, soln2, s);
		cerr << "Best-in heuristic soln (after minimalizing) = " << count(soln2.begin(), soln2.end(), true) << endl;
		vector<bool> initialSolution = soln2;
		TotaltimeinHeuristic = (double)(clock() - start1) / CLOCKS_PER_SEC;

		for (long i = 0; i < g.n; i++)
		{
			if (initialSolution[i])
			{
				X[i].set(GRB_DoubleAttr_Start, 1.0);
			}
			else
			{
				X[i].set(GRB_DoubleAttr_Start, 0.0);
			}
		}

		// fix articulation vertices in solution
		vector<bool> ArticulationVertices; // (g.n, false);
		FindBiconnectedComponents(g, ArticulationVertices);

		for (int i = 0; i < g.n; i++) 
		{
			if (ArticulationVertices[i])
			{
				X[i].set(GRB_DoubleAttr_LB, 1.0);
			}
		}

		cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		model.update();

		//End of initial solution and articulation vertex cut

		cerr << "Adding lazy constraints (the vertex cut inequalities)\n";
		
		LazyConstraints cb = LazyConstraints (X, g, s);	// tell Gurobi which function generates the lazy cuts.
		model.setCallback(&cb);

		cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;
		
		cerr << endl;
		cerr << "Number of callbacks : " << LazyConstraints::numCallbacks << endl;
		cerr << "Time in callbacks : " << LazyConstraints::TotalCallbackTime << endl;
		cerr << "Number of lazy cuts : " << LazyConstraints::numLazyCuts << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;

		for (int i = 0; i<g.n; i++)
		{
			if (X[i].get(GRB_DoubleAttr_X)>0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				cds.push_back(i);
			}
		}

		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << count(initialSolution.begin(), initialSolution.end(), true) << " " << cds.size() << " " << lb << " " << NumBBNodes << " " << LazyConstraints::numLazyCuts << " " << TotaltimeinHeuristic << " " << LazyConstraints::TotalCallbackTime << " ";

		delete[] X;
		

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return cds;
}

vector<long> solveMCDSrel(KGraph &g, long s)
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}
	try {
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);

		model.update();

		cerr << "Adding objective function" << endl;
		for (int i = 0; i<g.n; i++)
		{
			X[i].set(GRB_DoubleAttr_Obj, 1);
		}
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();

		cerr << "Adding domination constraints" << endl;
		long v;
		vector<long> minimalCut;
		string MinimalizeSetting = "Fast";
		for (int i = 0; i<g.n; i++) // for each vertex i, make sure you choose at least one vertex from N(i). This is a vertex cut!
		{
			if (MinimalizeSetting == "None")
				minimalCut = MinimalizeNone(g, g.adj[i], s);
			else if (MinimalizeSetting == "Basic")
				minimalCut = MinimalizeBasic(g, g.adj[i], s);
			else if (MinimalizeSetting == "Fast")
				minimalCut = Minimalize(g, g.adj[i], s);   // a minimal subset of N(i) that is a length-s cut.
			else cerr << "ERROR: not a supported value for MinimalizeSetting." << endl;
			GRBLinExpr expr = 0;
			for (long j = 0; j < minimalCut.size(); j++)
			{
				v = minimalCut[j];
				expr += X[v];
			}
			model.addConstr(expr >= 1);
		}
		model.update();


		cerr << "Relaxing" << endl;

		model.relax().optimize();

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return cds;
}

vector<long> solveVBMCDS(KGraph &g, long s)
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}
	try {
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		//env.set(GRB_IntParam_Presolve, 0);
		//env.set(GRB_IntParam_Cuts, 0);
		GRBModel model = GRBModel(env);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar ***Y = new GRBVar**[g.n];
		for (long i = 0; i<g.n; i++)
		{
			GRBVar **Y_temp = new GRBVar*[g.n];
			for (long j = 0; j<g.n; j++)
				Y_temp[j] = model.addVars(s - 1, GRB_BINARY);
			Y[i] = Y_temp;
		}
		GRBVar ***Z = new GRBVar**[g.n];
		for (long i = 0; i<g.n; i++)
		{
			GRBVar **Z_temp = new GRBVar*[g.n];
			for (long j = 0; j<g.n; j++)
				Z_temp[j] = model.addVars(s - 2, GRB_BINARY);
			Z[i] = Z_temp;
		}
		model.update();

		cerr << "Adding objective function" << endl;
		for (int i = 0; i<g.n; i++)
		{
			X[i].set(GRB_DoubleAttr_Obj, 1);
		}

		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();

		cerr << "Adding s=2 constraints" << endl;
		for (long k = 0; k < g.n; k++)
		{
			for (long i = 0; i < g.n; i++)
			{
				if (i == k) continue;
				vector<long> CN = g.CommonNeighborsList(i, k);

				for (long v = 0; v < CN.size(); v++)
				{
					long j = CN[v];
					model.addConstr(X[j] <= Y[i][k][0]);   // constraints (2)
				}

				GRBLinExpr expr = 0;
				for (long v = 0; v < CN.size(); v++)
				{
					long j = CN[v];
					expr += X[j];
				}
				model.addConstr(Y[i][k][0] <= expr);    // constraint (3)
			}
		}

		cerr << "Adding s>=3 constraints" << endl;
		long k;
		for (long t = 1; t < s - 1; t++)
		{
			for (long j = 0; j < g.n; j++)
			{
				for(long p = 0; p < g.degree[j]; p++)
				{
					k = g.adj[j][p];
					for (long i = 0; i < g.n; i++)
					{
						if (i == j || i == k) continue;
						model.addConstr(Y[i][j][t-1] + X[j] <= Y[i][k][t] + 1); // constraint (4)
					}
				}
				
			}
		}

		for (long t = 1; t < s - 1; t++)
		{
			for (long k = 0; k < g.n; k++)
			{
				for (long i = 0; i < g.n; i++)
				{
					if (i == k) continue;
					GRBLinExpr expr = 0;
					for (long p = 0; p < g.degree[k]; p++)
					{
						long j = g.adj[k][p];
						expr += Z[i][j][t - 1];
					}
					model.addConstr(Y[i][k][t] <= expr);	// constraint (5)	
				}
			}
		}

		for (long t = 1; t < s - 1; t++)
		{
			for (long j = 0; j < g.n; j++)
			{
				for (long i = 0; i < g.n; i++)
				{
					if (i == j) continue;
					model.addConstr(Z[i][j][t-1] <= Y[i][j][t-1]);	// constraint (6)	
					model.addConstr(Z[i][j][t-1] <= X[j]);		// constraint (7)	
					model.addConstr(Y[i][j][t-1] + X[j] <= Z[i][j][t - 1] + 1);	// constraint (8)
				}
			}
		}

		long f;
		
		for (long j = 0; j < g.n; j++)
		{
			vector<bool> neighbors(g.n, false);
			neighbors[j] = true;
			for (long p = 0; p < g.degree[j]; p++)
			{
				f = g.adj[j][p];
				neighbors[f] = true;
			}
			for (long i = 0; i < g.n; i++)
			{
				if (neighbors[i]) continue;
				GRBLinExpr expr = 0;
				for (long t = 0; t < s-1; t++)
				{
					expr += Y[i][j][t];
				}
				model.addConstr(expr >= 1);	// constraint (9)
			}
		}
		for (long t = 0; t < s - 1; t++)
		{
			for (long j = 0; j < g.n; j++)
			{
				for (long i = 0; i < g.n; i++)
				{
					if (i == j) continue;
					model.addConstr( Y[i][j][t] == Y[j][i][t] );	// we know we can set y^t_{ij}=y^t_{ji}
				}
			}
		}

		model.update();

		double TotaltimeinHeuristic = 0;
		time_t start1 = clock();
		vector<bool> soln2 = HeuristicLCDSBestIn(g, s);
		cerr << "Best-in heuristic soln = " << count(soln2.begin(), soln2.end(), true) << endl;
		soln2 = HeuristicLCDS(g, soln2, s);
		cerr << "Best-in heuristic soln (after minimalizing) = " << count(soln2.begin(), soln2.end(), true) << endl;
		vector<bool> initialSolution = soln2;
		TotaltimeinHeuristic = (double)(clock() - start1) / CLOCKS_PER_SEC;
		cerr << "Initial solution = ";
		for (long i = 0; i < g.n; i++)
		{
			if (initialSolution[i])
			{
				X[i].set(GRB_DoubleAttr_Start, 1.0);
				cerr << i << " ";
			}
			else
			{
				X[i].set(GRB_DoubleAttr_Start, 0.0);
			}
		}
		cerr << endl;

		cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;

		cerr << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;

		for (int i = 0; i<g.n; i++)
		{
			if (X[i].get(GRB_DoubleAttr_X)>0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				cds.push_back(i);
			}
		}

		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << count(initialSolution.begin(), initialSolution.end(), true) << " " << cds.size() << " " << lb << " " << NumBBNodes << " " << "NA" << " " << TotaltimeinHeuristic << " " << "NA" << " ";

		delete[] X;

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return cds;
}

vector<long> solveVBMCDS2(KGraph &g, long s)
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}
	try {
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		//env.set(GRB_IntParam_Presolve, 0);
		//env.set(GRB_IntParam_Cuts, 0);
		GRBModel model = GRBModel(env);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar ***Y = new GRBVar**[g.n];
		for (long i = 0; i<g.n; i++)
		{
			GRBVar **Y_temp = new GRBVar*[g.n];
			for (long j = 0; j<g.n; j++)
				Y_temp[j] = model.addVars(s - 1, GRB_BINARY);
			Y[i] = Y_temp;
		}
		GRBVar ***Z = new GRBVar**[g.n];
		for (long i = 0; i<g.n; i++)
		{
			GRBVar **Z_temp = new GRBVar*[g.n];
			for (long j = 0; j<g.n; j++)
				Z_temp[j] = model.addVars(s - 2, GRB_BINARY);
			Z[i] = Z_temp;
		}
		model.update();
		cerr << "Adding objective function" << endl;
		for (int i = 0; i<g.n; i++)
		{
			X[i].set(GRB_DoubleAttr_Obj, 1);
		}

		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

		model.update();

		cerr << "Adding s=2 constraints" << endl;
		for (long k = 0; k < g.n; k++)
		{
			for (long i = 0; i < g.n; i++)
			{
				if (i == k) continue;
				vector<long> CN = g.CommonNeighborsList(i, k);

				for (long v = 0; v < CN.size(); v++)
				{
					long j = CN[v];
					model.addConstr(X[j] <= Y[i][k][0]);   // constraints (2)
				}

				GRBLinExpr expr = 0;
				for (long v = 0; v < CN.size(); v++)
				{
					long j = CN[v];
					expr += X[j];
				}
				model.addConstr(Y[i][k][0] <= expr);    // constraint (3)
			}
		}

		cerr << "Adding s>=3 constraints" << endl;
		long k;
		for (long t = 1; t < s - 1; t++)
		{
			for (long j = 0; j < g.n; j++)
			{
				for (long p = 0; p < g.degree[j]; p++)
				{
					k = g.adj[j][p];
					for (long i = 0; i < g.n; i++)
					{
						if (i == j || i == k) continue;
						model.addConstr(Y[i][j][t - 1] + X[j] <= Y[i][k][t] + 1); // constraint (4)
					}
				}

			}
		}

		for (long t = 1; t < s - 1; t++)
		{
			for (long k = 0; k < g.n; k++)
			{
				for (long i = 0; i < g.n; i++)
				{
					if (i == k) continue;
					GRBLinExpr expr = 0;
					for (long p = 0; p < g.degree[k]; p++)
					{
						long j = g.adj[k][p];
						expr += Z[i][j][t - 1];
					}
					model.addConstr(Y[i][k][t] <= expr);	// constraint (5)	
				}
			}
		}

		for (long t = 1; t < s - 1; t++)
		{
			for (long j = 0; j < g.n; j++)
			{
				for (long i = 0; i < g.n; i++)
				{
					if (i == j) continue;
					model.addConstr(Z[i][j][t - 1] <= Y[i][j][t - 1]);	// constraint (6)	
					model.addConstr(Z[i][j][t - 1] <= X[j]);		// constraint (7)	
					model.addConstr(Y[i][j][t - 1] + X[j] <= Z[i][j][t - 1] + 1);	// constraint (8)
				}
			}
		}

		long f;

		for (long j = 0; j < g.n; j++)
		{
			vector<bool> neighbors(g.n, false);
			neighbors[j] = true;
			for (long p = 0; p < g.degree[j]; p++)
			{
				f = g.adj[j][p];
				neighbors[f] = true;
			}
			for (long i = 0; i < g.n; i++)
			{
				if (neighbors[i]) continue;
				GRBLinExpr expr = 0;
				for (long t = 0; t < s - 1; t++)
				{
					expr += Y[i][j][t];
				}
				model.addConstr(expr >= 1);	// constraint (9)
			}
		}
		for (long t = 0; t < s - 1; t++)
		{
			for (long j = 0; j < g.n; j++)
			{
				for (long i = 0; i < g.n; i++)
				{
					if (i == j) continue;
					model.addConstr(Y[i][j][t] == Y[j][i][t]);	// we know we can set y^t_{ij}=y^t_{ji}
				}
			}
		}

		cerr << "Adding domination constraints" << endl;
		
		for (long i = 0; i < g.n; i++)
		{
			GRBLinExpr expr2 = 0;
			long v;
			for (long j = 0; j < g.degree[i]; j++)
			{
				v = g.adj[i][j];
				expr2 += X[v];
			}
			model.addConstr(expr2 >= 1);
		}


		model.update();

		double TotaltimeinHeuristic = 0;
		time_t start1 = clock();
		vector<bool> soln2 = HeuristicLCDSBestIn(g, s);
		cerr << "Best-in heuristic soln = " << count(soln2.begin(), soln2.end(), true) << endl;
		soln2 = HeuristicLCDS(g, soln2, s);
		cerr << "Best-in heuristic soln (after minimalizing) = " << count(soln2.begin(), soln2.end(), true) << endl;
		vector<bool> initialSolution = soln2;
		TotaltimeinHeuristic = (double)(clock() - start1) / CLOCKS_PER_SEC;
		cerr << "Initial solution = ";
		for (long i = 0; i < g.n; i++)
		{
			if (initialSolution[i])
			{
				X[i].set(GRB_DoubleAttr_Start, 1.0);
				cerr << i << " ";
			}
			else
			{
				X[i].set(GRB_DoubleAttr_Start, 0.0);
			}
		}
		cerr << endl;

		cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;

		for (int i = 0; i<g.n; i++)
		{
			if (X[i].get(GRB_DoubleAttr_X)>0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				cds.push_back(i);
			}
		}

		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << count(initialSolution.begin(), initialSolution.end(), true) << " " << cds.size() << " " << lb << " " << NumBBNodes << " " << "NA" << " " << TotaltimeinHeuristic << " " << "NA" << " ";

		delete[] X;

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return cds;
}


vector<long> solveMCDSweighted(KGraph &g, long s, vector<long> W) {
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}
	try {
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		model.update();

		cerr << "Adding objective function" << endl;
		for (int i = 0; i<g.n; i++)
		{
			X[i].set(GRB_DoubleAttr_Obj, 1);
		}
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		//model.set(GRB_DoubleParam_MIPGapAbs, 0.9);
		model.update();

		cerr << "Adding domination constraints" << endl;
		long v;
		vector<long> minimalCut;
		for (int i = 0; i<g.n; i++) // for each vertex i, make sure you choose at least one vertex from N(i). This is a vertex cut!
		{
			minimalCut = MinimalizeWeighted(g, g.adj[i], s, W);   // a minimal subset of N(i) that is a length-s cut.
			GRBLinExpr expr = 0;
			for (long j = 0; j < minimalCut.size(); j++)
			{
				v = minimalCut[j];
				expr += X[v];
			}
			model.addConstr(expr >= 1);
		}
		model.update();


		// fix articulation vertices in solution
		vector<bool> ArticulationVertices; 
		FindBiconnectedComponents(g, ArticulationVertices);

		for (int i = 0; i < g.n; i++)
		{
			if (ArticulationVertices[i])
			{
				X[i].set(GRB_DoubleAttr_LB, 1.0);
			}
		}
		cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		model.update();

		//End of initial solution and articulation vertex cut

		cerr << "Adding lazy constraints (the vertex cut inequalities)\n";

		LazyConstraints3 cb = LazyConstraints3(X, g, s, W);	// tell Gurobi which function generates the lazy cuts.
		model.setCallback(&cb);

		cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;

		cerr << endl;
		cerr << "Number of callbacks : " << LazyConstraints3::numCallbacks << endl;
		cerr << "Time in callbacks : " << LazyConstraints3::TotalCallbackTime << endl;
		cerr << "Number of lazy cuts : " << LazyConstraints3::numLazyCuts << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;

		for (int i = 0; i<g.n; i++)
		{
			if (X[i].get(GRB_DoubleAttr_X)>0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				cds.push_back(i);
			}
		}

		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << "NA" << " " << cds.size() << " " << lb << " " << NumBBNodes << " " << LazyConstraints3::numLazyCuts << " " << "NA" << " " << LazyConstraints3::TotalCallbackTime << " ";
		delete[] X;

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return cds;

}

vector<long> solveReMCDS(KGraph &g, long s)
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}

	// check if the graph has any articulation vetices
	vector<bool> ArticulationVertices; // (g.n, false);
									   //vector< vector<long>> Z = FindBiconnectedComponents (g, ArticulationVertices);
	FindBiconnectedComponents(g, ArticulationVertices);
	long NumArticulationVertices = count(ArticulationVertices.begin(), ArticulationVertices.end(), true);
	if (NumArticulationVertices > 0)
	{
		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " ";
		cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		return cds;
	}

	else {
		try {
			GRBEnv env = GRBEnv();
			env.set(GRB_IntParam_OutputFlag, 0);
			env.set(GRB_DoubleParam_TimeLimit, 3600);
			//env.set(GRB_DoubleParam_MIPGapAbs, 0.01);
			GRBModel model = GRBModel(env);
			model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
			GRBVar *X = model.addVars(g.n, GRB_BINARY);
			model.update();

			cerr << "Adding objective function" << endl;
			for (int i = 0; i < g.n; i++)
			{
				X[i].set(GRB_DoubleAttr_Obj, 1);
			}
			model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
			model.update();

			cerr << "Adding domination constraints" << endl;
			long v;
			vector<long> minimalVertexCut;
			vector<long> minimalCut;
			for (int i = 0; i < g.n; i++) // for each vertex i, make sure you choose at least one vertex from N(i). This is a vertex cut!
			{
				//Adding vertex cut constraints
				minimalVertexCut = Minimalize(g, g.adj[i], g.n-1);
				GRBLinExpr expr = 0;
				for (long j = 0; j < minimalVertexCut.size(); j++)
				{
					v = minimalVertexCut[j];
					expr += X[v];
				}
				model.addConstr(expr >= 2);
				//Adding length-s cut constraints 
				minimalCut = Minimalize(g, g.adj[i], s);
				GRBLinExpr expr1 = 0;
				for (long j = 0; j < minimalCut.size(); j++)
				{
					v = minimalCut[j];
					expr1 += X[v];
				}
				model.addConstr(expr1 >= 1);
			}

			model.update();

			//End of initial solution and articulation vertex cut

			cerr << "Adding lazy constraints (the vertex cut inequalities)\n";

			LazyConstraints4 cb = LazyConstraints4(X, g, s);	// tell Gurobi which function generates the lazy cuts.
			model.setCallback(&cb);

			cerr << "Optimizing" << endl;
			model.optimize();

			long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
			double gap = model.get(GRB_DoubleAttr_MIPGap);
			cerr << "MIP gap = " << gap << endl;

			cerr << endl;
			cerr << "Number of callbacks : " << LazyConstraints2::numCallbacks << endl;
			cerr << "Time in callbacks : " << LazyConstraints2::TotalCallbackTime << endl;
			cerr << "Number of lazy cuts : " << LazyConstraints2::numLazyCuts << endl;
			cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

			cerr << "Extracting solution" << endl;
			int status = model.get(GRB_IntAttr_Status);
			cerr << statusNumtoString(status) << endl;

			for (int i = 0; i < g.n; i++)
			{
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
				{
					cds.push_back(i);
				}
			}

			cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << "NA" << " " << cds.size() << " " << gap << " " << NumBBNodes << " " << LazyConstraints4::numLazyCuts << " " << "NA" << " " << LazyConstraints4::TotalCallbackTime << " ";

			delete[] X;

		}

		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during optimization" << endl;
		}
		return cds;
	}
}

vector<long> solveROBMCDS(KGraph &g, long s)
{
	vector<long> cds;
	if (!g.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No CDS exists! Graph is not connected." << endl;
		return cds;
	}

	// check if the graph has any articulation vetices
	vector<bool> ArticulationVertices; // (g.n, false);
									   //vector< vector<long>> Z = FindBiconnectedComponents (g, ArticulationVertices);
	FindBiconnectedComponents(g, ArticulationVertices);
	long NumArticulationVertices = count(ArticulationVertices.begin(), ArticulationVertices.end(), true);
	if (NumArticulationVertices > 0) 
	{
		cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " " << "N/A" << " ";
		cerr << "***Number of articulation vertices = " << count(ArticulationVertices.begin(), ArticulationVertices.end(), true) << endl;
		return cds;
	}

	else {
		try {
			GRBEnv env = GRBEnv();
			env.set(GRB_IntParam_OutputFlag, 0);
			env.set(GRB_DoubleParam_TimeLimit, 3600);
			//env.set(GRB_DoubleParam_MIPGapAbs, 0.01);
			GRBModel model = GRBModel(env);
			model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
			GRBVar *X = model.addVars(g.n, GRB_BINARY);
			model.update();

			cerr << "Adding objective function" << endl;
			for (int i = 0; i < g.n; i++)
			{
				X[i].set(GRB_DoubleAttr_Obj, 1);
			}
			model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
			model.update();

			cerr << "Adding domination constraints" << endl;
			long v;
			vector<long> minimalVertexCut;
			vector<long> minimalCut;
			for (int i = 0; i < g.n; i++) // for each vertex i, make sure you choose at least one vertex from N(i). This is a vertex cut!
			{
				//Adding length-s cut constraints 
				minimalCut = Minimalize(g, g.adj[i], s);
				GRBLinExpr expr1 = 0;
				for (long j = 0; j < minimalCut.size(); j++)
				{
					v = minimalCut[j];
					expr1 += X[v];
				}
				model.addConstr(expr1 >= 2);
			}
			model.update();

			//initail solution and articulation vertex cut
			vector <bool> initialSolution(g.n, false);
			double TotaltimeinHeuristic = 0;
			time_t start1 = clock();
			vector<bool> soln2 = HeuristicLCDSBestIn(g, s);
			cerr << "Best-in heuristic soln (for r=1) = " << count(soln2.begin(), soln2.end(), true) << endl;
			soln2 = HeuristicRLCDS(g, soln2, s);
			cerr << "Best-in heuristic soln (for r=2, after augmenting the solution) = " << count(soln2.begin(), soln2.end(), true) << endl;
			soln2 = MinimalizeRLCDS(g, soln2, s);
			cerr << "Best-in heuristic soln (for r=2, after minimalizaing the solution) = " << count(soln2.begin(), soln2.end(), true) << endl;
			initialSolution = soln2;
			TotaltimeinHeuristic = (double)(clock() - start1) / CLOCKS_PER_SEC;
			cerr << "Is it actually a feasible solution? " << IsLatencyConstrainedRCDS(g, initialSolution, s) << endl;
			for (long i = 0; i < g.n; i++)
			{
				if (initialSolution[i])
				{
					X[i].set(GRB_DoubleAttr_Start, 1.0);
				}
				else
				{
					X[i].set(GRB_DoubleAttr_Start, 0.0);
				}
			}


			//cerr << "***Number of length-s cut vertices = " << lengthsAV << endl;
			model.update();

			//End of initial solution and articulation vertex cut

			cerr << "Adding lazy constraints (the vertex cut inequalities)\n";

			LazyConstraints2 cb = LazyConstraints2(X, g, s);	// tell Gurobi which function generates the lazy cuts.
			model.setCallback(&cb);

			cerr << "Optimizing" << endl;
			model.optimize();

			long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
			double lb = model.get(GRB_DoubleAttr_ObjBound);
			cerr << "MIP LB = " << lb << endl;

			cerr << endl;
			cerr << "Number of callbacks : " << LazyConstraints2::numCallbacks << endl;
			cerr << "Time in callbacks : " << LazyConstraints2::TotalCallbackTime << endl;
			cerr << "Number of lazy cuts : " << LazyConstraints2::numLazyCuts << endl;
			cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

			cerr << "Extracting solution" << endl;
			int status = model.get(GRB_IntAttr_Status);
			cerr << statusNumtoString (status) << endl;
			
			for (int i = 0; i < g.n; i++)
			{
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
				{
					cds.push_back(i);
				}
			}
			

			cout << g.name << " " << g.n << " " << g.m << " " << g.DiameterUnweighted() << " " << s << " " << count(initialSolution.begin(), initialSolution.end(), true) << " " << cds.size() << " " << lb << " " << NumBBNodes << " " << LazyConstraints2::numLazyCuts << " " << TotaltimeinHeuristic << " " << LazyConstraints2::TotalCallbackTime << " ";

			delete[] X;

		}

		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during optimization" << endl;
		}
		return cds;
	}
}
