/*
 * mwm.h
 *
 *  Created on: 10.06.2015
 *      Author: droschin
 */

#ifndef MWM_H
#define MWM_H
#define N_MATCHING_DEBUG

#include <global.h>
#pragma warning(push, 0)
#include <ogdf/basic/Graph.h>
#pragma warning(pop)
#include <set>
#include <vector>
#include <stack>
#ifdef GRAPHICS
#include <SFML/Graphics.hpp>
#endif

using ogdf::Graph;
using ogdf::edge;
using ogdf::node;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using std::vector;
using std::stack;

class MWM {
private:
	unsigned m_nodesA, m_nodesJ; // Number of vertices in each bipartite set
	Graph m_G; // bipartite graph
	vector<node> m_agent,m_job; // bipartite vertex sets
	NodeArray<node> m_vMatchingTovOriginal; // mapping of matching vertices to original (input) vertices
	EdgeArray<wType> m_egdeweight; // weights of the edges
	NodeArray<node> m_mate; // adjacent matching vertices
	NodeArray<wType> m_dual; // dual variables
	NodeArray<node> m_copy; // vertex copies
	enum NodeType {AGENT, JOB, AGENTCOPY, JOBCOPY};
	NodeArray<NodeType> m_nodeType;
#ifdef GRAPHICS
	sf::Font m_font_Arial;
#endif

	int m_computedMatchings = -1; // -1: no matching; 0..m_nodesJ-1: sub matching only; m_nodesJ: all matchings
	vector<wType> m_weight; // calculated weights of the MWMs

	bool m_withEnumeration; // true, if enumeration is to be handled
	vector<Graph> m_EnumerationGraph; // enumeration graphs
	vector<NodeArray<node>> m_vEnumTovOriginal; // Original vertices of the enumeration nodes
	vector<EdgeArray<bool>> m_EnumMatched; // true, if the edge of the enumeration graph is a matching edge
	vector<NodeArray<NodeType>> m_EnumNodeType; // Type of vertex
	vector<NodeArray<node>> m_EnumCopy; // Vertex copies
	vector<vector<node>> m_EnumNonIsolatedVertices; // all the non isolated vertices at the start of the recursion
	vector<unsigned> m_EnumNumberofNonIsolatedVertices;  // only the first m_EnumNumberofNonIsolatedVertices are non isolated in the current recursion
	struct EnumRestore // needed to restore changes made during the computation of the MWMs
	{
		//Graph::HiddenEdgeSetHandle hesh;
		Graph::HiddenEdgeSet hes;
		vector<edge> shortcutedges; // newly inserted edges
		vector<edge> circle; // alternating cycle
		vector<edge> trimmedEdges; // trimmed edges
		vector<edge> edgeSet2; // edges of second set
		//vector<node> isolatedVertices;
		unsigned removedIsolatedVertices=0;
		int additionalFixedMatchingEdges=0;
		int f=0; // f value of graph
		EnumRestore(Graph &G):hes(G)  {
			//hesh=G.newHiddenEdgeSet();
			}
		~EnumRestore() {};
	};
	vector<vector<EnumRestore*>> m_EnumRestore; // Information to restore graph affected by trimming and alternating cycle
	vector<EdgeArray<bool>> m_EnumDeg2TrimmedEdges; // Edges on a degree 2 path, that were trimmed
	vector<EdgeArray<vector<edge>*>> m_EnumDeg2ReplacedPath; // vector of edges, which are represented by the edge in the EdgeArray
#ifdef GRAPHICS
	vector<sf::RenderWindow*> m_EnumRenderWindow; // Graphical display of the enumeration matchings
#endif
	vector<vector<edge>> m_EnumFixedMatching; // currently fixed matching edges (matching edges on no directed cycle)
	enum StackOP {ENUM, PREP_M1, PREP_M2, UNDO};
	vector<stack<StackOP>> m_EnumStack;
	vector<NodeArray<short>> m_EnumDFSstatus; // Dfs Status. status 0=not visited, 1=on stack, 2=visited, 3=completed (used in computeSCCs) -initialized here, to avoid O(n) initialization time in each step of computeSSCs
	vector<NodeArray<unsigned>> m_EnumSCC; // Number of strongly connected component of the nodes
	vector<EdgeArray<edge>> m_EnumPreviousEdgeOnCircle; // Number of strongly connected component of the nodes
	vector<unsigned> m_EnumMaxSCCNumber; // Number of strongly connected component of the nodes
	vector<int> m_EnumEdgeID; // keep edge IDs as low as possible

	bool m_EnumStoreMatchings=true; // want to store matchings
	vector<bool> m_EnumStoredMatchings; // all matchings computed and stored
	vector<unsigned> m_EnumCurStoredMatching;
	vector<vector<vector<edge>*>> m_EnumStoredMatching;

	bool m_graph_copied=false; // true, if vertices have been copied
	vector<unsigned> m_agentsInComp, m_jobsInComp; // number of vertices in each component of the full matching graph

public:
	struct AgentJobPair
	{
		node matchingAgent;
		node matchingJob;
		AgentJobPair(node matchingAgent,node matchingJob):matchingAgent(matchingAgent), matchingJob(matchingJob) {}
	};
private:
	vector<vector<AgentJobPair>> m_matchingAgentAndJob; // computed initial matching edges of all bipartite graph instances
	vector<unsigned> m_refToMatchings; // Reference i for the above variables, e.g.,  weight[refToMatchings[i]]

	void addGraphCopy(); // Create duplicates of given nodes, and edges between them, if necessary
	void compute(unsigned jobIndex); // if jobIndex=m_nodesJ, computes main and sub MWMs; else if m_computedMatchings=-2 compute sub MWM else compute main and sub MWMs

	void constructEqualitySubgraph(const node jobToRemove, unsigned nodesInEachSet); // construct equality subgraph based on m_G and m_dual, appends it in m_EnumerationGraph
	// compute strongly connected components and hide edges between different components as well as isolated vertices
	//void trim(Graph &G, const EdgeArray<bool>& matched, const NodeArray<NodeType>& nodeType, const NodeArray<node>& copy, vector<node>& nonIsolatedVertex, vector<Graph::HiddenEdgeSetHandle>& HESH);
	void trim(unsigned jobIndex, const vector<edge>* edgeToRemove=nullptr); // trims the graph according to current matching as desribed in Uno's algorithm. optional edgeToRemove are edges of E1 or E2
	int addEgdeToFixedMatching(unsigned jobIndex, const edge e, bool matchingFlag);
	void computeSCCs(unsigned jobIndex);
	void getOtherEdgeAndNode(const node n, const edge e, node &othernode, edge &otheredge); // get other edge of node n, which is not e. Requires deg(n)=2. Other node and edge is stored in othernode/edge
	void getOtherEdgeAndNodeWithDifferentMatchedValue(const EdgeArray<bool>& matched, const node n, const edge e, node &othernode, edge &otheredge); // get another edge of node n, which is not e. Other node and edge is stored in othernode/edge. ensures matched[e]!=matched[otheredge]


public:
	MWM(bool withEnumeration=false);
	~MWM();

	void addAgent(node agent); // Add agent to bipartite graph.
	void addJob(node job); // Add job to bipartite graph.
	void addEdge(unsigned IndexAgent, unsigned IndexJob, wType weight); // Adds edge with given weight. Do this after creating duplicates only! todo: remove this restriction

	wType getWeight(unsigned jobIndex);// { return m_weight[m_refToMatchings[jobIndex]]; }
	wType getFullWeight() { return getWeight(m_nodesJ); }
	bool computedAllMatchings() const { return m_computedMatchings== (int) m_nodesJ; } // true if all matchings have been computed
	int computedMatchings() const { return m_computedMatchings; }

	const vector<AgentJobPair>& getMatching(unsigned jobIndex) { getWeight(jobIndex); return m_matchingAgentAndJob[m_refToMatchings[jobIndex]];  }
	const vector<AgentJobPair>& getFullMatching() { return getMatching(m_nodesJ);  }
	void displayEnum(unsigned jobIndex); // graphical display of enumeration graph
	vector<edge>* enumNextMatching(unsigned jobindex); // enumerate next MWM. Return pointer to matching edge vector. If there is none, return nullptr. Then next call returns first matching
	const NodeArray<node>& getvEnumTovOriginal(unsigned jobIndex) const { return m_vEnumTovOriginal[jobIndex]; }
};

#endif /* MWM_H */
