#ifndef BBPMCSI_H
#define BBPMCSI_H

#include <global.h>
#include <mwm.h>
#pragma warning(push, 0)
#include <ogdf/decomposition/BCTree.h>
#pragma warning(pop)
#include <vector>
#include <forward_list>
#ifdef GRAPHICS
#include <SFML/Graphics.hpp>
using namespace sf;
#endif
#pragma warning(push, 0)
#include <ogdf/planarity/PlanarizationGridLayout.h>
#pragma warning(pop)
#include <graphVisualization.h>
#include <MCIS2.h>
#include <auxiliary.h>

using namespace ogdf;
using namespace std;

enum LaWeCSExpand
{
	NO_EXPANSION,
	SKIPPED_VERTEX,
	MWM_EXPANSION,
	MWM_AND_SKIPPED_VERTEX
};

class SkippedVertex
{
private:
	vector<wType> m_weight; // calculated weights for the skipped vertices
	node m_best_aH = nullptr, m_2nd_best_aH = nullptr; // child of H with best and 2nd best solution
	wType m_best_weight = WEIGHT_NOT_COMPATIBLE, m_2nd_best_weight = WEIGHT_NOT_COMPATIBLE; // their weights
	Mapping m_best_mapping, m_2nd_best_mapping; // their mappings

	int m_computed_instance = -1; // -1: no instance; 0..k: sub instance k; -2: all sub instances; only for children of H
	node m_skipped_child = nullptr;// skipped block bH of the first sub instance of skipped vertices

public:
	SkippedVertex(unsigned neighbours);
	~SkippedVertex() { }
	void updateStoredBestSolutions(node aH, wType weight, const Mapping& mapping);
	wType getWeight(unsigned jobIndex) const;
	void setWeight(unsigned jobIndex, wType weight);
	int getComputedInstance() const { return m_computed_instance; }
	void setComputedInstance(int computed_instance) { m_computed_instance = computed_instance; }
	node getSkippedChild() const { return m_skipped_child; }
	void setSkippedChild(node skippedaH) { m_skipped_child = skippedaH; }
	wType getBestWeightAndMapping(node excluded, Mapping &mapping) const;
};

class BBP_MCSI
{
private:
	LabelFunction &m_labelFunction;
	InputGraph *m_ig_G,*m_ig_H;
	BCTree *m_BC_G,*m_BC_H; // BC-Trees of G and H
	const Graph *m_G,*m_H,*m_BCG,*m_BCH,*m_AG,*m_AH; // Pointers to the graphs, the BC-trees and the auxiliary graphs
	node m_b_initial; // Initial b-node of BC_G
	const NodeArray<labelType> *m_nodeLabelSimpleG=nullptr,*m_nodeLabelSimpleH=nullptr;
	const EdgeArray<labelType> *m_edgeLabelSimpleG=nullptr,*m_edgeLabelSimpleH=nullptr;
	const vector<string> *m_SimpleLabelToString;
	bool m_simpleLabelG, m_simpleLabelH, m_enumerate;

	NodeArray<vector<node>>& m_VbcG; // vertices of the auxiliary graph in a B oder C-node of BCG
	NodeArray<vector<node>>& m_VbcH; // vertices of the auxiliary graph in a B oder C-node of BCH
	NodeArray<Graph>& m_bGAG; // the single parts of the auxiliary graph of G
	NodeArray<Graph>& m_bHAH; // the single parts of the auxiliary graph of H
	NodeArray<node>& m_aGtoSPaG; // mapping from aG to the nodes of the single parts of AG
	NodeArray<node>& m_aHtoSPaH; // mapping from aH to the nodes of the single parts of AH
	NodeArray<NodeArray<node>>& m_SPaGtoaG; // mapping from the nodes of the single parts of AG to aG
	NodeArray<NodeArray<node>>& m_SPaHtoaH; // mapping from the nodes of the single parts of AH to aG
	NodeArray<EdgeArray<edge>>& m_SPeGtoaeG; // mapping from the edges of the single parts of AG to the edges aG
	NodeArray<EdgeArray<edge>>& m_SPeHtoaeH; // mapping from the edges of the single parts of AH to the edges aH
	NodeArray<bool>& m_bGOuterplanar; // true if bG is outerplanar
	NodeArray<bool>& m_bHOuterplanar; // true if bH is outerplanar
	NodeArray<StaticPlanarSPQRTree*>& m_bGToSPQRTree; // bG to SPQR tree
	NodeArray<StaticPlanarSPQRTree*>& m_bHToSPQRTree; // bH to SPQR tree
	NodeArray<node> m_parent_BC_G; // parent node of each vertex in the BC tree of G, where m_parent_BC_G[m_b_initial]=m_b_initial

	NodeArray<NodeArray<MaxComPart2Tree*>> m_mcp2txG; // mappings between two blocks without given mapping
	NodeArray<NodeArray<MaxComPart2Tree*>> m_mcp2tvG; // mappings between two blocks with given mapping vG to some other CV

	vector<node> m_CVH; // Contains all the CVs of H
	wType m_weightBBPIso; // Computed weight of the isomorphism
	clock_t m_maxCliqueTime; // maximum time per clique mcs computation
	unsigned m_unfinishedCliqueComputations=0; // number of stopped clique computations due to exceeded time
	wType m_distancePenalty; // penalty for skipped vertices WEIGHT_NOT_COMPATIBLE

#ifdef GRAPHICS
	// graphical display
	GraphVisualization *m_graphVisualizationG=nullptr,*m_graphVisualizationH=nullptr;
	bool m_graphical_display;
	double m_lastscaling=-1;
	RenderWindow m_window;
#endif

	// Maximum weight matchings
	NodeArray<NodeArray<MWM*>> m_Matching; // Set of MWMs for a pair of CVs. Parameters are vG, vH.

	//void computeMatchings(const node bG, const node xG=nullptr); // computes all occurring Maximum Weight Matchings
	wType getMatchingValue(const node vG, const node aH); // return weight of a MWM for CVs of G and H. B-Node of aH is excluded, if aH is a vertex of a B-node.
	const vector<MWM::AgentJobPair>& getMatching(const node vG, const node aH) const; // return MWM for CVs of G and H. B-Node of aH is excluded, if aH is a vertex of a B-node.

	// Skipped vertices
	NodeArray<NodeArray<SkippedVertex*>> m_SkippedVertex; // Structure analog to MWMs for skipped vertices.

	 // return weight of skipped vertex LAWECS for CVs of G and H. B-Node of aH is excluded, if aH is a vertex of a B-node.
	wType getSkippedVertexValue(const node vG, const node aH);

	// compute the weight of a matching edge for blocks bG, bH with fixed mapping vG->vH. Separate weights for expansion through matching and skipped vertex are also computed
	wType computeMatchingEdgeWeight(node bG, node bH, node vG, node vH, wType &wMatching, wType &wSkippedVertex);

	NodeArray<NodeArray<LaWeCSExpand>> m_LaWeCSExpand; // stores expansion type along pair of vertices vG, aH - only if for corresponding blocks bG, bH are bridges!
	NodeArray<NodeArray<Mapping>> m_LaWeCSMapping; // stores the closest non skipped vertices of the maximum skipped solution and if further expansion necessary

	NodeArray<wType> m_LaWeCS_SkippedRootWeight; // stores the weight of a LAWECS with first vertex vG as inner root vertex
	NodeArray<vector<Mapping>> m_LaWeCS_SkippedRootMapping; // stores two end points and possible expansion for above LAWECS.

	// Structure for sub instance before main instance computation (see LAWECSu paper)
	NodeArray<NodeArray<unsigned>> m_submatching_index; // index of skipped block bH of the first sub matching. Parameters are vG, vH.
	NodeArray<NodeArray<node>> m_submatching_block; // skipped block bH of the first sub instance. Parameters are vG, vH.
	NodeArray<unsigned> m_instance_ID; // Internal references for the main and sub instances. Parameter is aH.

	// Computed weights of BBP_Edge and BBP_SingleVertex
	NodeArray<NodeArray<wType>> m_Weight_BBPE_bGTobH; // Computed weights of BBP_Egde (at least one edge), where vG, vG are nullptrs (call from SetSx). Parameters are bG, bH.
	NodeArray<NodeArray<wType>> m_Weight_BBPE_aGToaH; // Computed weights of BBP_Egde, where vG, vH are set. Parameters are auxiliary vertices aG, aH.
	NodeArray<NodeArray<wType>> m_Weight_BBPSV_vGTovH; // Computed weights of BBP_SingleVertex. Parameters are vG,vH.

	// Computed local maximum parts of BBP_Egde, where vG, vH are nullptr (call from SetSx). Parameters are bG, bH. Each list entry contains a maximum local mapping. Without enumeration only one is stored.
	NodeArray<NodeArray<forward_list<vector<Mapping>>>> m_LocalIso_BBPE_bGTobH;
	NodeArray<NodeArray<bool>> m_LocalIso_BBPE_bGTobH_first;
	NodeArray<NodeArray<forward_list<vector<Mapping>>>> m_LocalIso_BBPE_aGToaH; // LocalIso of BBP_Egde, where vG, vH are set. Parameters are auxiliary vertices aG, aH.

	// a) nullptr: single vertex mapping with possible expansion into other blocks; b) Pointer to a maximum m_LocalIso_BBPE_aGToaH vector
	// Note: Enumeration of different aH blocks through code
	NodeArray<NodeArray<vector<Mapping>*>> m_LocalIso_BBPSV_vGTovH;

	// Enumeration stuff
	enum class expansionType {BLOCK_BRIDGE_BG, BLOCK_BRIDGE_AG, MATCHING};
	struct EnumStackElement
	{
		expansionType type;
		node nG,nH; // depending on expansionType. bG,bH for BLOCK_BRIDGE_BG; aG,aH for BLOCK_BRIDGE_AG; vG,aH for MATCHING
		unsigned num_elements_in_hierarchy; // number of elements in the current hierarchy. Only set for first element
		EnumStackElement(expansionType type, node nG, node nH, unsigned num_elements_in_hierarchy):type(type),nG(nG),nH(nH),num_elements_in_hierarchy(num_elements_in_hierarchy) {}
	};
	stack<EnumStackElement> m_EnumStackInProgressToMaximize;
	stack<EnumStackElement> m_EnumStackToEnumerate;
	vector<node> m_EnumMappingSourceNode, m_EnumMappingTargetNode; // current BBP-MCSI vertex mapping during enumeration
	vector<edge> m_EnumMappingSourceEdge, m_EnumMappingTargetEdge; // current BBP-MCSI edge mapping during enumeration
	NodeArray<NodeArray<forward_list<vector<Mapping>>::iterator>> m_EnumBlockBridgeMappingIterator_bGTobH; // list iterator of vector<Mapping> bGtobH
	NodeArray<NodeArray<forward_list<vector<Mapping>>::iterator>> m_EnumBlockBridgeMappingIterator_aGToaH; // list iterator of vector<Mapping> aGtoaH
	//NodeArray<NodeArray<vector< forward_list<vector<Mapping>> * >::iterator>> m_EnumBlockBridgeSV_Iterator; // vector iterator to different m_LocalIso_BBPE_aGToaH forward_lists - not needed, enumeration through code

	unsigned m_mcs_number=0;
	//const vector<edge>* enumNextMatching(const node vG, const node aH); // enumerate next matching

	void compute_BC_G_parents(node v); // compute the DFS parents in the m_BCG with root m_b_initial; needed to compute the MWMs
	wType setSx(); // Computes weight of maximum isomorphism on CC(V(G)\{xG}, V(bG))
	// Weight of maximum isomorphism phi, where at least one edge of B-Node bG is mapped to B-Node bH,
	//  restricted to CC(V(G)\{xG},V(bG)) and phi(vG)=vH
	// bG, bH are B-Nodes; xG is a possible excluded vertex in G; vG, vH are vertices of G,H
	wType bbpEdge(const node bG, const node bH, const node xG, const node vG=nullptr, const node vH=nullptr);
	// auxiliary function for bbpEdge, since there are 8 possibilities for mappings with skipped vertices
	// bG and bH are the blocks
	// mapping is eG.source()->eH.source(), if !flipEdges, else eG.source()->eH.target(), other edge vice versa
	// skipvG is true, if that mapping from eG.Source() is skipped, analog for skipTwinG
	// if the computed solution is larger than the previous, it is stored in m_LocalIso_BBPE_bGTobH[bG][bH] and m_Weight_BBPE_bGTobH[bG][bH]
	void bbpEdgeBridgebGbH(const node bG, const node bH, const edge eG, const edge eH, bool flipEdges, bool skipvG, bool skipTwinG);

	// Weight of maximum isomorphism phi, where exactly one vertex of B-Node bG is mapped, phi(vG)=vH
	// bG is a B-Node; vG, vH are vertices of G,H
	wType bbpSingleVertex(const node bG, const node vG, const node vH);
	void outputIsomorphism(vector<Mapping>& localMapping);
	void maximizeIsomorphism(); // maximize current solution based on m_EnumStackInProgressToMaximize
	void outputIsomorphism(); // output computed solution stored in m_EnumMappingSourceNode, m_EnumMappingTargetNode, (m_EnumMappingSourceEdge, m_EnumMappingTargetEdge), as well as graphical display
	void enumerateIsomorphisms(); // enumerate all the next solutions based on ... m_EnumStackFinished (?)
	bool enumerateBBPE_bGTobH(const node bG, const node bH); // enumerates next bridge/block bG to bH mapping with at least one mapped edge. Returns false if there is none.
	bool enumerateBBPE_aGToaH(const node aG, const node aH); // enumerates next bridge/block aG to aH mapping with at least one mapped edge. Returns false if there is none.
	void enumerateBbpEdgeaG(const node aG, const node aH); // enumerates next bridge  - maybe replace through previous method
	void enumerateBbpSingleVertex(const node vG, const node vH);
	void removeEnumMapping(int numNodes);

	// Weight of mapping between two vertices / edges of the graphs G,H
	wType w(const node vG, const node vH) const;
	wType w(const edge eG, const edge eH) const;

	// Check if mapping between edges and vertices is compatible
	bool compatible(const node vG, const node vH) const { return w(vG,vH)>=0; }
	bool compatible(const edge eG, const edge eH) const { return w(eG,eH)>=0; }

	void display(); // display the graphs, called from computeIsomorphism

public:
	wType getWeight(const node aG, const node aH) const { return w(m_BC_G->original(aG),m_BC_H->original(aH)); } // only for MCIS2
	wType getWeightPlusCVMatching(const node vG, const node aH, bool &expand); // get weight of matching (if there is one, then expand is set to true, else false) plus the weight of the mapped vertices vG->aH.
	wType getWeight(const edge eG, const edge eH) const { return w(m_BC_G->original(eG),m_BC_H->original(eH)); } // only for MCIS2
	unsigned getNumberOfUnfinishedCliqueComputations() const { return  m_unfinishedCliqueComputations; }
	// Construct BC-Trees, select initial B-Node, computes embedding into the plane
	BBP_MCSI(LabelFunction& labelFunction, InputGraph &iG, InputGraph &iH, bool enumerate=false, bool graphical_display=false, const vector<string>* simpleLabelToString=nullptr, clock_t maxCliqueTime=0, wType distancePenalty= WEIGHT_NOT_COMPATIBLE);

	~BBP_MCSI();

	BCTree& getBCG() const { return *m_BC_G; }
	BCTree& getBCH() const { return *m_BC_H; }
	wType computeSize(); // computes the weight of a BBPMCSI between G and H
	wType getSize() { return m_weightBBPIso; } // returns the weight of a BBPMCSI between G and H; use only after it has been computed!
	void computeIsomorphism(); // computes an isomorphism of weight weightBBPIso
	InputGraph* get_ig_G() const { return m_ig_G; }
	InputGraph* get_ig_H() const { return m_ig_H; }
};

#endif // BBPMCSI_H
