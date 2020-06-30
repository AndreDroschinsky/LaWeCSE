#ifndef MCIS2_H
#define MCIS2_H

#include <global.h>

#pragma warning(push, 0)
#include <ogdf/decomposition/BCTree.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>
#include <ogdf/basic/Graph.h>
#pragma warning(pop)
#include <vector>
#include <forward_list>
#include <queue>
#include <stack>

class BBP_MCSI;

using namespace ogdf;
using namespace std;

class NodeBijection
{
protected:
    NodeArray<node> m_map;
    NodeArray<node> m_invMap;
    std::vector<node> m_mapped;

public:
    NodeBijection (const Graph& G, const Graph& H):m_map(G,NULL), m_invMap(H,NULL) { };
    ~NodeBijection() {};

    node map(node v) { return m_map[v]; };
    node invMap(node v) { return m_invMap[v]; };
    void addMap(node u, node v) { m_map[u] = v; m_invMap[v] = u; m_mapped.push_back(u); };

    void clear() {
        for (node u:m_mapped)
        {
            m_invMap[m_map[u]]=NULL;
            m_map[u]=NULL;
        }
        m_mapped.clear();
        //m_map.fill(NULL); m_invMap.fill(NULL);
    };

};

struct ProductnodeNoExpand {
    int nodeG; int nodeH;
    ProductnodeNoExpand(int nG, int nH) : nodeG(nG),nodeH(nH) {}
    ProductnodeNoExpand() : nodeG(0),nodeH(0) {}
    friend ostream& operator<<(ostream &os, ProductnodeNoExpand &pn)
    {
        os << "("<<pn.nodeG<<","<<pn.nodeH<<")";
//		if (pn.expand)
//			os << "E";
        os << " ";
        return os;
    }

};

struct Productnode:ProductnodeNoExpand {
    bool expand=false;
    Productnode(int nG, int nH):ProductnodeNoExpand(nG,nH) {};
    Productnode():ProductnodeNoExpand() {};
    friend ostream& operator<<(ostream &os, Productnode &pn)
    {
        os << "("<<pn.nodeG<<","<<pn.nodeH<<")";
        if (pn.expand)
            os << "E";
        os << " ";
        return os;
    }
};


#define ProductNodeNeighborsMAXNEIGHBORS 4 // maximum degree in molecular graphs
class ProductNodeNeighbors {
    // +0 neighbor of nG     +1 his quantity      +2 neighbor of nH     +3 his quantity
    vector<int> m_neighborsAndQuantity;
    size_t m_curNeighborsG = 0, m_curNeighborsH = 0;
    bool m_neighborship_overfull=false; // if into neighborship more data than possible needs to be stored, nodes will not be sorted out during clique computations
    bool m_wasRemoved=true; // if true, this vertex has been removed; necessary for recursive removal and hasSufficientAmountOfNeighbors; true initially, because nothing was added

public:
    // initialize vector to store neighborship
    ProductNodeNeighbors();

    // adds single target product node to neighborship
    void addNeighbor(ProductnodeNoExpand &target);

    //removes single target product node from neighborship; return true iff this product node had sufficient neighbors before, but not any more
    bool removeNeighbor(ProductnodeNoExpand &target);

    void setwasRemoved(bool wasRemoved) { m_wasRemoved = wasRemoved; }
    bool getwasRemoved() { return m_wasRemoved; }

    bool hasSufficientAmountOfNeighbors();
};

enum { NO_EDGE, C_EDGE, D_EDGE };

class MaxComPart2Tree
{
protected:
    const Graph &m_G, &m_H;
    StaticPlanarSPQRTree *m_SPG, *m_SPH;
    EdgeArray<EdgeArray<wType*> > m_D;
    EdgeArray<EdgeArray<char> > m_detectCaseCache;
    BBP_MCSI *m_bbp=nullptr;
    NodeBijection	*m_map_snode;
    bool m_enumerate=false;

    Graph m_MCS;
    NodeArray<node> m_vGtovH, m_vGtoMCS, m_vMCStovG;
    EdgeArray<edge> m_eGtoeH,  m_eMCStoeG;
    const NodeArray<node> &m_vGtoaG, &m_vHtoaH;
    const EdgeArray<edge> &m_eGtoaeG, &m_eHtoaeH;
    NodeArray<wType> m_vMCStoWeight;
    EdgeArray<wType> m_eMCStoWeight;
    NodeArray<bool> m_vMCStoExpand; // Possible expansion. Set during applyWeightsToMCS().
    NodeArray<std::vector<node>> m_VbcMCS; // vertices of the blocks and bridges of the MCS

    void applyWeightsToMCS(); // computes the weights for each edge and node of m_MCS and defaults WEIGHT_NOT_COMPATIBLE to the edges of D
    void computeSplitIsomorphisms(std::forward_list<std::vector<Mapping>>& maximumMappingList, wType &maxWeight); // computes and stores the maximum split isomorphisms and returns maximum weight
    void deleteCurrentSolution(); // deletes the stored current solution
    node m_xG,m_vG; // excluded vertex or given mapping

    void addNodeMappingToMCS(node vG, node vH); // add a node mapping to the MCS
    void addEdgeMappingToMCS(edge eG, edge eH); // add a edge mapping to the MCS

    node getCopy	(const node v, const Skeleton &S);
    void matchCycle	(const node uG, const node vG, const node uH, const node vH, const node sG, const node sH);
    void matchEdge	(const node uG, const node vG, const node uH, const node vH, const node sG, const node sH, const edge eG, const edge eH);
    edge getRealEdge(const edge sNodeSkeletonEdge, const Skeleton &skel, const StaticSPQRTree *sprqtree);
    int detectCase(edge eG, edge eH);

    node original	(const Skeleton &skel, const node u){return u == NULL ? NULL : skel.original(u);};
    enum class GraphSelect  {GraphG, GraphH};
    const Graph&	getOriginalGraph	(GraphSelect g){return g == GraphSelect::GraphG ? m_G : m_H;};
    StaticSPQRTree*	getSPGraph			(GraphSelect g){return g == GraphSelect::GraphG ? m_SPG : m_SPH;};


    // additional part for non outerplanar 2MCS computation
    bool m_outerplanar; // set to true / false during initialization
    int m_nodesProductGraph=0;
    int m_sizeG=0, m_sizeH=0;
    //std::vector<bool> *adjG, *adjH;
    vector<edge> *m_adjG, *m_adjH; // triangular adjacency matrix (from lower IDs to higher IDs)
    node *m_nodesG, *m_nodesH;
    wType *m_maxWeightNonOuterPtr=nullptr;
    forward_list<vector<Mapping>>* m_maximumMappingListPtr=nullptr;
    NodeArray<wType> *m_maxWeightArrayNonOuterPtr=nullptr;
    NodeArray<forward_list<vector<Mapping>>>* m_maximumMappingArrayListPtr=nullptr;
    vector<ProductnodeNoExpand> m_R;
    clock_t m_cliqueTimeStart=0;
    clock_t m_maxCliqueTime;
    wType *m_productNodeWeight=nullptr; // store weight between node pairs
    vector<bool> m_productNodeExpand; // store possible expansion through product nodes
    ProductNodeNeighbors *m_productNodeNeighbors=nullptr; // neighbors and quantity

    void c_Clique(vector<ProductnodeNoExpand> &P, vector<ProductnodeNoExpand> &Q, bool biconnected=false, wType curWeight=0);
    void ReportClique(wType cliqueWeight);
    int getEdgeType(ProductnodeNoExpand from, ProductnodeNoExpand to); // edge type in product graph (enum { NO_EDGE, C_EDGE, D_EDGE };)
    wType productNodeWeight(ProductnodeNoExpand &n); // return weight of a vertex in the product graph
    wType productEdgeWeight(ProductnodeNoExpand &from, ProductnodeNoExpand &to); // return weight of and edge in the product graph
    void eliminateNodesOfDegreeOneAndZeroInProductGraph(); // determine vertices in product graph, which cannot be connected through two C_EDGEs
    void addNeighborship(ProductnodeNoExpand& source); // Add vertex source into the neighborship of all of sources neighbors
    // Recursively remove vertex source from the neighborship of all of sources neighbors
    // Puts removed vertices on the stack sInsufficientNeighbors if not nullptr
    void removeNeighborship(ProductnodeNoExpand source, stack<ProductnodeNoExpand>* sInsufficientNeighbors_copy=nullptr);
    // additional part for checking biconnectivity
    int *m_bi_num=nullptr, *m_bi_low=nullptr;
    int m_bi_count=0;
    bool hasAC(int u, int parent);

unsigned m_cliqueRecursions=0;

public:
    MaxComPart2Tree	(InputGraph& ig_G, InputGraph& ig_H, node bG, node bH, BBP_MCSI* bbp, node xG, node vG, clock_t maxCliqueTime)
        :m_G(ig_G.m_bGAG[bG]), m_H(ig_H.m_bGAG[bH]), m_SPG(ig_G.m_bGToSPQRTree[bG]), m_SPH(ig_H.m_bGToSPQRTree[bH]), m_vGtoaG(ig_G.m_SPaGtoaG[bG]),
         m_vHtoaH(ig_H.m_SPaGtoaG[bH]), m_eGtoaeG(ig_G.m_SPeGtoaeG[bG]), m_eHtoaeH(ig_H.m_SPeGtoaeG[bH]), m_xG(xG), m_vG(vG),
         m_outerplanar((ig_G.m_bGOuterplanar[bG] && ig_H.m_bGOuterplanar[bH])?true:false), m_maxCliqueTime(maxCliqueTime)  { m_bbp=bbp; init(); };

    ~MaxComPart2Tree();

    void init();
    unsigned computeBlock(std::forward_list<std::vector<Mapping>>& maximumMappingList, wType &maxWeight, const bool &enumerate); // for mappings with possible excluded vertex
    unsigned computeFixedMapping(NodeArray<std::forward_list<std::vector<Mapping>>> &aGmaximumMappingList, NodeArray<wType> &maxWeight, const bool &enumerate); // for mappings with given mapping from m_vG
};

#endif //MCIS2_H
