#ifndef _GLOBAL_
#define _GLOBAL_
#define DNDEBUG
#ifdef DEBUG
#define DNDEBUG
#else
//#define GRAPHICS
#endif


#define LAWECSU_NDEBUG

#pragma warning(push, 0)
#include <ogdf/basic/Graph_d.h>
#include <ogdf/decomposition/BCTree.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>
#pragma warning(pop)
#include <unordered_map>

using namespace ogdf;


// OGDF Makros from 2015-06-29

#define forall_edges(e,G) for((e)=(G).firstEdge(); (e); (e)=(e)->succ())
#define forall_nodes(v,G) for((v)=(G).firstNode(); (v); (v)=(v)->succ())
#define forall_adj(adj,v) for((adj)=(v)->firstAdj(); (adj); (adj)=(adj)->succ())

typedef double wType;
typedef int16_t labelType;
typedef int32_t labelPairType;
#define LABEL_MULTIPLYER 65536

#define EQDIST 0.0000001
#define WEIGHT_UNDEF -2147483647
#define WEIGHT_NOT_COMPATIBLE -2147483648
#define LABEL_NOT_COMPATIBLE -32767

#ifdef MAIN
  #define EXTERN
#else
  #define EXTERN extern
#endif

//EXTERN int cntnb;

// Computed local maximum part of an isomorphism
struct Mapping
{
	ogdf::node vG;
	ogdf::node aH;
	bool expand;
	Mapping(ogdf::node vG,ogdf::node aH, bool expand):vG(vG),aH(aH),expand(expand) {}
	Mapping():vG(nullptr),aH(nullptr),expand(false) {}
};

struct BGxGPair
{
	ogdf::node bG;
	ogdf::node xG;
	BGxGPair(ogdf::node bG, ogdf::node xG):bG(bG),xG(xG) {}
};

struct ComputationTime
{
	unsigned mintime;
	unsigned maxtime;
	unsigned totaltime;
	unsigned instances;
	ComputationTime(unsigned time):mintime(time),maxtime(time),totaltime(time),instances(1) {}
	unsigned average() { return totaltime/instances; }
	void addTime(unsigned time) { mintime=min(time,mintime); maxtime=max(time,maxtime); totaltime+=time; ++instances; }
};

struct LabelFunction
{
	wType sameNodeLabel=1.0;
	wType differentNodeLabel=WEIGHT_NOT_COMPATIBLE;
	wType sameEdgeLabel=1.0;
	wType differentEdgeLabel=WEIGHT_NOT_COMPATIBLE;
	std::unordered_map<labelPairType,wType> weightTable;
};

struct InputGraph
{
	Graph graph;
	string graphLabel;
	ogdf::BCTree *bcTree=nullptr;

	// Mappings between auxiliary graph in the BC tree and the single parts (graphs)
	NodeArray<std::vector<node>> m_VbcG; // vertices of the auxiliary graph in a B oder C-node of BCG
	NodeArray<Graph> m_bGAG; // the single parts of the auxiliary graph of G
	NodeArray<node> m_aGtoSPaG; // mapping from aG to the nodes of the single parts of AG
	NodeArray<NodeArray<node>> m_SPaGtoaG; // mapping from the nodes of the single parts of AG to aG
	NodeArray<EdgeArray<edge>> m_SPeGtoaeG; // mapping from the edges of the single parts of AG to the edges aG
	NodeArray<bool> m_bGOuterplanar; // true, if block is outerplanar
	bool mOuterplanar=true; // true, if whole graph is outerplanar

	NodeArray<labelType>* nodeLabel=nullptr;
	EdgeArray<labelType>* edgeLabel=nullptr;
	bool simpleLabel=true;
	unsigned graphIndexToDBIndex=0;

	NodeArray<StaticPlanarSPQRTree*> m_bGToSPQRTree; // map from block to SPRQ tree
	wType size=0; // weight of mapping of vertices to themselves
	int fogNodeEdgeOffset=0; // Fog format starts vertices and edges at id 1

	~InputGraph()
	{
		if (nodeLabel != nullptr)
			delete nodeLabel;
		if (edgeLabel != nullptr)
			delete edgeLabel;

		if (bcTree != nullptr)
		{
			//node b;
			//forall_nodes(b,bcTree->bcTree())
			for(node b:bcTree->bcTree().nodes)
			{
				//if (m_VbcG[b].size()>2)
				if (m_bGToSPQRTree[b] != nullptr)
				{
					delete m_bGToSPQRTree[b];
				}
			}
			delete bcTree;
		}
	}
};

#endif
