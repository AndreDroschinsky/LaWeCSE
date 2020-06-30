#ifndef TREEGEN_H
#define TREEGEN_H

#pragma warning(push, 0)
#include <ogdf/basic/Graph.h>
#pragma warning(pop)
#include <global.h>

using namespace ogdf;

class OuterPlanarGen
{
private:
	int		m_maxComponentSize;
	int		m_maxBridgeLength;
	bool	m_biConnected;
public:
	OuterPlanarGen(): m_maxComponentSize(-1), m_maxBridgeLength(-1), m_biConnected(false){};
	OuterPlanarGen(int CompSize, int BridgeLength, bool biConnected): m_maxComponentSize(CompSize), m_maxBridgeLength(BridgeLength), m_biConnected(biConnected){};

	void generate(const int size, Graph& Graph);
};

#endif // TREEGEN_H
