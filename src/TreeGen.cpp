#include <TreeGen.h>

/*void OuterPlanarGen::generate(const int size, const float density, ogdf::Graph& G)
{
	this->generate(size, G);
}*/

void OuterPlanarGen::generate(const int size, ogdf::Graph& G)
{
	// Initialize random seed:
//	time_t time_seed;
//	time(&time_seed);
//    srand((unsigned int)time_seed);

	int freeVert = size;
	int bridgeLen = m_maxBridgeLength;
	if(bridgeLen == -1)
		bridgeLen = freeVert;

	int compSize = m_maxComponentSize;
	if(compSize == -1)
		compSize = freeVert;
	
	ogdf::List<ogdf::edge> edgeCandidates;

	G.clear();

	while(freeVert > 0)
	{
		// Decide whether to build a bridge or a component:
		// Build a component, i.e. a ring which can either be attached to the end of a bridge or 
		// connected to another component with one common edge:
		if(rand() % 2 || m_biConnected)
		{
			// Minimum size of a component is 3:
			int i = (rand() % compSize) + 1;

			// Nothing generated yet: 
			ogdf::node start;
			if(edgeCandidates.empty() && G.empty())
			{
				start = G.newNode();
				ogdf::node cur = start;

				for(int j = 1; j < i+2; j++)
				{
					ogdf::node next = G.newNode(); freeVert--;
					edgeCandidates.pushBack(G.newEdge(cur, next));

					cur = next;
				}
				
				edgeCandidates.pushBack(G.newEdge(cur, start));
			}// if(edgeCandidates.empty() && G.empty())
			else
			{
				// Decide whether to choose an edge or a node as root:
				int r = rand() % edgeCandidates.size() + G.numberOfNodes();

				if(r < G.numberOfNodes())
				{
					start = G.chooseNode();
					ogdf::node cur = start;

					for(int j = 1; j <= i; j++)
					{
						ogdf::node next = G.newNode(); freeVert--;
						edgeCandidates.pushBack(G.newEdge(cur, next));

						cur = next;
					}
				
					edgeCandidates.pushBack(G.newEdge(cur, start));
				} //if(r < G.numberOfNodes())
				else
				{
					int index = r - G.numberOfNodes();
					start = (*edgeCandidates.get(index))->source();
					ogdf::node cur = start;

					for(int j = 1; j <= i; j++)
					{
						ogdf::node next = G.newNode(); freeVert--;
						edgeCandidates.pushBack(G.newEdge(cur, next));

						cur = next;
					}
				
					edgeCandidates.pushBack(G.newEdge(cur, (*edgeCandidates.get(index))->target()));
					edgeCandidates.del(edgeCandidates.get(index));
				} // else (r < G.numberOfNodes())
			} // else(edgeCandidates.empty() && G.empty())
		} // if(freeVert > 3 && rand() % 2)
		else
		{
			int i = (rand() % bridgeLen) + 1;

			ogdf::node start;
			
			if(edgeCandidates.empty() && G.empty())
			{
				start = G.newNode(); freeVert--;
			}
			else
				start = G.chooseNode();

			ogdf::node cur = start;
			for(int j = 1; j <= i; j++)
			{
				ogdf::node next = G.newNode(); freeVert--;
				edgeCandidates.pushBack(G.newEdge(cur, next));

				cur = next;
			}	
		} //else(freeVert > 3 && rand() % 2)
	} // while
}
