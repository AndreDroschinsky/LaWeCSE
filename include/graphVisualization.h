/*
 * graphVisualization.h
 *
 *  Created on: 08.07.2015
 *      Author: droschin
 */

#ifndef GRAPHVISUALIZATION_H_
#define GRAPHVISUALIZATION_H_

using namespace std;
#include <global.h>

#ifdef GRAPHICS
#include <ogdf/decomposition/BCTree.h>
//#include <ogdf/planarity/PlanarizationGridLayout.h> // NOT working in ogdf 2018
#include <ogdf/planarity/PlanarizationLayout.h>

#include <SFML/Graphics.hpp>

using namespace sf;
using namespace ogdf;

class GraphVisualization {
private:
	Graph* m_G=nullptr; // original graph
	const BCTree* m_BC_G=nullptr; // BC-Tree
	GraphAttributes m_GraphAttributes;
 	//PlanarizationGridLayout m_layout;
 	PlanarizationLayout m_layout;
 	NodeArray<sf::CircleShape> m_vertexShape;
	Font m_font_ArialBold,m_font_Arial;
	NodeArray<sf::Text> m_MCSidText;
	NodeArray<sf::Text> m_VertexIdText;
	bool m_showVertexId=false;
	EdgeArray<sf::Text> m_BlockIdText;
	bool m_showBlockId=false;
	bool m_BlockColors=true;
	bool m_showEdgeId=false;
	bool m_showLabels=true;
	NodeArray<int> m_mcs_id; //
	int max_mcs_id=0;
	sf::Color m_basicVertexColor;
	EdgeArray<sf::VertexArray> m_edgeShape;
	RenderTexture m_RenderTexture;
	Sprite m_texture;
	int m_sizex=0,m_sizey=0;
	bool m_rotate=false;
	double m_vertexRadius=0; // vertex radius
	double m_scaling=1.0;
	const NodeArray<labelType> *m_nodeLabel=nullptr;
	const EdgeArray<labelType> *m_edgeLabel=nullptr;
	const vector<string>* m_SimpleLabelToString=nullptr;
	bool m_SimpleLabel=true;
	int fogNodeEdgeOffset=0;

public:
	GraphVisualization() {};
	~GraphVisualization() {};
	void init(Graph &originalGraph, BCTree &bcTree, const sf::Color &vertexColor, const NodeArray<labelType> *nodeLabel=nullptr, const EdgeArray<labelType> *edgeLabel=nullptr,
			const vector<string>* simpleLabelToString=nullptr, bool simpleLabel=true, int fogOffset=0);
	int getSizex() const { return m_sizex; };
	int getSizey() const { return m_sizey; };
	void setTexturePosition(int x, int y) { m_texture.setPosition(x,y); };
	const Sprite& getTexture() { return m_texture; };
	void computeTextureAndCoords(double scaling);
	void toogleShowVertices() { m_showVertexId=!m_showVertexId; }
	void toogleShowEdges() { m_showEdgeId=!m_showEdgeId; }
	void toogleShowBlocks() { m_showBlockId=!m_showBlockId; }
	void toogleBlockColors() { m_BlockColors=!m_BlockColors; }
	void toogleLabels() { m_showLabels=!m_showLabels; }
	void setVertexMcsID(node n, int id);
	void applyColors(); // applies the computed colors of the isomorphism to the shapes
	void drawAll(); // Draws the complete graph onto a texture
};

#endif /* GRAPHICS */

#endif /* GRAPHVISUALIZATION_H_ */


