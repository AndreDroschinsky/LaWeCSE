/*
 * graphVisualization.cpp
 *
 *  Created on: 08.07.2015
 *      Author: droschin
 */


#include "graphVisualization.h"
#ifdef GRAPHICS
#include <ogdf/decomposition/BCTree.h>
//#include <string>
//#include <SFML/System/String.hpp>
#include <iostream>


void GraphVisualization::init(Graph &originalGraph, BCTree &bcTree, const sf::Color &vertexColor, const NodeArray<labelType> *nodeLabel, const EdgeArray<labelType> *edgeLabel,
		const vector<string>* simpleLabelToString, bool simpleLabel, int fogOffset)
{
	fogNodeEdgeOffset=fogOffset;
	//fogNodeEdgeOffset=0; //todo: Remove this
	m_basicVertexColor=vertexColor;
	m_G=&originalGraph;
	m_BC_G=&bcTree;
	m_nodeLabel=nodeLabel;
	m_edgeLabel=edgeLabel;
	m_SimpleLabelToString=simpleLabelToString;
	m_SimpleLabel=simpleLabel;
	if (!m_font_ArialBold.loadFromFile("Arial_Bold.ttf"))
	{
		cerr << "Font Arial_Bold.ttf nicht gefunden!"<<endl;
		exit(1);
	}
	if (!m_font_Arial.loadFromFile("arial.ttf"))
	{
		cerr << "Font arial.ttf nicht gefunden!"<<endl;
		exit(1);
	}
	// compute layout
	m_GraphAttributes.init(*m_G,GraphAttributes::edgeType | GraphAttributes::nodeType | GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics);
	edge e;
	//Graph::HiddenEdgeSetHandle hes=m_G->newHiddenEdgeSet();
	Graph::HiddenEdgeSet hes(*m_G);
	if (edgeLabel!=nullptr)
		forall_edges(e, *m_G)
			if ((*m_edgeLabel)[e]==LABEL_NOT_COMPATIBLE) // edges which only purpose is to connect unconnected parts will not be shown
			{
				//m_G->hideEdge(hes,e);
				hes.hide(e);
			}
	m_layout.call(m_GraphAttributes);
	//m_G->restoreEdges(hes);
	hes.restore();

#ifndef LAWECSU_NDEBUG
	//GraphIO::writeGML(graphAttributes,"GraphGA.gml");
 	//GraphIO::writeGML(graphAttributes,cout);
#endif

	// compute borders of representation in graphAttributes
	DRect bBox=m_GraphAttributes.boundingBox();
	m_vertexRadius=m_GraphAttributes.width(m_G->firstNode())/2;
	m_sizex=bBox.width()+2*m_vertexRadius;
	m_sizey=bBox.height()+2*m_vertexRadius;
	if (m_sizex>m_sizey)
	{
		m_rotate=true;
		swap(m_sizex,m_sizey);
	}
	m_vertexShape.init(*m_G);
	m_MCSidText.init(*m_G);
	m_VertexIdText.init(*m_G);
	m_edgeShape.init(*m_G);
	m_mcs_id.init(*m_G,0);
	m_BlockIdText.init(*m_G);

	// basic properties of graphical shapes
	node n;
	forall_nodes(n,*m_G)
	{
		m_vertexShape[n].setOutlineThickness(1);
		m_MCSidText[n].setFont(m_font_ArialBold);
		m_MCSidText[n].setFillColor(sf::Color::White);
		m_VertexIdText[n].setFont(m_font_Arial);
		m_VertexIdText[n].setFillColor(sf::Color(128,128,128,255));
	}
	forall_edges(e,*m_G)
	{
		m_edgeShape[e].setPrimitiveType(sf::LinesStrip); // Quads
		m_BlockIdText[e].setFont(m_font_Arial);
		m_BlockIdText[e].setStyle(sf::Text::Italic);
	}
}

void GraphVisualization::computeTextureAndCoords(double scaling)
{
	String str;
//	string str;
	m_scaling=scaling;
	DRect bBox=m_GraphAttributes.boundingBox();
	int left=bBox.p1().m_x;
	int top=bBox.p1().m_y;
	// transform vertices into graphical shapes
	node n; edge e;
	int offsetx=-left;
	int offsety=-top;
	forall_nodes(n,*m_G)
	{
		m_vertexShape[n].setRadius(m_vertexRadius*m_scaling);
		if (m_rotate)
			m_vertexShape[n].setPosition((m_GraphAttributes.y(n)+offsety)*m_scaling,(m_GraphAttributes.x(n)+offsetx)*m_scaling);
		else
			m_vertexShape[n].setPosition((m_GraphAttributes.x(n)+offsetx)*m_scaling,(m_GraphAttributes.y(n)+offsety)*m_scaling);
	}

	// add MCS and real vertex numbers
	forall_nodes(n,*m_G)
	{
		m_MCSidText[n].setCharacterSize(1.3*m_vertexRadius*m_scaling);
		m_VertexIdText[n].setCharacterSize(1.2*m_vertexRadius*m_scaling); // *2 weg
		str="";
		if (m_showVertexId)
		{
			str=to_string(n->index()+fogNodeEdgeOffset);
			if (m_showLabels && m_nodeLabel!=nullptr)
				str += " ";
		}
		if (m_showLabels && m_nodeLabel!=nullptr)
		{
			if (m_SimpleLabel)
				str += to_string((*m_nodeLabel)[n]);
			else
				str += (*m_SimpleLabelToString)[(*m_nodeLabel)[n]];
		}
		m_VertexIdText[n].setString(str);
	}


	if (m_rotate)
	{
		offsety+=2;
		offsetx+=1;
		forall_nodes(n,*m_G)
		{
			m_MCSidText[n].setPosition((m_GraphAttributes.y(n)+offsety)*m_scaling,(m_GraphAttributes.x(n)+offsetx)*m_scaling);
			if (m_mcs_id[n]<10 && m_mcs_id[n]>0)
				m_MCSidText[n].move(m_vertexRadius/2.1*m_scaling,0);
			m_VertexIdText[n].setPosition((m_GraphAttributes.y(n)+offsety+m_vertexRadius*1.2)*m_scaling,(m_GraphAttributes.x(n)+offsetx-m_vertexRadius*1.2)*m_scaling);
		}
	}
	else
	{
		offsetx+=2;
		offsety+=1;
		forall_nodes(n,*m_G)
		{
			m_MCSidText[n].setPosition((m_GraphAttributes.x(n)+offsetx)*m_scaling,(m_GraphAttributes.y(n)+offsety)*m_scaling);
			if (m_mcs_id[n]<10 && m_mcs_id[n]>0)
				m_MCSidText[n].move(m_vertexRadius/2.1*m_scaling,0);
			m_VertexIdText[n].setPosition((m_GraphAttributes.x(n)+offsetx+m_vertexRadius*1.2)*m_scaling,(m_GraphAttributes.y(n)+offsety-m_vertexRadius*1.2)*m_scaling);
		}
	}

	// transform edges into graphical shapes and add block and/or block IDs
	offsetx=-left+m_vertexRadius;
	offsety=-top+m_vertexRadius;
	float x,y,prevx=0,prevy=0;

	forall_edges(e,*m_G)
	{
		m_BlockIdText[e].setCharacterSize(1.2*m_vertexRadius*m_scaling);
		m_edgeShape[e].clear();
		float textBendPoint=0.5f+m_GraphAttributes.bends(e).size()/2.0f, curBendPoint=0;
		for (DPoint& dPoint:m_GraphAttributes.bends(e))
		{
			++curBendPoint;

			if (m_rotate)
			{
				x=(dPoint.m_y+offsety)*m_scaling;
				y=(dPoint.m_x+offsetx)*m_scaling;
			}
			else
			{
				x=(dPoint.m_x+offsetx)*m_scaling;
				y=(dPoint.m_y+offsety)*m_scaling;
			}
			m_edgeShape[e].append(Vertex(Vector2f(x,y)));
			//if (m_rotate)
			//	m_edgeShape[e].append(Vertex(Vector2f((dPoint.m_y+offsety)*m_scaling,(dPoint.m_x+offsetx)*m_scaling)));
			//else
			//	m_edgeShape[e].append(Vertex(Vector2f((dPoint.m_x+offsetx)*m_scaling,(dPoint.m_y+offsety)*m_scaling)));
			if (curBendPoint==textBendPoint) // display block ID it the center bend point
				m_BlockIdText[e].setPosition(x,y-m_vertexRadius*m_scaling*1.25);
			else if (curBendPoint==textBendPoint+0.5f)// display block ID in the middle of the center edge
				m_BlockIdText[e].setPosition((x+prevx)/2,(y+prevy)/2-m_vertexRadius*m_scaling*1.25);
			prevx=x;
			prevy=y;
			str="";
			if (m_showEdgeId)
			{
				str = to_string(e->index()+fogNodeEdgeOffset);
				if (m_showBlockId || (m_showLabels && m_nodeLabel!= nullptr))
					str +=" ";
			}
			if (m_showLabels && m_edgeLabel != nullptr)
			{
				if (m_SimpleLabel)
					str+=to_string((*m_edgeLabel)[e]);
				else
					str+=(*m_SimpleLabelToString)[(*m_edgeLabel)[e]];
				if (m_showBlockId)
					str += " ";
			}
			if (m_showBlockId)
				str += to_string(m_BC_G->bcproper(e)->index());
			m_BlockIdText[e].setString(str);
		}
	}
	m_texture.setTextureRect(IntRect(0, m_sizey*m_scaling, m_sizex*m_scaling, -m_sizey*m_scaling));

	// create renderTexture
	m_RenderTexture.create(m_sizex*m_scaling,m_sizey*m_scaling);
	m_texture.setTexture(m_RenderTexture.getTexture());
}

void GraphVisualization::setVertexMcsID(node n, int id)
{
//	bool newid=false;
	//if (m_mcs_id[n] != id)
		//newid=true;
	m_mcs_id[n]=id;
	max_mcs_id=max(max_mcs_id,id);
	String label=to_string(id);
	//if (newid==false)
		//label="";
	if (id<10)
		m_MCSidText[n].move(m_vertexRadius/2.1*m_scaling,0);
	m_MCSidText[n].setString(label);
}

void GraphVisualization::applyColors()
{
	node n;
	forall_nodes(n,*m_G)
	{
		if (m_mcs_id[n]==0) // not mapped
		{
			m_vertexShape[n].setFillColor(m_basicVertexColor);
			m_MCSidText[n].setString("");
			m_vertexShape[n].setOutlineColor(sf::Color(210,210,210,255));
		}
		else
		{
			int color=765*(m_mcs_id[n]-1)/max_mcs_id;
			Uint8 red=color<256?255-color:(color>510?color-510:0);
			Uint8 green=color<256?color:(color<510?510-color:0);
			Uint8 blue=color>510?765-color:(color>255?color-255:0);
			green*=0.6;
			red*=0.8;
			m_vertexShape[n].setFillColor(sf::Color(red,green,blue,255));
			m_vertexShape[n].setOutlineColor(sf::Color::Black);
		}
	}
	edge e;
	forall_edges(e,*m_G)
	{
		sf::Color edgecolor;
		if (m_BlockColors)
		{
			int color=2*765*(m_BC_G->bcproper(e)->index())/(m_BC_G->numberOfBComps()+m_BC_G->numberOfCComps());
			color+=(color/765)*(765/2/(m_BC_G->numberOfBComps()+m_BC_G->numberOfCComps()));

			Uint8 red=color<256?255-color:(color>510?color-510:0);
			Uint8 green=color<256?color:(color<510?510-color:0);
			Uint8 blue=color>510?765-color:(color>255?color-255:0);
			Uint8 alpha=255;
			green*=0.6;
			red*=0.8;
			if (!(m_mcs_id[e->target()]!=0 && m_mcs_id[e->source()]!=0))
				alpha=100;
			edgecolor=sf::Color(red,green,blue,alpha);
			m_BlockIdText[e].setFillColor(sf::Color(red,green,blue,alpha));
		}
		else
		{
			if (m_mcs_id[e->target()]!=0 && m_mcs_id[e->source()]!=0)
				edgecolor=sf::Color(0,0,0,255);
			else
				edgecolor=sf::Color(210,210,210,255);
			m_BlockIdText[e].setFillColor(sf::Color(180,180,180,255));
		}

		for (unsigned i=0; i<m_edgeShape[e].getVertexCount(); ++i)
			m_edgeShape[e][i].color=edgecolor;
	}
}

void GraphVisualization::drawAll()
{
	m_RenderTexture.clear(sf::Color(240,240,240));
 	edge e;
	forall_edges(e,*m_G)
 		if ((*m_edgeLabel)[e]!=LABEL_NOT_COMPATIBLE)
 			m_RenderTexture.draw(m_edgeShape[e]);

	node n;
	forall_nodes(n,*m_G)
	{
		m_RenderTexture.draw(m_vertexShape[n]);
		m_RenderTexture.draw(m_MCSidText[n]);
		if (m_showVertexId || m_showLabels)
			m_RenderTexture.draw(m_VertexIdText[n]);
	}

	if (m_showBlockId || m_showEdgeId || m_showLabels)
 		forall_edges(e,*m_G)
			if ((*m_edgeLabel)[e]!=LABEL_NOT_COMPATIBLE)
				m_RenderTexture.draw(m_BlockIdText[e]);
}

#endif /* GRAPHICS */
