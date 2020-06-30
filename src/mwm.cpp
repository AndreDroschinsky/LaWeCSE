/*
 * mwm.cpp
 *
 *  Created on: 10.06.2015
 *      Author: droschin
 */

#include <mwm.h>
#include <hungarian.h>
#include <stack>
#include <algorithm>
#pragma warning(push, 0)
#include <ogdf/basic/simple_graph_alg.h>
#pragma warning(pop)

MWM::MWM(bool withEnumeration):m_nodesA(0),m_nodesJ(0),m_withEnumeration(withEnumeration)
{
    m_egdeweight.init(m_G);
    m_copy.init(m_G);
    m_nodeType.init(m_G);
    m_vMatchingTovOriginal.init(m_G);
    m_mate.init(m_G);
    m_dual.init(m_G);
}

MWM::~MWM()
{
    // Free reserved memory of first trim-call
    if (m_withEnumeration)
    {
//#ifndef N_MATCHING_DEBUG
        if (m_EnumDeg2ReplacedPath.size() != m_EnumRestore.size())
            cout << "size diff";
//#endif
        for (unsigned i=0; i<m_EnumRestore.size(); ++i)
        {
            while (!m_EnumRestore[i].empty())
            {
                for (edge e:m_EnumRestore[i].back()->shortcutedges)
                {
                    delete m_EnumDeg2ReplacedPath[i][e];
                    //m_EnumDeg2ReplacedPath[e]=nullptr;
                }
                delete m_EnumRestore[i].back();
                m_EnumRestore[i].pop_back();
            }
        }
    }

    if (m_EnumStoreMatchings)
    {
        for (unsigned i=0; i<m_EnumStoredMatchings.size(); ++i)
        {
            if (m_EnumStoredMatchings[i])
            {
                vector<vector<edge>*>& storedMatching=m_EnumStoredMatching[i];
                for (vector<edge>* ptr:storedMatching)
                {
                    delete ptr;
                }
            }
        }
    }
}

void MWM::addAgent(node agent)
{
    const node newAgent=m_G.newNode();
    m_vMatchingTovOriginal[newAgent]=agent;
    m_agent.push_back(newAgent);
    m_nodeType[newAgent]=NodeType::AGENT;
    ++m_nodesA;
#ifndef N_MATCHING_DEBUG
    cout << "A"<<m_nodesA<<":"<<agent<<"  ";
#endif
}

void MWM::addJob(node job)
{
    const node newJob=m_G.newNode();
    m_vMatchingTovOriginal[newJob]=job;
    m_job.push_back(newJob);
    m_nodeType[newJob]=NodeType::JOB;
    ++m_nodesJ;
#ifndef N_MATCHING_DEBUG
    cout << "J"<<m_nodesJ<<":"<<job<<"  ";
#endif
}

void MWM::addGraphCopy()
{
    if (m_graph_copied) // only copy once
        return;
    /*bool regularCompute=false;
    if (!m_withEnumeration) // fast matching computation only without enumeration
    {
        NodeArray<int> component(m_G);
        int numComponents = connectedComponents(m_G,component);
        //cout << "numc:"<<numComponents<<endl;
        m_agentsInComp.resize(numComponents);
        m_jobsInComp.resize(numComponents);
        //cout << "A: "<<m_nodesA<<" J:"<<m_nodesJ<<"  ";
        for (node n:m_G.nodes)
        {
            if (m_nodeType[n]==NodeType::AGENT)
            {
                if (++m_agentsInComp[component[n]]>1)
                    regularCompute=true;
            }
            else
            {
                ++m_jobsInComp[component[n]];
            }
        }
        //for (int i=0; i<numComponents; ++i)
            //cout <<m_agentsInComp[i]<< " "<<m_jobsInComp[i]<< "  ";
        //cout << endl;
        if (!regularCompute)
            return;
    }*/

    // edges of the graph without auxiliary edges
    vector<edge> allEdges;
    for (edge e : m_G.edges)
        allEdges.push_back(e);

    // copy agents and connect them
    node cp;
    for (unsigned i=0; i<m_nodesA; ++i)
    {
        cp = m_G.newNode();
        m_job.push_back(cp);
        edge e = m_G.newEdge(m_agent[i],cp);
        m_egdeweight[e]=0;
        m_copy[m_agent[i]]=cp;
        m_copy[cp]=m_agent[i];
        m_nodeType[cp]=NodeType::AGENTCOPY;
#ifndef N_MATCHING_DEBUG
        cout << "AC"<<i+1<<":"<<m_vMatchingTovOriginal[m_agent[i]]<<"  ";
#endif
    }
    // copy jobs and connect them

    for (unsigned i=0; i<m_nodesJ; ++i)
    {
        cp=m_G.newNode();
        m_agent.push_back(cp);
        edge e = m_G.newEdge(cp,m_job[i]);
        m_egdeweight[e]=0;
        m_copy[m_job[i]]=cp;
        m_copy[cp]=m_job[i];
        m_nodeType[cp]=NodeType::JOBCOPY;
#ifndef N_MATCHING_DEBUG
        cout << "JC"<<i+1<<":"<<m_vMatchingTovOriginal[m_job[i]]<<"  ";
#endif
    }
#ifndef N_MATCHING_DEBUG
    cout <<endl;
#endif

    // copy edges of the graph
    for (edge e:allEdges)
    {
        edge e2=m_G.newEdge(m_copy[e->target()],m_copy[e->source()]);
        m_egdeweight[e2]=m_egdeweight[e];
    }
    m_graph_copied = true;
}

void MWM::addEdge(unsigned indexAgent, unsigned indexJob, wType weight_)
{
    edge e = m_G.newEdge(m_agent[indexAgent],m_job[indexJob]);
    m_egdeweight[e]=weight_;
    //if (!m_withEnumeration && m_nodesA==1)
        //return;
    if (m_graph_copied)
    {
        e = m_G.newEdge(m_copy[m_job[indexJob]],m_copy[m_agent[indexAgent]]);
        m_egdeweight[e]=weight_;
    }
}

void MWM::constructEqualitySubgraph(const node jobToRemove,unsigned nodesInEachSet)
{
    // nodesInEachSet: generally m_nodesA+m_nodesJ vertices. Similar graphs with one job removed have one less, unless there was no iteration of the Hungarian method
    //cout << "param: "<<nodesInEachSet;
    //nodesInEachSet=m_nodesA+m_nodesJ-min(m_EnumerationGraph.size(),1UL)+(jobToRemove==nullptr?0:1);
    //cout << " comp param: "<<nodesInEachSet<<endl;

    //unsigned nodesInEachSet=m_nodesA+m_nodesJ; //-min(m_EnumerationGraph.size(),1UL)+(jobToRemove==nullptr?0:1);
    m_EnumerationGraph.push_back(Graph());
    Graph &enumG=m_EnumerationGraph.back();
    m_vEnumTovOriginal.push_back(NodeArray<node>(enumG));
    NodeArray<node>& vEnumTovOriginal=m_vEnumTovOriginal.back();
    NodeArray<node> vEnumTomG(enumG);
    m_EnumMatched.push_back(EdgeArray<bool>(enumG,false));
    EdgeArray<bool>& enumMatched=m_EnumMatched.back();
    m_EnumNodeType.push_back(NodeArray<NodeType>(enumG));
    NodeArray<NodeType>& enumNodeType=m_EnumNodeType.back();
    m_EnumCopy.push_back(NodeArray<node>(enumG));
    NodeArray<node>& enumCopy=m_EnumCopy.back();
    m_EnumNonIsolatedVertices.push_back(vector<node>());
    vector<node>& enumNonIsolatedVertices=m_EnumNonIsolatedVertices.back();
    enumNonIsolatedVertices.reserve(nodesInEachSet*2);
    m_EnumRestore.push_back(vector<EnumRestore*>());
    m_EnumDeg2TrimmedEdges.push_back(EdgeArray<bool>(enumG,false));
    m_EnumDeg2ReplacedPath.push_back(EdgeArray<vector<edge>*>(enumG,nullptr));
    NodeArray<node> vm_G_to_enumG(m_G);
    m_EnumFixedMatching.push_back(vector<edge>());
    m_EnumFixedMatching.back().reserve(min(m_nodesA,m_nodesJ));
    m_EnumStack.push_back(stack<StackOP>({StackOP::ENUM}));
    m_EnumSCC.push_back(NodeArray<unsigned>(enumG));
    m_EnumMaxSCCNumber.push_back(0);
    m_EnumDFSstatus.push_back(NodeArray<short>(enumG));
    m_EnumPreviousEdgeOnCircle.push_back(EdgeArray<edge>(enumG));
    if (m_EnumStoreMatchings)
    {
        m_EnumStoredMatching.push_back(vector<vector<edge>*>());
        m_EnumStoredMatchings.push_back(false);
        m_EnumCurStoredMatching.push_back(0);
    }

#ifndef N_MATCHING_DEBUG
//	cout <<"Equality Subgraph:\nNodes: ";
#endif

    // add vertices to G
    node matchingNode,enumNode;
    for (unsigned nodeIndex=0; nodeIndex<nodesInEachSet; ++nodeIndex)
    {
        matchingNode=m_agent[nodeIndex];
        if (jobToRemove==nullptr || matchingNode != m_copy[jobToRemove])
        {
            enumNode=enumG.newNode(matchingNode->index());
            vm_G_to_enumG[matchingNode]=enumNode;
            vEnumTomG[enumNode]=matchingNode;
            vEnumTovOriginal[enumNode]=m_vMatchingTovOriginal[matchingNode];
            enumNodeType[enumNode]=m_nodeType[matchingNode];
            enumNonIsolatedVertices.push_back(enumNode);
#ifndef N_MATCHING_DEBUG
            //cout <<	enumNode->index()<< " ";
#endif
        }

        matchingNode=m_job[nodeIndex];
        if (matchingNode != jobToRemove)
        {
            enumNode=enumG.newNode(matchingNode->index());
            vm_G_to_enumG[matchingNode]=enumNode;
            vEnumTomG[enumNode]=matchingNode;
            vEnumTovOriginal[enumNode]=m_vMatchingTovOriginal[matchingNode];
            enumNodeType[enumNode]=m_nodeType[matchingNode];
            enumNonIsolatedVertices.push_back(enumNode);
#ifndef N_MATCHING_DEBUG
            //cout <<	enumNode->index()<< " ";
#endif
        }
    }
    // node copies
    forall_nodes(enumNode,enumG)
        enumCopy[enumNode]=vm_G_to_enumG[m_copy[vEnumTomG[enumNode]]];

#ifndef N_MATCHING_DEBUG
//	cout << endl;
#endif
    // add equality edges to G
    edge e,eEnum;
    forall_edges(e,m_G)
    {
        // both end points may not be the possibly excluded job
        if (jobToRemove==nullptr || (e->source()!=m_copy[jobToRemove] && e->target()!=jobToRemove))
        {
            // edge part of the equality subgraph?
            if (abs(m_dual[e->source()]+m_dual[e->target()]+m_egdeweight[e]) <= EQDIST)
            {
                // must also be true for copy of the edge
                if (abs(m_dual[m_copy[e->source()]]+m_dual[m_copy[e->target()]]+m_egdeweight[e]) <= EQDIST)
                {
                    // at least one node must be an original node
                    //if (e->source()->index()<int(m_nodesA+m_nodesJ) || e->target()->index()<int(m_nodesA+m_nodesJ))
                    if (m_nodeType[e->source()]<= NodeType::JOB || m_nodeType[e->target()]<= NodeType::JOB)
                    {
                        eEnum=enumG.newEdge(vm_G_to_enumG[e->source()],vm_G_to_enumG[e->target()]);
                        if (m_mate[e->source()]==e->target())
                        {
                            enumMatched[eEnum]=true;
    #ifndef N_MATCHING_DEBUG
        //				cout << "M";
    #endif
                        }
    #ifndef N_MATCHING_DEBUG
        //				cout << "("<<e->source()->index()<<","<<e->target()->index()<<") ";
    #endif
                    }
    #ifndef N_MATCHING_DEBUG
    //				else
    //				{
    //					cout << "X("<<e->source()->index()<<","<<e->target()->index()<<") ";
    //				}
    #endif
                }
            }
        }
    }
#ifndef N_MATCHING_DEBUG
    //cout <<endl;
#endif

    // hide edges not part of a directed cycle using strongy connected components
    //trim(enumG,enumMatched,enumNodeType,enumCopy,enumNonIsolatedVertices,enumHESH);
    m_EnumEdgeID[m_EnumerationGraph.size()-1]=enumG.numberOfEdges()-1;
    m_EnumNumberofNonIsolatedVertices[m_EnumerationGraph.size()-1]= static_cast<unsigned>(enumNonIsolatedVertices.size());
    trim(static_cast<unsigned int>(m_EnumerationGraph.size())-1);
}
int MWM::addEgdeToFixedMatching(unsigned jobIndex, const edge e, bool matchingFlag)
{
#ifndef N_MATCHING_DEBUG
    cout << "Added fixed matching edges: ";
#endif
    EdgeArray<vector<edge>*>& deg2ReplacedPath=m_EnumDeg2ReplacedPath[jobIndex];
    vector<edge>& fixedMatching=m_EnumFixedMatching[jobIndex];
    const NodeArray<NodeType>& nodeType=m_EnumNodeType[jobIndex];
    int replacedEdges=0;
    if (deg2ReplacedPath[e]->size()%2==0) // Last edge in vector is assigned matched[e], therefore swap matchingFlag, if size is even
        matchingFlag = !matchingFlag;
    for (edge edgeInReplacedPath:*deg2ReplacedPath[e])
    {
#ifndef N_MATCHING_DEBUG
        cout << edgeInReplacedPath<<endl;
#endif
        if (deg2ReplacedPath[edgeInReplacedPath]==nullptr)
        {
            if (matchingFlag && nodeType[edgeInReplacedPath->source()]==NodeType::AGENT && nodeType[edgeInReplacedPath->target()]==NodeType::JOB)
            {
                fixedMatching.push_back(edgeInReplacedPath);
                ++replacedEdges;
            }
        }
        else
        {
            replacedEdges += addEgdeToFixedMatching(jobIndex, edgeInReplacedPath, matchingFlag);
        }
        matchingFlag = !matchingFlag;
    }
    return replacedEdges;
}

void MWM::trim(unsigned jobIndex, const vector<edge>* edgeToRemove)
{
    Graph &G=m_EnumerationGraph[jobIndex];
    EdgeArray<bool>& matched=m_EnumMatched[jobIndex];
    const NodeArray<NodeType>& nodeType=m_EnumNodeType[jobIndex];
    const NodeArray<node>& copy = m_EnumCopy[jobIndex];
    vector<node>& nonIsolatedVertex=m_EnumNonIsolatedVertices[jobIndex];
    unsigned& numberofNonIsolatedVertices = m_EnumNumberofNonIsolatedVertices[jobIndex];
    vector<EnumRestore*>& restore=m_EnumRestore[jobIndex];
    NodeArray<unsigned>& scc=m_EnumSCC[jobIndex];
    vector<edge>& fixedMatching=m_EnumFixedMatching[jobIndex];
    edge e;
//	if (edgeToRemove != nullptr)
    //	cout << "pre emplace: "<<edgeToRemove->size()<<" edges in E1/2:";
    restore.emplace_back(new EnumRestore(G));
    Graph::HiddenEdgeSet& hes=restore.back()->hes;
    vector<edge>& shortcutedges=restore.back()->shortcutedges;
    int &additionalFixedMatchingEdges=restore.back()->additionalFixedMatchingEdges;
    vector<edge>& trimmedEdges=restore.back()->trimmedEdges;
    EdgeArray<bool>& deg2TrimmedEdges=m_EnumDeg2TrimmedEdges[jobIndex];
    EdgeArray<vector<edge>*>& deg2ReplacedPath=m_EnumDeg2ReplacedPath[jobIndex];
    // hide edges in edgeToRemove and add matching edges therein to m_EnumFixedMatching
    if (edgeToRemove != nullptr)
    {
#ifndef N_MATCHING_DEBUG
        cout << "Hiding "<<edgeToRemove->size()<<" edges in E1/2:";
#endif
        for(edge e:*edgeToRemove)
        {
#ifndef N_MATCHING_DEBUG
            cout << e;
#endif
            if (deg2ReplacedPath[e]==nullptr)
            {
                if (matched[e] && nodeType[e->source()]==NodeType::AGENT && nodeType[e->target()]==NodeType::JOB)
                {
                    fixedMatching.push_back(e);
                    ++ additionalFixedMatchingEdges;
                }
            }
            else
                additionalFixedMatchingEdges += addEgdeToFixedMatching(jobIndex,e,matched[e]);
            //G.hideEdge(hesh,e);
            hes.hide(e);
        }
#ifndef N_MATCHING_DEBUG
        cout <<endl;
#endif
    }
    /*// remove vertices of degree 0
    for (unsigned i=nonIsolatedVertex.size(); i>0; --i)
    {
        if (nonIsolatedVertex[i-1]->degree()==0 && copy[nonIsolatedVertex[i-1]]->degree()==0)
        {
            restore.back()->isolatedVertices.push_back(nonIsolatedVertex[i-1]);
            if (i != nonIsolatedVertex.size())
                swap(nonIsolatedVertex[i-1],nonIsolatedVertex.back());
            nonIsolatedVertex.pop_back();
        }
    }*/
    // remove vertices of degree 0
    for (int i=numberofNonIsolatedVertices-1; i>=0; --i)
    {
        if (nonIsolatedVertex[i]->degree()==0 && copy[nonIsolatedVertex[i]]->degree()==0)
        {
            ++(restore.back()->removedIsolatedVertices);
            if (i+1 != (int)numberofNonIsolatedVertices)
                swap(nonIsolatedVertex[i],nonIsolatedVertex[numberofNonIsolatedVertices-1]);
            --numberofNonIsolatedVertices;
        }
    }

    // compute SCC numbers
    computeSCCs(jobIndex);
    // Hide edges between different strongly connected components
    vector<edge> edgesToRemove;
    forall_edges(e,G) // remaining not hidden edges
    {
        node s=e->source();
        node t=e->target();
        if (scc[s]!=scc[t])
        {
            edgesToRemove.push_back(e);
            if (deg2ReplacedPath[e]==nullptr)
            {
                if (matched[e] && nodeType[s]==NodeType::AGENT && nodeType[t]==NodeType::JOB)
                {
                    fixedMatching.push_back(e);
                    ++ additionalFixedMatchingEdges;
                }
            }
            else
            {
                additionalFixedMatchingEdges += addEgdeToFixedMatching(jobIndex,e,matched[e]);
            }
        }
    }
    for (edge e:edgesToRemove)
    {
        //G.hideEdge(hesh,e);
        hes.hide(e);
    }
    edgesToRemove.clear();

#ifndef N_MATCHING_DEBUG
    cout << "removed isolated vertices, E and edges between SCCs"<<endl;
#endif
    //displayEnum(jobIndex);
#ifndef N_MATCHING_DEBUG
    cout <<"Shrinking paths of degree 2\n";
#endif
    // Shrink consecutive degree 2 vertices
    //vector<edge> agentToAgentCopyToRemove, JobCopyToJobToRemove;
    vector<edge> edgeBetweenGraphAndCopyOnCircle,edgeInOriginalGraphOnCircle;
    forall_edges(e,G) // remaining unhidden edges
    {
        node s=e->source();
        node t=e->target();
        if (copy[t]==s)
        {
            if (nodeType[s]==NodeType::AGENT && s->degree()==2)
                edgeBetweenGraphAndCopyOnCircle.push_back(e);
            else if (nodeType[s]==NodeType::JOBCOPY && t->degree()==2)
                edgeBetweenGraphAndCopyOnCircle.push_back(e);
        }
        else
        {
            if (t->degree()==2 && s->degree()==2)
                edgeInOriginalGraphOnCircle.push_back(e);
        }
    }

    // search paths/circles of length at least 4 and shrink them - starting from edges between original and copy
    for (edge e:edgeBetweenGraphAndCopyOnCircle)
    {
        edge neighbourEdge;
        node neighbourNode=nullptr;
        node pathNode=e->source();
        node startNode=e->target();
        if (nodeType[pathNode]==NodeType::JOBCOPY)
            swap(pathNode,startNode);
        if (pathNode->degree()!=2 || startNode->degree()!=1 || deg2TrimmedEdges[e]) // edge or adjacent edges already hidden
            continue;
        edge pathEdge=e;
        getOtherEdgeAndNode(pathNode,pathEdge,neighbourNode,neighbourEdge);
        if (neighbourNode == startNode) // circle of length 2
            continue;
        vector<edge>* newPath=new vector<edge>;
        newPath->push_back(e);
        //G.hideEdge(hesh,e);
        hes.hide(e);
        deg2TrimmedEdges[e]=true;
        trimmedEdges.push_back(e);
        // construct the path from startNode to a node of degree >2 or a copy node
        while (copy[pathNode]!=neighbourNode && neighbourNode->degree()==2)
        {
            newPath->push_back(neighbourEdge);
            //G.hideEdge(hesh,neighbourEdge);
            hes.hide(neighbourEdge);
            deg2TrimmedEdges[neighbourEdge]=true;
            trimmedEdges.push_back(neighbourEdge);
            pathNode=neighbourNode;
            pathEdge=neighbourEdge;
            getOtherEdgeAndNode(pathNode,pathEdge,neighbourNode,neighbourEdge);
        }
        // new shortcut edge
        if (nodeType[neighbourNode]==NodeType::AGENT || nodeType[neighbourNode]==NodeType::JOBCOPY)
            shortcutedges.push_back(G.newEdge(neighbourNode,copy[neighbourNode],++m_EnumEdgeID[jobIndex]));
        else
            shortcutedges.push_back(G.newEdge(copy[neighbourNode],neighbourNode,++m_EnumEdgeID[jobIndex]));

        if (neighbourNode->degree()>2) // path leading to a vertex of degree 3 or larger.
        {
            newPath->push_back(neighbourEdge);
            //G.hideEdge(hesh,neighbourEdge);
            hes.hide(neighbourEdge);
            deg2TrimmedEdges[neighbourEdge]=true;
            trimmedEdges.push_back(neighbourEdge);
            matched[shortcutedges.back()]=matched[neighbourEdge];
        }
        else // cycle between graph and copy.
        {
            matched[shortcutedges.back()]=!matched[neighbourEdge];
        }
        deg2ReplacedPath[shortcutedges.back()]=newPath;
        // mark edges on path as hidden
        //for (edge eOnPath:*newPath)
    //		deg2TrimmedEdges[eOnPath]=true;
    }

    // search paths/circles of length at least 4 and shrink them - starting from edges on a circle in the original graph
    for (edge e:edgeInOriginalGraphOnCircle)
    {
        edge neighbourEdge;
        node neighbourNode=nullptr;
        node pathNode=e->target();
        node startNode=e->source();
        if (pathNode->degree()!=2 || startNode->degree()!=2 || deg2TrimmedEdges[e]) // edge already hidden
            continue;
        edge pathEdge=e;
        getOtherEdgeAndNode(pathNode,pathEdge,neighbourNode,neighbourEdge);
        if (neighbourNode == startNode) // circle of length 2
            continue;
        vector<edge>* newPath=new vector<edge>;
        newPath->push_back(e);
        //G.hideEdge(hesh,e);
        hes.hide(e);
        deg2TrimmedEdges[e]=true;
        trimmedEdges.push_back(e);
        // construct the path from startNode to a node of degree >2 or startNode again
        while (neighbourNode->degree()==2 && neighbourNode != startNode)
        {
            newPath->push_back(neighbourEdge);
            //G.hideEdge(hesh,neighbourEdge);
            hes.hide(neighbourEdge);
            deg2TrimmedEdges[neighbourEdge]=true;
            trimmedEdges.push_back(neighbourEdge);
            pathNode=neighbourNode;
            pathEdge=neighbourEdge;
            getOtherEdgeAndNode(pathNode,pathEdge,neighbourNode,neighbourEdge);
        }
        if (neighbourNode == startNode) // circle, therefore add shortcut between startNode and last pathNode for the path (of uneven length) traveled
        {
            shortcutedges.push_back(G.newEdge(startNode,pathNode,++m_EnumEdgeID[jobIndex]));
            matched[shortcutedges.back()]=matched[pathEdge];
        }
        else // reached node of degree 3 or larger. Therefore inverse traveled path and extend into other direction
        {
            // add edge leading to the node of degree 3 or larger
            newPath->push_back(neighbourEdge);
            //G.hideEdge(hesh,neighbourEdge);
            hes.hide(neighbourEdge);
            // reverse path
            reverse(newPath->begin(), newPath->end());
            // set new startNode
            startNode=neighbourNode;
            // extend into other direction
            pathNode=e->source();
            pathEdge=e;
            getOtherEdgeAndNode(pathNode,pathEdge,neighbourNode,neighbourEdge);
            while (neighbourNode->degree()==2)
            {
                newPath->push_back(neighbourEdge);
                //G.hideEdge(hesh,neighbourEdge);
                hes.hide(neighbourEdge);
                deg2TrimmedEdges[neighbourEdge]=true;
                trimmedEdges.push_back(neighbourEdge);
                pathNode=neighbourNode;
                pathEdge=neighbourEdge;
                getOtherEdgeAndNode(pathNode,pathEdge,neighbourNode,neighbourEdge);
            }
            // path length must be uneven. Therefore, we add edge leading to deg >2 node exactly if current path has even length
            if (newPath->size()%2==0)
            {
                newPath->push_back(neighbourEdge);
                //G.hideEdge(hesh,neighbourEdge);
                hes.hide(neighbourEdge);
                deg2TrimmedEdges[neighbourEdge]=true;
                trimmedEdges.push_back(neighbourEdge);
                if (nodeType[startNode]==NodeType::AGENT)
                    shortcutedges.push_back(G.newEdge(startNode,neighbourNode,++m_EnumEdgeID[jobIndex]));
                else
                    shortcutedges.push_back(G.newEdge(neighbourNode,startNode,++m_EnumEdgeID[jobIndex]));
                matched[shortcutedges.back()]=matched[neighbourEdge];
            }
            else // path length already uneven. Therefore do not add last edge. Connect startNode with pathNode
            {
                if (nodeType[startNode]==NodeType::AGENT)
                    shortcutedges.push_back(G.newEdge(startNode,pathNode,++m_EnumEdgeID[jobIndex]));
                else
                    shortcutedges.push_back(G.newEdge(pathNode,startNode,++m_EnumEdgeID[jobIndex]));
                matched[shortcutedges.back()]=matched[pathEdge];
            }
        }
        deg2ReplacedPath[shortcutedges.back()]=newPath;
        // hide all the edges on the path
        //for (edge eOnPath:*newPath)
        //	deg2TrimmedEdges[eOnPath]=true;
    }

    // remove vertices of degree 0
    for (int i=numberofNonIsolatedVertices-1; i>=0; --i)
    {
        if (nonIsolatedVertex[i]->degree()==0 && copy[nonIsolatedVertex[i]]->degree()==0)
        {
            ++restore.back()->removedIsolatedVertices;
            if (i+1 != (int) numberofNonIsolatedVertices)
                swap(nonIsolatedVertex[i],nonIsolatedVertex[numberofNonIsolatedVertices-1]);
            --numberofNonIsolatedVertices;
        }
    }

    // show current fixed matching
#ifndef N_MATCHING_DEBUG
    cout << "Fixed matching: ";
    for (auto it=m_EnumFixedMatching[jobIndex].begin(); it!= m_EnumFixedMatching[jobIndex].end();++it)
    {
        cout << *it<<" ";
    }
    cout <<endl;
#endif

    // recompute scc Numbers
    computeSCCs(jobIndex);
    //displayEnum(m_EnumerationGraph.size()-1);
}

void MWM::computeSCCs(unsigned jobIndex)
{
    //Graph &G=m_EnumerationGraph[jobIndex];
    vector<node>& nonIsolatedVertex=m_EnumNonIsolatedVertices[jobIndex];
    NodeArray<short> &dfsstatus=m_EnumDFSstatus[jobIndex];
    unsigned& numberofNonIsolatedVertices = m_EnumNumberofNonIsolatedVertices[jobIndex];
    for (unsigned i=0; i<numberofNonIsolatedVertices; ++i)
        dfsstatus[nonIsolatedVertex[i]]=0;
//	for (node n:nonIsolatedVertex)
    //	dfsstatus[n]=0;
    stack<node> dfsstack;
    node stacknode,originalStacknode,child,n;
    bool isCopy,requireMatchingEdge;
    const NodeArray<node>& copy = m_EnumCopy[jobIndex];
    vector<node> fstack; // Nodes inserted with descending f-numbers. Thus we can compute second dfs run from last to first element
    fstack.reserve(numberofNonIsolatedVertices);
    adjEntry adj;
    EdgeArray<bool>& matched=m_EnumMatched[jobIndex];
    NodeArray<unsigned>& scc=m_EnumSCC[jobIndex];
    const NodeArray<NodeType>& nodeType=m_EnumNodeType[jobIndex];

// first dfs run
#ifndef N_MATCHING_DEBUG
    cout << endl<<"First run: ";
#endif
    for (unsigned i=0; i<numberofNonIsolatedVertices; ++i)
    {
        node n=nonIsolatedVertex[i];
        // invariant: stack empty. Thus start new dfs search from next unvisited node
        if (dfsstatus[n]==0)
        {
            dfsstack.push(n);
#ifndef N_MATCHING_DEBUG
            cout << "p"<<n<< " ";
#endif
            //++dfsstatus[n];
        }

        // if stack is empty, insert the next not visited node
        while (!dfsstack.empty())
        {
            stacknode=dfsstack.top();
            if (nodeType[stacknode]<=NodeType::JOB) // agent or job
            {
                originalStacknode=stacknode;
                isCopy=false;
            }
            else
            {
                isCopy=true;
                originalStacknode=copy[stacknode];
            }

            if (dfsstatus[stacknode]==1) // recursion for this node finished
            {
                ++dfsstatus[stacknode];
                fstack.push_back(stacknode);
#ifndef N_MATCHING_DEBUG
                cout << "f"<<stacknode->index()<< " ";
#endif
                dfsstack.pop();
            }
            else if (dfsstatus[stacknode]==0) // first time visited . Thus add children to stack, which were not visited yet on stack
            {
                ++dfsstatus[stacknode];
                forall_adj(adj,originalStacknode)
                {
                    if (nodeType[stacknode]==NodeType::AGENT || nodeType[stacknode]==NodeType::JOBCOPY) // first dfs run: from agent or jobcopy non matching edge; else matching
                        requireMatchingEdge=false;
                    else
                        requireMatchingEdge=true;
                    child=adj->twinNode();
                    if (isCopy)
                        child=copy[child];
                    if (matched[adj->theEdge()]==requireMatchingEdge)
                    {
                        if (dfsstatus[child]==0)
                        {
#ifndef N_MATCHING_DEBUG
                            cout << "f"<<stacknode->index()<<","<<child->index()<<" ";
#endif
                            dfsstack.push(child);
                        }
#ifndef N_MATCHING_DEBUG
                        else if (dfsstatus[child]==1)
                            cout << "b"<<stacknode->index()<<","<<child->index()<<" ";
                        else
                            cout << "c"<<stacknode->index()<<","<<child->index()<<" ";
#endif
                    }
                }
            }
            else
            {
#ifndef N_MATCHING_DEBUG
                cout << "fb"<<stacknode->index()<< " ";
#endif
                dfsstack.pop();
            }
        }
    }
#ifndef N_MATCHING_DEBUG
    cout << endl<<"Second run: ";
#endif
    // second dfs run
    int scc_number=-1;
    for (unsigned i=0; i<numberofNonIsolatedVertices; ++i)
        dfsstatus[nonIsolatedVertex[i]]=0;
    //for (node n:nonIsolatedVertex)
        //dfsstatus[n]=0;
    while (!fstack.empty())
    {
        n=fstack.back();
        // invariant: stack empty. Thus start new dfs search from next unvisited node
        if (dfsstatus[n]==0) // not visited
        {
            dfsstack.push(n);
#ifndef N_MATCHING_DEBUG
            cout << "p"<<n<< " ";
#endif
            ++scc_number; // new strongly connected component
        }

        // if stack is empty, insert the next not visited node
        while (!dfsstack.empty())
        {
            stacknode=dfsstack.top();
            if (nodeType[stacknode]<=NodeType::JOB) // agent or job
            {
                isCopy=false;
                originalStacknode=stacknode;
            }
            else
            {
                isCopy=true;
                originalStacknode=copy[stacknode];
            }

            if (dfsstatus[stacknode]==1) // recursion for this node finished
            {
                ++dfsstatus[stacknode];
#ifndef N_MATCHING_DEBUG
                cout << "f"<<stacknode->index()<< " ";
#endif
                dfsstack.pop();
            }
            else if (dfsstatus[stacknode]==0) // first time visited. Thus add children to stack, which were not visited yet
            {
                scc[stacknode]=scc_number;
                ++dfsstatus[stacknode];
                forall_adj(adj,originalStacknode)
                {
                    if (nodeType[stacknode]==NodeType::AGENT || nodeType[stacknode]==NodeType::JOBCOPY) // second dfs run: from agent to job matching edge; from J to A non matching
                        requireMatchingEdge=true;
                    else
                        requireMatchingEdge=false;
                    child=adj->twinNode();
                    if (isCopy)
                        child=copy[child];
                    if (matched[adj->theEdge()]==requireMatchingEdge)
                    {
                        if (dfsstatus[child]==0)
                        {
#ifndef N_MATCHING_DEBUG
                            cout << "f"<<stacknode->index()<<","<<child->index()<<" ";
#endif
                            dfsstack.push(child);
                        }
#ifndef N_MATCHING_DEBUG
                        else if (dfsstatus[child]==1)
                            cout << "b"<<stacknode->index()<<","<<child->index()<<" ";
                        else
                            cout << "c"<<stacknode->index()<<","<<child->index()<<" ";
#endif
                    }
                }
            }
            else
            {
#ifndef N_MATCHING_DEBUG
                cout << "fb"<<stacknode->index()<< " ";
#endif
                dfsstack.pop();
            }
        }
        fstack.pop_back();
    }
    m_EnumMaxSCCNumber[jobIndex]=scc_number;

#ifndef N_MATCHING_DEBUG
    cout << "computed SCCs\nNon Isolated Vertices: ";
    for (unsigned i=0; i<numberofNonIsolatedVertices; ++i)
    //for (node n:nonIsolatedVertex)
    {
        node n=nonIsolatedVertex[i];
        cout << n<<":"<<scc[n]<<" ";
    }
    cout << endl<<endl;
#endif
    //displayEnum(jobIndex);
}

void MWM::getOtherEdgeAndNode(const node n, const edge e, node& othernode, edge& otheredge)
{
    adjEntry adj;
    forall_adj(adj,n)
    {
        otheredge=adj->theEdge();
        if (otheredge==e) // need to find other edge
            continue;
        if (otheredge->source()==n)
            othernode=otheredge->target();
        else
            othernode=otheredge->source();
        return;
    }
}

void MWM::getOtherEdgeAndNodeWithDifferentMatchedValue(const EdgeArray<bool>& matched, const node n, const edge e, node& othernode, edge& otheredge)
{
    adjEntry adj;
    forall_adj(adj,n)
    {
        otheredge=adj->theEdge();
        if (otheredge==e) // need to find other edge
            continue;
        if (matched[e]==matched[otheredge])
            continue;
        if (otheredge->source()==n)
            othernode=otheredge->target();
        else
            othernode=otheredge->source();
        return;
    }
}

wType MWM::getWeight(unsigned jobIndex)
{
    if (m_computedMatchings == -1) // First MWM. Compute all matchings, if jobIndex=m_nodesJ, else only sub matching
    {
        //cout << "first call: "<<jobIndex<<"\n";
        compute(jobIndex);
    }
    else if (m_computedMatchings != (int) jobIndex && m_computedMatchings < (int) m_nodesJ) // 2nd not yet computed matching. Then compute all MWMs
    {
        //cout << "second call: "<<jobIndex<<"\n";
        compute(m_nodesJ);
    }
    return m_weight[m_refToMatchings[jobIndex]];
}

void MWM::compute(unsigned jobIndex)
{
    // copy graph if not done yet (reduction MWM -> MWPM)
    addGraphCopy();
    //m_mate.init(m_G);
    edge e;
    mco::HungarianMethod solver;
    if (m_computedMatchings == -1)
    {
        m_weight.reserve(m_nodesJ+1);
        m_refToMatchings.resize(m_nodesJ+1); // in total nodesJ+1 matchings
        m_matchingAgentAndJob.reserve(min(m_nodesA,m_nodesJ)+1); // at most min(m_nodesA,m_nodesJ)+1 of them are different
        // Data structures for enumeration. For each excluded vertex and the full instance there is a (different) trimmed equality subgraph
        if (m_withEnumeration)
        {
            m_EnumerationGraph.reserve(m_nodesJ+1);
            m_EnumMatched.reserve(m_nodesJ+1);
            m_EnumNodeType.reserve(m_nodesJ+1);
            m_EnumCopy.reserve(m_nodesJ+1);
            m_EnumNonIsolatedVertices.reserve(m_nodesJ+1);
    #ifdef GRAPHICS
            m_EnumRenderWindow.reserve(m_nodesJ+1);
    #endif
            m_EnumRestore.reserve(m_nodesJ+1);
            m_EnumFixedMatching.reserve(m_nodesJ+1);
            m_EnumDeg2TrimmedEdges.reserve(m_nodesJ+1);
            m_EnumDeg2ReplacedPath.reserve(m_nodesJ+1);
            m_EnumStack.reserve(m_nodesJ+1);
            m_EnumSCC.reserve(m_nodesJ+1);
            m_EnumMaxSCCNumber.reserve(m_nodesJ+1);
            m_EnumDFSstatus.reserve(m_nodesJ+1);
            m_EnumPreviousEdgeOnCircle.reserve(m_nodesJ+1);
            m_EnumEdgeID.resize(m_nodesJ+1);
            m_EnumNumberofNonIsolatedVertices.resize(m_nodesJ+1);
            if (m_EnumStoreMatchings)
            {
                m_EnumStoredMatching.reserve(m_nodesJ+1);
                m_EnumStoredMatchings.reserve(m_nodesJ+1);
                m_EnumCurStoredMatching.reserve(m_nodesJ+1);
            }
        }
    }
#ifndef N_MATCHING_DEBUG
    cout << "Kantengewichte  ";
    forall_edges(e,m_G)
    {
        //if (m_egdeweight[e]>0 && e->source()->index()<(signed) m_nodesA)
        //if (m_egdeweight[e]>0 && m_nodeType[e->source()]==NodeType::AGENT)
        if (m_nodeType[e->source()]==NodeType::AGENT)
            cout <<"A"<<" "<< e<<":"<<m_egdeweight[e]<< ", ";
        else
            cout <<"J"<<" "<< e<<":"<<m_egdeweight[e]<< ", ";
            //cout <<"("<< originalAgent[e->source()->index()]<<","<< originalJob[e->target()->index()]<<")"<<":"<<egdeweight[e]<< ", ";
    }
    cout << endl;
#endif

    m_matchingAgentAndJob.push_back(vector<AgentJobPair>());

    // remove edges if first matching is not full matching; in LAWESCU this is not necessary
    Graph::HiddenEdgeSet hes(m_G);
    if (jobIndex != m_nodesJ)
    {
        node nodeToRemove=m_job[jobIndex];
        for(adjEntry adj : nodeToRemove->adjEntries)
        {
            e=adj->theEdge();
            //cout << "hide 1 "<<m_nodeType[e->adjSource()->theNode()]<<" \n";
            if (m_nodeType[e->adjSource()->theNode()] == NodeType::AGENT)
            {
                //cout << "hide 2\n";
                hes.hide(e);
            }
        }
        for(adjEntry adj : m_agent[jobIndex+m_nodesA]->adjEntries)
        {
            e=adj->theEdge();
            //cout << "hide 3 "<<m_nodeType[e->adjTarget()->theNode()]<<" \n";
            if (m_nodeType[e->adjTarget()->theNode()] == NodeType::AGENTCOPY)
            {
                //cout << "hide 4\n";
                hes.hide(e);
            }
        }
    }



    m_refToMatchings[jobIndex]= static_cast<unsigned int>(m_weight.size()); // reference of first sub matching (=0) or of main matching (=1)

    // Initial mates and duals
    //m_dual.init(m_G,0);
    if (m_computedMatchings == -1) // first call of compute(...)
    {
        if (m_nodesA <= m_nodesJ)
        {
            for (unsigned i=0; i<m_nodesA; ++i)
            {
                node agent=m_agent[i];
                edge minimum = agent->firstAdj()->theEdge();
                for(adjEntry adj : agent->adjEntries)
                {
                    e=adj->theEdge();
                //forall_adj_edges(e, agent)
                    if (m_egdeweight[e]>m_egdeweight[minimum])
                    {
                        minimum = e;
                    }
                }
                m_dual[agent] = -m_egdeweight[minimum];
                m_dual[m_job[i+m_nodesJ]]= -m_egdeweight[minimum];
            }
            // compute partners and remaining duals
            for (unsigned i=0; i<m_nodesJ; ++i)
            {
                m_mate[m_agent[i+m_nodesA]]=m_job[i];
                m_mate[m_job[i]]=m_agent[i+m_nodesA];
                m_dual[m_job[i]]=0;
                m_dual[m_copy[m_job[i]]]=0;
            }
            // exponierte Knoten
            for (unsigned i=0; i<m_nodesA; ++i)
            {
                m_mate[m_agent[i]]=m_agent[i];
                m_mate[m_job[i+m_nodesJ]]=m_job[i+m_nodesJ];
            }
        }
        else
        {
            for (unsigned i=0; i<m_nodesJ; ++i)
            {
                node job=m_job[i];
                edge minimum = job->firstAdj()->theEdge();
                for(adjEntry adj : job->adjEntries)
                {
                    e=adj->theEdge();
                //forall_adj_edges(e, job)
                    if (m_egdeweight[e]>m_egdeweight[minimum])
                    {
                        minimum = e;
                    }
                }
                m_dual[job] = -m_egdeweight[minimum];
                m_dual[m_agent[i+m_nodesA]]= -m_egdeweight[minimum];
            }
            // Partner zuweisen
            for (unsigned i=0; i<m_nodesA; ++i)
            {
                m_mate[m_job[i+m_nodesJ]]=m_agent[i];
                m_mate[m_agent[i]]=m_job[i+m_nodesJ];
                m_dual[m_agent[i]]=0;
                m_dual[m_copy[m_agent[i]]]=0;
            }
            // exponierte Knoten
            for (unsigned i=0; i<m_nodesJ; ++i)
            {
                m_mate[m_job[i]]=m_job[i];
                m_mate[m_agent[i+m_nodesA]]=m_agent[i+m_nodesA];
            }
        }
        m_weight.push_back(solver.Solve_Max(m_G, m_egdeweight, m_agent, m_job, m_mate, m_dual, min(m_nodesA,m_nodesJ)));
    }
    else // else case only in second call of compute(unsigned jobIndex). Then 2nd MWM is a main matching
        // we need to adjust duals and expone the previously skipped jop and its copy
    {
        node nodeToInsert=m_job[m_computedMatchings];
        // expone job and its copy; they were previously connected
        m_mate[nodeToInsert]=nodeToInsert;
        m_mate[m_copy[nodeToInsert]]=m_copy[nodeToInsert];

        // adjust dual of skipped job
        edge minimum = nodeToInsert->firstAdj()->theEdge();
        for(adjEntry adj : nodeToInsert->adjEntries)
        {
            e=adj->theEdge();
        //forall_adj_edges(e, job)
            if (m_egdeweight[e]+m_dual[e->source()] > m_egdeweight[minimum]+m_dual[minimum->source()])
            {
                minimum = e;
            }
        }
        m_dual[nodeToInsert] = -m_egdeweight[minimum]-m_dual[minimum->source()];

        // adjust dual of copy of skipped job
        nodeToInsert=m_copy[nodeToInsert];
         minimum = nodeToInsert->firstAdj()->theEdge();
        for(adjEntry adj : nodeToInsert->adjEntries)
        {
            e=adj->theEdge();
        //forall_adj_edges(e, job)
            if (m_egdeweight[e]+m_dual[e->target()]>m_egdeweight[minimum]+m_dual[minimum->target()])
            {
                minimum = e;
            }
        }
        m_dual[nodeToInsert] = -m_egdeweight[minimum]-m_dual[minimum->target()];

        m_weight.push_back(solver.Solve_Max(m_G, m_egdeweight, m_agent, m_job, m_mate, m_dual, 1));
    }

    if (jobIndex != m_nodesJ)
        hes.restore();

    //m_refToMatchings[m_nodesJ]=0; // main matching has reference 0

#ifndef N_MATCHING_DEBUG
    cout << "Matchingkanten: ";
#endif
    // determine matching edges of input graph, i.e., the original edges of the computed matching
    for (unsigned i=0; i<m_nodesA; ++i)
    {
        if (m_nodeType[m_mate[m_agent[i]]]==NodeType::JOB)
        {
            m_matchingAgentAndJob.back().push_back(AgentJobPair(m_vMatchingTovOriginal[m_agent[i]],m_vMatchingTovOriginal[m_mate[m_agent[i]]]));
#ifndef N_MATCHING_DEBUG
            cout << "(" << m_agent[i] << "," << m_mate[m_agent[i]] << "):"<< -m_dual[m_agent[i]]-m_dual[m_mate[m_agent[i]]]<<", ";
#endif
        }
    }
#ifndef N_MATCHING_DEBUG
    cout << " Wert: "<<solver.value()<<endl<<endl;
#endif

    if (jobIndex < m_nodesJ) // if the first matching is a sub matching, no more computation for now.
    {
        //todo : enum code
        m_computedMatchings = jobIndex;
        return;
    }

    // Mates und Duale Variablen der vollständiges Lösung speichern
    NodeArray<node> mate_full(m_mate);
    NodeArray<wType> dual_full(m_dual);

    // Main Enumeration instance
    //if (m_withEnumeration)
        //constructEqualitySubgraph(nullptr,m_nodesJ+m_nodesA);

    // Compute sub matchings
    // if (nodesJ>1) // never happens, because CVs are connected to at least 2 B-Nodes
    unsigned refCounter= static_cast<unsigned int>(m_weight.size())-1;  // ref of main matching (also current last) depends on whether it was computed first or second.
    unsigned mainMWMref=refCounter;
    for (unsigned job=0; job<m_nodesJ; ++job)
    {
        if (job == jobIndex) // we computed this sub matching in the first step
            continue;
        node nodeToRemove=m_job[job];
        // Nicht mit Agent-Knoten verbunden
        //if (mate_full[nodeToRemove]->index()>=(signed) m_nodesA)
        if (m_nodeType[mate_full[nodeToRemove]] != NodeType::AGENT)

        {
#ifndef N_MATCHING_DEBUG
            cout << "Knoten "<<nodeToRemove<< " ist kein Teil des echten Matchings.\n\n";
#endif
            m_refToMatchings[job] = mainMWMref; // in this case full matching and matching without vertex m_job[job] are the same
            if (m_withEnumeration)
                constructEqualitySubgraph(nodeToRemove,m_nodesJ+m_nodesA); // but trimmed equality subgraph may differ
        }
        else
        {
#ifndef N_MATCHING_DEBUG
            cout << "Knoten "<<nodeToRemove<< " wird aus Graph entfernt.\n";
#endif
            // anderen Knoten der Matchingkante exponieren
            m_mate[m_mate[nodeToRemove]]=m_mate[nodeToRemove];
            // Knoten, der zu Hilfsknoten durch Matchingkante verbunden, exponieren
            m_mate[m_mate[m_agent[job+m_nodesA]]]=m_mate[m_agent[job+m_nodesA]];

            // alle Kanten, die zum entfernenden Knoten und dessen Hilfsknoten inzidient sind, verstecken
            //Graph::HiddenEdgeSetHandle HESH=m_G.newHiddenEdgeSet();
            Graph::HiddenEdgeSet HES(m_G);
            for(adjEntry adj : nodeToRemove->adjEntries)
            {
                e=adj->theEdge();
                HES.hide(e);
            }
            //forall_adj_edges(e, nodeToRemove)
                //m_G.hideEdge(HESH,e);
            for(adjEntry adj : m_agent[job+m_nodesA]->adjEntries)
            {
                e=adj->theEdge();
                HES.hide(e);
            }
            //forall_adj_edges(e, m_agent[job+m_nodesA])
                //m_G.hideEdge(HESH,e);

            // zu entfernenden Knoten und dessen Hilfsknoten aus den agents und jobs entfernen
            node removedAgent=m_agent[job+m_nodesA];
            node removedJob=m_job[job];
            m_agent[job+m_nodesA]=m_agent.back();
            m_agent.erase(m_agent.end()-1);
            m_job[job]=m_job.back();
            m_job.erase(m_job.end()-1);
#ifndef N_MATCHING_DEBUG
            cout << "\nAgents. ";
            for(node n : m_agent)
                std::cout << "(" << n << ", " << m_mate[n] << "), ";
            cout << "\nJobs. ";
            for(node n : m_job)
                std::cout << "(" << n << ", " << m_mate[n] << "), ";

    cout << "Kantengewichte  ";
    forall_edges(e,m_G)
    {
        //if (m_egdeweight[e]>0 && e->source()->index()<(signed) m_nodesA)
        //if (m_egdeweight[e]>0 && m_nodeType[e->source()]==NodeType::AGENT)
        if (m_nodeType[e->source()]==NodeType::AGENT)
            cout <<"A"<<" "<< e<<":"<<m_egdeweight[e]<< ", ";
        else
            cout <<"J"<<" "<< e<<":"<<m_egdeweight[e]<< ", ";
            //cout <<"("<< originalAgent[e->source()->index()]<<","<< originalJob[e->target()->index()]<<")"<<":"<<egdeweight[e]<< ", ";
    }
    cout << endl;
#endif

            // Eine Iteration der ungarischen Methode
            m_refToMatchings[job]=++refCounter; // New reference for new matching
            m_weight.push_back(solver.Solve_Max(m_G, m_egdeweight, m_agent, m_job, m_mate, m_dual,1));

            // Ergebnis
            m_matchingAgentAndJob.push_back(vector<AgentJobPair>());
#ifndef N_MATCHING_DEBUG
            cout << "Matchingkanten: ";
#endif
            for (unsigned i=0; i<m_nodesA; ++i)
            {
                if (m_nodeType[m_mate[m_agent[i]]]==NodeType::JOB)
                {
                    m_matchingAgentAndJob[refCounter].push_back(AgentJobPair(m_vMatchingTovOriginal[m_agent[i]],m_vMatchingTovOriginal[m_mate[m_agent[i]]]));
#ifndef N_MATCHING_DEBUG
                    cout << "(" << m_vMatchingTovOriginal[m_agent[i]] << "," << m_vMatchingTovOriginal[m_mate[m_agent[i]]] << "):"<< -m_dual[m_agent[i]]-m_dual[m_mate[m_agent[i]]]<<", ";
#endif
                }
            }
#ifndef N_MATCHING_DEBUG
            cout << " Wert: "<<solver.value()<<endl<<endl;
#endif

            if (m_withEnumeration)
                constructEqualitySubgraph(nullptr,m_nodesJ+m_nodesA-1);
            // Ursprünglichen Graphen wiederherstellen
            //m_G.restoreEdges(HESH);
            HES.restore();
            m_mate=mate_full;
            m_dual=dual_full;
            if (m_agent.size() == job + m_nodesA)
            {
                m_agent.push_back(removedAgent);
            }
            else 
            {
                m_agent.push_back(m_agent[job + m_nodesA]);
                m_agent[job + m_nodesA] = removedAgent;
            }
            m_job.push_back(m_job[job]);
            m_job[job]=removedJob;
        }
    }
    // Main Enumeration instance
    if (m_withEnumeration)
        constructEqualitySubgraph(nullptr,m_nodesJ+m_nodesA);

    m_computedMatchings = m_nodesJ; // all matchings computed
}

void MWM::displayEnum(unsigned jobIndex)
{
//	return;
#ifdef GRAPHICS
    Graph& G=m_EnumerationGraph[jobIndex];
    EdgeArray<bool>& enumMatched=m_EnumMatched[jobIndex];
    NodeArray<NodeType>& enumNodeType=m_EnumNodeType[jobIndex];
    EdgeArray<vector<edge>*> enumDeg2ReplacedPath=m_EnumDeg2ReplacedPath[jobIndex];
    float scaling=2.0;
    int imageSizex=max(m_nodesA,m_nodesJ)*60-20;
    int imageSizey=260;
    while (m_EnumRenderWindow.size()<=jobIndex)
    {
        string title="MWM ";
        title += to_string(jobIndex);
        m_EnumRenderWindow.push_back(new sf::RenderWindow (sf::VideoMode(scaling*imageSizex, scaling*imageSizey), title.c_str(), sf::Style::Titlebar | sf::Style::Close));
        m_EnumRenderWindow.back()->setPosition(sf::Vector2<int>(jobIndex*imageSizex*scaling+jobIndex*10,200));
    }
    sf::RenderWindow& window=*m_EnumRenderWindow[jobIndex];
    if (!window.isOpen())
    {
        string title="MWM ";
        title += to_string(jobIndex);
        window.create(sf::VideoMode(scaling*imageSizex, scaling*imageSizey), title.c_str(), sf::Style::Titlebar | sf::Style::Close);
        window.setPosition(sf::Vector2<int>(jobIndex*imageSizex*scaling+jobIndex*10,100));
    }
    window.clear(sf::Color(255,255,255));
    // draw borders for original and copied vertices
    sf::RectangleShape background;
    background.setFillColor(sf::Color(240,240,240));
    background.setPosition(5*scaling,20*scaling);
    background.setSize(sf::Vector2f((imageSizex-10)*scaling,(imageSizey/2-30)*scaling));
    window.draw(background);
    background.setPosition(5*scaling,(imageSizey/2+10)*scaling);
    window.draw(background);
    // draw not hidden edges
    edge e;
    forall_edges(e,G)
    {
        node s=e->source();
        node t=e->target();
        int offsetx=0;
        sf::VertexArray line(sf::Lines, 2);
        if (enumDeg2ReplacedPath[e])
        {
            if (enumMatched[e])
            {
                line[0].color=sf::Color(255,128,0,255); // replacement matching edges
                line[1].color=sf::Color(255,128,0,255);
                offsetx=4;
            }
            else
            {
                line[0].color=sf::Color(92,160,255,255); // replacement non matching
                line[1].color=sf::Color(92,160,255,255);
                offsetx=-4;
            }

        }
        else
        {
            if (enumMatched[e])
            {
                line[0].color=sf::Color(255,0,0,255); // matching edges red
                line[1].color=sf::Color(255,0,0,255);
            }
            else
            {
                line[0].color=sf::Color(92,92,255,255); // non matching blue
                line[1].color=sf::Color(92,92,255,255);
            }
        }
        if (enumNodeType[s]==NodeType::AGENT)
        {
            if (enumNodeType[t]==NodeType::JOB) // edges between agent and job vertices
            {
                for (int x=0; x<=1; ++x)
                {
                    for (int y=0; y<=1; ++y)
                    {
                        line[0].position = sf::Vector2f((s->index()*60+20+offsetx)*scaling+x,40*scaling+y); // in original part
                        line[1].position = sf::Vector2f(((t->index()-m_nodesA)*60+20+offsetx)*scaling+x,100*scaling+y);
                        window.draw(line);
                    }
                }
                line[0].position = sf::Vector2f((s->index()*60+20+offsetx)*scaling,220*scaling); // in copied part
                line[1].position = sf::Vector2f(((t->index()-m_nodesA)*60+20+offsetx)*scaling,160*scaling);
                window.draw(line);
            }
            else // edges between agent and agent copy
            {
                line[0].position = sf::Vector2f((s->index()*60+20+offsetx)*scaling,40*scaling);
                line[1].position = sf::Vector2f((s->index()*60+20+offsetx)*scaling,10*scaling);
                window.draw(line);
                line[0].position = sf::Vector2f((s->index()*60+20+offsetx)*scaling,220*scaling);
                line[1].position = sf::Vector2f((s->index()*60+20+offsetx)*scaling,250*scaling);
                window.draw(line);
            }
        }
        else // edges between job copy and job
        {
            line[0].position = sf::Vector2f(((t->index()-m_nodesA)*60+20+offsetx)*scaling,100*scaling);
            line[1].position = sf::Vector2f(((t->index()-m_nodesA)*60+20+offsetx)*scaling,160*scaling);
            window.draw(line);
        }
    }
    // draw fixed matching
    for (edge e:m_EnumFixedMatching[jobIndex])
    {
        node s=e->source();
        node t=e->target();
        sf::VertexArray line(sf::Lines, 2);
        line[0].color=sf::Color(250,120,120,255); // matching edges red, a bit lighter than regular
        line[1].color=sf::Color(250,120,120,255);
        for (int x=0; x<=1; ++x)
        {
            for (int y=0; y<=1; ++y)
            {
                line[0].position = sf::Vector2f((s->index()*60+20)*scaling+x,40*scaling+y); // in original part
                line[1].position = sf::Vector2f(((t->index()-m_nodesA)*60+20)*scaling+x,100*scaling+y);
                window.draw(line);
            }
        }
    }
    // draw vertices
    sf::CircleShape vertexShape;
    sf::Text nodeIndexText;
    string nodeIndexString;
    vertexShape.setRadius(10*scaling);
    vertexShape.setOutlineThickness(1*scaling);
    nodeIndexText.setCharacterSize(12*scaling);
    nodeIndexText.setFont(m_font_Arial);
    nodeIndexText.setFillColor(sf::Color::Black);
    node n;
    forall_nodes(n,G)
    {
        unsigned char alpha=255;
        if (n->degree()==0)
            alpha=100;
        vertexShape.setOutlineColor(sf::Color(0,0,0,alpha));

        if (enumNodeType[n]==NodeType::AGENT)
        {
            vertexShape.setPosition((n->index()*60+10)*scaling,30*scaling);
            if (n->degree()==0)
            {
                vertexShape.setFillColor(sf::Color(240,240,240,255)); // black background, to overdraw fixed matching edges
                vertexShape.setOutlineColor(sf::Color(240,240,240,255));
                window.draw(vertexShape);
                vertexShape.setOutlineColor(sf::Color(0,0,0,alpha));
            }
            vertexShape.setFillColor(sf::Color(0,255,0,alpha)); // agent
            window.draw(vertexShape);
            nodeIndexString=to_string(n->index());
            nodeIndexText.setPosition((n->index()*60+16)*scaling,32*scaling);
            nodeIndexText.setString(nodeIndexString);
            nodeIndexText.setFillColor(sf::Color(0,0,0,alpha));
            if (n->index()>9)
                nodeIndexText.move(-8,0);
            window.draw(nodeIndexText);

            vertexShape.setPosition((n->index()*60+10)*scaling,210*scaling);
            if (n->degree()==0)
            {
                vertexShape.setFillColor(sf::Color(240,240,240,255)); // black background, to overdraw fixed matching edges
                vertexShape.setOutlineColor(sf::Color(240,240,240,255));
                window.draw(vertexShape);
                vertexShape.setOutlineColor(sf::Color(0,0,0,alpha));
            }
            vertexShape.setFillColor(sf::Color(0,255,0,alpha)); // copy
            window.draw(vertexShape);
            nodeIndexString=to_string(n->index()+m_nodesA+m_nodesJ);
            nodeIndexText.setPosition((n->index()*60+16)*scaling,212*scaling);
            nodeIndexText.setString(nodeIndexString);
            nodeIndexText.setFillColor(sf::Color(0,0,0,alpha));
            if (n->index()+m_nodesA+m_nodesJ>9)
                nodeIndexText.move(-8,0);
            window.draw(nodeIndexText);

        }
        else if (enumNodeType[n]==NodeType::JOB)
        {
            vertexShape.setPosition(((n->index()-m_nodesA)*60+10)*scaling,90*scaling);
            if (n->degree()==0)
            {
                vertexShape.setFillColor(sf::Color(240,240,240,255)); // black background, to overdraw fixed matching edges
                vertexShape.setOutlineColor(sf::Color(240,240,240,255));
                window.draw(vertexShape);
                vertexShape.setOutlineColor(sf::Color(0,0,0,alpha));
            }
            vertexShape.setFillColor(sf::Color(0,200,0,alpha)); // job
            window.draw(vertexShape);
            nodeIndexString=to_string(n->index());
            nodeIndexText.setPosition(((n->index()-m_nodesA)*60+16)*scaling,92*scaling);
            nodeIndexText.setString(nodeIndexString);
            nodeIndexText.setFillColor(sf::Color(0,0,0,alpha));
            if (n->index()>9)
                nodeIndexText.move(-7,0);
            window.draw(nodeIndexText);

            vertexShape.setPosition(((n->index()-m_nodesA)*60+10)*scaling,150*scaling);
            if (n->degree()==0)
            {
                vertexShape.setFillColor(sf::Color(240,240,240,255)); // black background, to overdraw fixed matching edges
                vertexShape.setOutlineColor(sf::Color(240,240,240,255));
                window.draw(vertexShape);
                vertexShape.setOutlineColor(sf::Color(0,0,0,alpha));
            }
            vertexShape.setFillColor(sf::Color(0,200,0,alpha)); // copy
            window.draw(vertexShape);
            nodeIndexString=to_string(n->index()+m_nodesA+m_nodesJ);
            nodeIndexText.setPosition(((n->index()-m_nodesA)*60+16)*scaling,152*scaling);
            nodeIndexText.setString(nodeIndexString);
            nodeIndexText.setFillColor(sf::Color(0,0,0,alpha));
            if (n->index()+m_nodesA+m_nodesJ>9)
                nodeIndexText.move(-8,0);
            window.draw(nodeIndexText);

        }
    }
    while (true)
    {
        for (unsigned i=0; i<m_EnumRenderWindow.size(); ++i)
        {
            if (!m_EnumRenderWindow[i]->isOpen())
                continue;
            sf::Event event;
            if (m_EnumRenderWindow[i]->pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                {
                    m_EnumRenderWindow[i]->close();
                    if (i==jobIndex)
                        return;
                }
                else if (event.type ==  sf::Event::KeyPressed)
                {
                    return;
                }
            }
            m_EnumRenderWindow[i]->display();
        }
    }
#endif
}

vector<edge>* MWM::enumNextMatching(unsigned jobIndex)
{
    if (m_EnumStoreMatchings && m_EnumStoredMatchings[jobIndex]) // stored enumerated matchings? return address of computed matching
    {
        //vector<vector<vector<edge>*>> m_EnumStoredMatching;
        vector<vector<edge>*>& storedMatching=m_EnumStoredMatching[jobIndex];
        unsigned &curStoredMatching=m_EnumCurStoredMatching[jobIndex];
        if (curStoredMatching < storedMatching.size()) // there is at least one more stored matching
        {
            return storedMatching[curStoredMatching++];
        }
        else // start new
        {
            curStoredMatching=0;
            return nullptr;
        }
    }

    Graph &G=m_EnumerationGraph[jobIndex];
    stack<StackOP>& enumStack=m_EnumStack[jobIndex];
    NodeArray<unsigned>& scc=m_EnumSCC[jobIndex];
    const NodeArray<node>& copy = m_EnumCopy[jobIndex];
    vector<node>& nonIsolatedVertex=m_EnumNonIsolatedVertices[jobIndex];
    unsigned& numberofNonIsolatedVertices = m_EnumNumberofNonIsolatedVertices[jobIndex];
    const NodeArray<NodeType>& nodeType=m_EnumNodeType[jobIndex];
    EdgeArray<bool>& matched=m_EnumMatched[jobIndex];
    vector<edge>& fixedMatching=m_EnumFixedMatching[jobIndex];

    //node n;
    edge e;

    while (!enumStack.empty())
    {
        EnumRestore& enumRestore=*(m_EnumRestore[jobIndex].back());
        StackOP stackElement=enumStack.top();
        enumStack.pop();
        if (stackElement==StackOP::ENUM)
        {
#ifndef N_MATCHING_DEBUG
            cout << "OP: ENUM"<<endl;
#endif
            // if G is empty, then we have computed a MWM and return it
            if (G.numberOfEdges()==0)
            {
#ifndef N_MATCHING_DEBUG
                cout << "Matching complete. ";
#endif
                //displayEnum(jobIndex);
                if (m_EnumStoreMatchings) // Store computed matchings?
                {
                    vector<edge>* newMatchingToStore= new vector<edge>;
                    newMatchingToStore->reserve(fixedMatching.size());
                    for(edge e:fixedMatching)
                    {
                        newMatchingToStore->push_back(e);
                    }
                    m_EnumStoredMatching[jobIndex].push_back(newMatchingToStore);
                }
                return &fixedMatching;
            }

            // else put recursive operations on stack
            enumStack.push(StackOP::UNDO);
            enumStack.push(StackOP::ENUM);
            enumStack.push(StackOP::PREP_M2);
            enumStack.push(StackOP::ENUM);
            enumStack.push(StackOP::PREP_M1);
        }

        else if (stackElement==StackOP::PREP_M1)
        {
#ifndef N_MATCHING_DEBUG
            cout << "OP: PREP_M1"<<endl;
#endif
            // do the following: compute strongly connected component, in which to find an alternating circle C; change matching with respect to C; determine E1,E2, Trim(G/E1);

            // compute strongly connected component Z with minimum value F(Z)=|E(Z)|-|V(Z)|+cc(G)
            unsigned numSCCs=m_EnumMaxSCCNumber[jobIndex]+1;
            vector<int> f(numSCCs,0);
            vector<bool> hasOriginalNodes(numSCCs,false);
            vector<bool> hasCopyNodes(numSCCs,false);
            NodeArray<short>& nodeVisited=m_EnumDFSstatus[jobIndex];
            EdgeArray<edge>& previusEdgeOnCircle=m_EnumPreviousEdgeOnCircle[jobIndex];
            forall_edges(e,G) // +1 for edges
            {
                ++f[scc[e->source()]];
                if (nodeType[e->source()]<=NodeType::JOB && nodeType[e->target()]<=NodeType::JOB)
                    ++f[scc[copy[e->source()]]];
                previusEdgeOnCircle[e]=nullptr; // for computing the circle
            }
            for (unsigned i=0; i<numberofNonIsolatedVertices; ++i)
//			for(node n:nonIsolatedVertex) // -1 for nodes
            {
                node n=nonIsolatedVertex[i];
                --f[scc[n]];
                if (nodeType[n] <= NodeType::JOB)
                    hasOriginalNodes[scc[n]]=true;
                else
                    hasCopyNodes[scc[n]]=true;
                nodeVisited[n]=false; // for computing the circle
            }
            size_t minfIndex=0; // find minimum
            for (size_t fIndex=1; fIndex<f.size(); ++fIndex)
            {
                if (f[fIndex]<f[minfIndex])
                    minfIndex=fIndex;
            }
            //cout << "MinfIndex="<<minfIndex<<" has value "<<f[minfIndex]<<endl;

            // compute circle C with at least one edge in original graph
            if (!hasOriginalNodes[minfIndex]) // only vertices in copy; then compute scc of same component in original graph
            {
                for (unsigned i=0; i<numberofNonIsolatedVertices; ++i)
                //for(node n:nonIsolatedVertex)
                {
                    node n=nonIsolatedVertex[i];
                    if (scc[n]==minfIndex)
                    {
                        minfIndex=scc[copy[n]];
                        break; // for(node n:nonIsolatedVertex)
                    }
                }
            }
            edge startedge=nullptr;
            if (!hasCopyNodes[minfIndex]) // only vertices in original graph. Then select any edge of this component
            {
                forall_edges(e,G)
                {
                    if (scc[e->source()] == minfIndex)
                    {
                        startedge=e;
                        break; // forall_edges(e,G)
                    }
                }
            }
            else // vertices in both parts. Then select edge between copy and original graph.
            {
                forall_edges(e,G)
                {
                    const node s=e->source();
                    const node t=e->target();
                    if (scc[s] == minfIndex && ((nodeType[s] == NodeType::AGENT && nodeType[t] == NodeType::AGENTCOPY) || (nodeType[s] == NodeType::JOBCOPY && nodeType[t] == NodeType::JOB)))
                    {
                        startedge=e;
                        break; // forall_edges(e,G)
                    }
                }
            }
            edge pathEdge=startedge,neighbourEdge; // traverse graph, until a node is visited twice or edge to copy is found. With that information a circle can be constructed
            node pathNode=startedge->source(),neighbourNode=nullptr;
            if (nodeType[pathNode] == NodeType::JOBCOPY)
                pathNode=startedge->target();
            nodeVisited[pathNode]=true;
            getOtherEdgeAndNodeWithDifferentMatchedValue(matched, pathNode, pathEdge, neighbourNode, neighbourEdge);
            while (!nodeVisited[neighbourNode] && nodeType[neighbourNode] <= NodeType::JOB) // stop, if next node is a copy node or an already visited node.
            {
                previusEdgeOnCircle[neighbourEdge]=pathEdge;
                nodeVisited[neighbourNode]=true;
                pathNode=neighbourNode;
                pathEdge=neighbourEdge;
                getOtherEdgeAndNodeWithDifferentMatchedValue(matched, pathNode, pathEdge, neighbourNode, neighbourEdge);
            }
            vector<edge> &circle=enumRestore.circle;
            if (nodeVisited[neighbourNode]) // last node already visited. Then build circle by iterating back until node is part of the edge again.
            {
                previusEdgeOnCircle[neighbourEdge]=pathEdge;
                do
                {
                    circle.push_back(neighbourEdge);
                    //cout << neighbourEdge;
                    neighbourEdge=previusEdgeOnCircle[neighbourEdge];
                }
                while (neighbourEdge->source() != neighbourNode && neighbourEdge->target() != neighbourNode); // neighbourNode is start of circle
                circle.push_back(neighbourEdge);
            }
            else  // last node is of graph copy. Then build circle by iterating back to starting edge
            {
                circle.push_back(neighbourEdge);
                circle.push_back(pathEdge);
                while (previusEdgeOnCircle[pathEdge] != nullptr)
                {
                    //cout << pathEdge;
                    pathEdge=previusEdgeOnCircle[pathEdge];
                    circle.push_back(pathEdge);
                }
            }
#ifndef N_MATCHING_DEBUG
            cout << "Circle: ";
            for (edge e:circle)
            {
                cout << e<< " ";
            }
            cout << endl;
#endif
            // find matched edge in circle
            edge e1=circle.front();
            if (!matched[e1])
                e1=circle[1];

            // Compute edge sets E1,E2 of Uno's Algorithm
            vector<edge> edgeSet1 (1,e1); // Edge set E1 of Uno's Algorithm
            node v=e1->source(); // arbitrarily node v of the edge in E1
            if (nodeType[v]>=AGENTCOPY)
                v=e1->target(); // make sure v is vertex of original graph. This simplifies computation of E2
            vector<edge>& edgeSet2 = enumRestore.edgeSet2;
            edgeSet2.reserve(v->degree()-1);
            adjEntry adj;
            forall_adj(adj,v) // construct E2
            {
                e=adj->theEdge();
                if (e1!=e) // edge is not the edge in E1
                {
                    edgeSet2.push_back(e);
                }
            }
            // todo: additional tests and possible new computation of E1,E2, as described in the paper.

#ifndef N_MATCHING_DEBUG
            cout << "v:"<<v<<" E1:"<<e1<<" E2:";
            for(edge e:edgeSet2)
                cout << e;
            cout <<endl;
#endif

            // apply alternating circle to current matching
            for (edge e:circle)
                matched[e]=!matched[e];
            //displayEnum(jobIndex);

            // remove edges in E1 and trim graph
            trim(jobIndex,&edgeSet1);
        }

        else if (stackElement==StackOP::PREP_M2)
        {
#ifndef N_MATCHING_DEBUG
            cout << "OP: PREP_M2"<<endl;
#endif
            // do the following: undo fixed matching, reverse trimming, undo alternating circle, Trim(G/E2)

            // undo fixed matching
#ifndef N_MATCHING_DEBUG
            cout << "undo "<<enumRestore.additionalFixedMatchingEdges<<" fixed matching edges: ";
#endif
            for (int i=0; i<enumRestore.additionalFixedMatchingEdges; ++i)
            {
#ifndef N_MATCHING_DEBUG
                cout << fixedMatching.back();
#endif
                fixedMatching.pop_back();
            }
            // reverse trimming
#ifndef N_MATCHING_DEBUG
            cout <<endl<< "del "<< enumRestore.shortcutedges.size()<<" shortcuts:";
#endif
            EdgeArray<vector<edge>*>& deg2ReplacedPath=m_EnumDeg2ReplacedPath[jobIndex];
            m_EnumEdgeID[jobIndex] -= static_cast<int>(enumRestore.shortcutedges.size());
            for (edge e:enumRestore.shortcutedges)
            {
#ifndef N_MATCHING_DEBUG
                cout << e;
#endif
                delete deg2ReplacedPath[e];
                deg2ReplacedPath[e]=nullptr;
                G.delEdge(e);
            }
            //G.restoreEdges(enumRestore.hesh);
            enumRestore.hes.restore();
            EdgeArray<bool>& deg2TrimmedEdges=m_EnumDeg2TrimmedEdges[jobIndex];
#ifndef N_MATCHING_DEBUG
            cout << endl << "undo "<<enumRestore.trimmedEdges.size()<<" trimmed Edges:";
#endif
            for (edge e:enumRestore.trimmedEdges)
            {
#ifndef N_MATCHING_DEBUG
                cout << e;
#endif
                deg2TrimmedEdges[e]=false;
            }
#ifndef N_MATCHING_DEBUG
            cout << endl<<"Restoring "<<enumRestore.removedIsolatedVertices << " isolated vertices:";
#endif
    /*		for (node n:enumRestore.isolatedVertices)
            {
#ifndef N_MATCHING_DEBUG
                cout << n<< " ";
#endif
                nonIsolatedVertex.push_back(n);
            }*/
            numberofNonIsolatedVertices += enumRestore.removedIsolatedVertices;
            delete m_EnumRestore[jobIndex].back(); // remove saved (and restored) data
            m_EnumRestore[jobIndex].pop_back();
            EnumRestore* newEnumRestore=m_EnumRestore[jobIndex].back();
            // undo alternating circle
#ifndef N_MATCHING_DEBUG
            cout <<endl<< "undo " << newEnumRestore->circle.size()<<" circle edges: ";
#endif
            for (edge e:newEnumRestore->circle)
            {
#ifndef N_MATCHING_DEBUG
                cout << e;
#endif
                matched[e] = !matched[e];
            }
#ifndef N_MATCHING_DEBUG
            cout <<endl;
#endif
            //displayEnum(jobIndex);

            // Trim (G/E2)
            trim(jobIndex, &(newEnumRestore->edgeSet2));

        }

        else // UNDO
        {
#ifndef N_MATCHING_DEBUG
            cout << "OP: UNDO"<<endl;
#endif
            // do the following: undo fixed matching, reverse trimming

            // undo fixed matching
#ifndef N_MATCHING_DEBUG
            cout << "undo "<<enumRestore.additionalFixedMatchingEdges<<" fixed matching edges: ";
#endif
            for (int i=0; i<enumRestore.additionalFixedMatchingEdges; ++i)
            {
#ifndef N_MATCHING_DEBUG
                cout << fixedMatching.back();
#endif
                fixedMatching.pop_back();
            }
            // reverse trimming
#ifndef N_MATCHING_DEBUG
            cout <<endl<< "del "<< enumRestore.shortcutedges.size()<<" shortcuts:";
#endif
            EdgeArray<vector<edge>*>& deg2ReplacedPath=m_EnumDeg2ReplacedPath[jobIndex];
            m_EnumEdgeID[jobIndex] -= static_cast<int>(enumRestore.shortcutedges.size());
            for (edge e:enumRestore.shortcutedges)
            {
#ifndef N_MATCHING_DEBUG
                cout << e;
#endif
                delete deg2ReplacedPath[e];
                deg2ReplacedPath[e]=nullptr;
                G.delEdge(e);
            }
            //G.restoreEdges(enumRestore.hesh);
            enumRestore.hes.restore();
            EdgeArray<bool>& deg2TrimmedEdges=m_EnumDeg2TrimmedEdges[jobIndex];
#ifndef N_MATCHING_DEBUG
            cout << endl << "undo "<<enumRestore.trimmedEdges.size()<<" trimmed Edges:";
#endif
            for (edge e:enumRestore.trimmedEdges)
            {
#ifndef N_MATCHING_DEBUG
                cout << e;
#endif
                deg2TrimmedEdges[e]=false;
            }
#ifndef N_MATCHING_DEBUG
            cout << endl<<"Restoring "<<enumRestore.removedIsolatedVertices << " isolated vertices:";
#endif
            /*for (node n:enumRestore.isolatedVertices)
            {
#ifndef N_MATCHING_DEBUG
                cout << n<< " ";
#endif
                nonIsolatedVertex.push_back(n);
            }*/
            numberofNonIsolatedVertices += enumRestore.removedIsolatedVertices;

            delete m_EnumRestore[jobIndex].back(); // remove saved (and restored) data
            m_EnumRestore[jobIndex].pop_back();

            //displayEnum(jobIndex);
        }
    }

    // enumeration finished. Start from beginning with next call of function
    m_EnumRestore[jobIndex].back()->circle.clear();
    m_EnumRestore[jobIndex].back()->edgeSet2.clear();
    computeSCCs(jobIndex);
    enumStack.push(StackOP::ENUM);
    if (m_EnumStoreMatchings)
        m_EnumStoredMatchings[jobIndex]=true;
    return nullptr;


    /* EnumMatchings(G,M)
         * ||G||=0 -> output "matching"
         * E1,E2 bestimmen:
         *   ZHK nach Uno bestimmen
         *   Falls ZHK nur aus G1-Knoten besteht: Suche beliebigen Kreis darin
         *   Falls ZHK nur aus G2-Knoten besteht: Nehme stattdessen die entsprechende ZHK in dem anderen Teilgraphen
         *   sonst:
         *    // (Falls in G1 keine Matchingkanten sind: wähle beliebige Kante e in G1, sowie die 2 Matchingkanten, die zu den Endpunkten von e inzident sind, sowie die Kopie von e in G2. Ist ein Kreis.) auslassen?
         *      Starte mit einer Kante, die von G2 zu G1 verläuft (muss existieren), Tiefensuche, bis bereits besuchten Knoten gefunden -> Kreis, oder bis Kante von G1 zu G2 gewählt.
         *        Dann Ergibt sich Kreis aus den bereits gewählten Kanten und den Kopien der Kanten aus G1 in G2.
         * Setze E1=beliebige Matchingkante des Kreises, v=beliebiger Knoten der Kante
         * Setze E2=alle anderen Kanten, die mit v verbunden sind
         * bestimme M' durch M und C
         * (Hier müsste noch E1,E2 nach Uno geprüft und eventuell neu bestimmt werden)
         *
         * G:=Trim(G\E1,M')
         * output all edges of M' not in G (i.e. store them in array)
         * EnumMatchings(G,Trim(M'))
         * output "delete" and all edges of M' not in G (i.e. remove from array)
         * Reverse Trimming of G
         *
         * same as above, with E2,M instead of E1,M'*/
}

