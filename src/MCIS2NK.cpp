#include <global.h>
#include <MCIS2.h>
#pragma warning(push, 0)
#include <ogdf/fileformats/GraphIO.h>
#pragma warning(pop)
#include <bbpmcsi.h>

#define CASE_NOT_COMPUTED 127

#ifdef _DEBUG
    #include <assert.h>
#endif

using namespace ogdf;

ProductNodeNeighbors::ProductNodeNeighbors()
{
    m_neighborsAndQuantity.resize(ProductNodeNeighborsMAXNEIGHBORS * 4); // memory to store neighbors of nG and nH and their quantity
}

void ProductNodeNeighbors::addNeighbor(ProductnodeNoExpand &target)
{
    // already overfull
    if (m_neighborship_overfull)
        return;

    // all neighbors of nG
    bool found = false;
    for (size_t neighbor = 0; neighbor < m_curNeighborsG ; ++neighbor)
    {
        if (m_neighborsAndQuantity[neighbor*4] == target.nodeG) // this neighbor was found previously
        {
            ++m_neighborsAndQuantity[neighbor*4 + 1];
            found = true;
            break;
        }
    }
    if (!found) // insert new neighbor
    {
        if (m_curNeighborsG < ProductNodeNeighborsMAXNEIGHBORS)
        {
            m_neighborsAndQuantity[m_curNeighborsG*4]=target.nodeG;
            m_neighborsAndQuantity[m_curNeighborsG*4+1]=1;
            ++m_curNeighborsG;
        }
        else
        {
            m_neighborship_overfull = true;
        }
    }

    // all neighbors of nH
    found = false;
    for (size_t neighbor = 0; neighbor < m_curNeighborsH ; ++neighbor)
    {
        if (m_neighborsAndQuantity[neighbor*4 + 2] == target.nodeH) // this neighbor was found previously
        {
            ++m_neighborsAndQuantity[neighbor*4 + 3];
            found = true;
            break;
        }
    }
    if (!found) // insert new neighbor
    {
        if (m_curNeighborsH < ProductNodeNeighborsMAXNEIGHBORS)
        {
            m_neighborsAndQuantity[m_curNeighborsH*4 + 2]=target.nodeH;
            m_neighborsAndQuantity[m_curNeighborsH*4 + 3]=1;
            ++m_curNeighborsH;
        }
        else
        {
            m_neighborship_overfull = true;
        }
    }
}


bool ProductNodeNeighbors::removeNeighbor(ProductnodeNoExpand &target)
{
    // if already removed, do not do it again, prevent queuing by returning false
    /*if (m_wasRemoved)
    {
        cerr << "not here???";
        exit(1);
        return false;
    }*/

    // was  overfull, therefore return as "still enough neighbors"
    if (m_neighborship_overfull)
        return false;

    bool hasSufficientNeighborsBefore=hasSufficientAmountOfNeighbors();
    int nodesRemoved=0;
    // remove neighbors
    for (size_t neighbor = 0; neighbor < m_curNeighborsG ; ++neighbor)
    {
        if (m_neighborsAndQuantity[neighbor*4] == target.nodeG) // found corresponding neighbor
        {
            if (--m_neighborsAndQuantity[neighbor*4 + 1] == 0) // decrement amount, and move last neighbor to this position, if there are more
            {
                --m_curNeighborsG;
                if (neighbor < m_curNeighborsG) // there are neighbors left after this position
                {
                    m_neighborsAndQuantity[neighbor*4] = m_neighborsAndQuantity[m_curNeighborsG*4];
                    m_neighborsAndQuantity[neighbor*4 + 1] = m_neighborsAndQuantity[m_curNeighborsG*4 + 1];
                }
            }
            ++nodesRemoved;
            break;
        }
    }
    for (size_t neighbor = 0; neighbor < m_curNeighborsH ; ++neighbor)
    {
        if (m_neighborsAndQuantity[neighbor*4 + 2] == target.nodeH) // found corresponding neighbor
        {
            if (--m_neighborsAndQuantity[neighbor*4 + 3] == 0) // decrement amount, and move last neighbor to this position, if there are more
            {
                --m_curNeighborsH;
                if (neighbor < m_curNeighborsH) // there are neighbors left after this position
                {
                    m_neighborsAndQuantity[neighbor*4 + 2] = m_neighborsAndQuantity[m_curNeighborsH*4 + 2];
                    m_neighborsAndQuantity[neighbor*4 + 3] = m_neighborsAndQuantity[m_curNeighborsH*4 + 3];
                }
            }
            ++nodesRemoved;
            break;
        }
    }

    if (nodesRemoved != 2)
    {
        cerr << "MCIS2NK: Tried to remove non existing nodes. This should not happen. Please contact program author.\n";
        exit(1);
    }
    if (hasSufficientNeighborsBefore && !hasSufficientAmountOfNeighbors())
    {
        return true;
    }
    return false;
}

bool ProductNodeNeighbors::hasSufficientAmountOfNeighbors()
{
    //if (m_wasRemoved)
        //return false;
    if (m_neighborship_overfull || (m_curNeighborsG >= 2 && m_curNeighborsH >= 2))
        return true;
    return false;
}







MaxComPart2Tree::~MaxComPart2Tree()
{
    if (m_outerplanar)
    {
        edge e,f;
        forall_edges(e, m_G) {
            forall_edges(f, m_H) {
                delete[] m_D[e][f];
            }
        }
        delete m_map_snode;
    }
    else
    {
        delete[] m_bi_low;
        delete[] m_bi_num;
        delete m_adjH;
        delete m_adjG;
        delete[] m_productNodeWeight;
        delete[] m_productNodeNeighbors;
        delete[] m_nodesH;
        delete[] m_nodesG;
    }
}

// add a node mapping to the MCS
void MaxComPart2Tree::addNodeMappingToMCS(node vG, node vH)
{
    node n=m_MCS.newNode();
    m_vGtovH[vG]=vH;
    m_vGtoMCS[vG]=n;
    m_vMCStovG[n]=vG;
}

// add a edge mapping to the MCS
void MaxComPart2Tree::addEdgeMappingToMCS(edge eG, edge eH)
{
    node vMCS=m_vGtoMCS[eG->source()];
    node uMCS=m_vGtoMCS[eG->target()];
    edge e=m_MCS.newEdge(vMCS, uMCS);
    m_eGtoeH[eG]=eH;
    //m_eGtoMCS[eG]=e;
    m_eMCStoeG[e]=eG;
}

void MaxComPart2Tree::applyWeightsToMCS()
{
    node nMCS;
    edge eMCS,eG,eH;
    size_t caseCache;
    // default to WEIGHT_NOT_COMPATIBLE, overwrite from the computed split isomorphisms
    forall_edges(eMCS, m_MCS) {
        eG = m_eMCStoeG[eMCS];
        eH = m_eGtoeH[eG];
        caseCache=(m_detectCaseCache[eG][eH]=detectCase(eG,eH));
#ifndef LAWECSU_NDEBUG
        if (m_D[eG][eH][caseCache] != WEIGHT_UNDEF)
            cerr << "MCIS2 m_D calculated twice!"<<endl;
#endif
        m_D[eG][eH][caseCache] = WEIGHT_NOT_COMPATIBLE;
    }

    // edge weights
    forall_edges(eMCS,m_MCS)
        m_eMCStoWeight[eMCS]= m_bbp->getWeight(m_eGtoaeG[m_eMCStoeG[eMCS]],m_eHtoaeH[m_eGtoeH[m_eMCStoeG[eMCS]]]);

    // node weights
    forall_nodes(nMCS,m_MCS)
    {
        if (m_vMCStovG[nMCS]==m_xG) // check for excluded vertex
            m_vMCStoWeight[nMCS] = WEIGHT_NOT_COMPATIBLE;
        else if (m_vMCStovG[nMCS]==m_vG) // along vertex of fixed mapping no expansion into other B-nodes allowed
        {
            m_vMCStoWeight[nMCS]= m_bbp->getWeight(m_vGtoaG[m_vMCStovG[nMCS]],m_vHtoaH[m_vGtovH[m_vMCStovG[nMCS]]]);
            m_vMCStoExpand[nMCS]= NO_EXPANSION;
        }
        else
            m_vMCStoWeight[nMCS]= m_bbp->getWeightPlusCVMatching(m_bbp->getBCG().original(m_vGtoaG[m_vMCStovG[nMCS]]),m_vHtoaH[m_vGtovH[m_vMCStovG[nMCS]]],m_vMCStoExpand[nMCS]);

        if (m_vMCStoWeight[nMCS] == WEIGHT_NOT_COMPATIBLE) // incident edges of incompatible vertices are also incompatible
        {
            for(adjEntry adj : nMCS->adjEntries)
            {
                edge eMCS=adj->theEdge();
                //forall_adj_edges(eMCS,nMCS) //
                m_eMCStoWeight[eMCS] = WEIGHT_NOT_COMPATIBLE;
            }
        }
    }
}

void MaxComPart2Tree::computeSplitIsomorphisms(forward_list<vector<Mapping>>& maximumMappingList, wType &maxWeight)
{
    const bool &enumerate=m_enumerate;
    // new Version based on SP-Tree
    StaticPlanarSPQRTree spSPQRTree(m_MCS);
    const Graph &sptree=spSPQRTree.tree();

    // For each not compatible real edge in a P-Node mark the P-Node and adjacent S-Nodes as not compatible and for each not compatible real edge in a S-Node mark the S-Node as not compatible
    NodeArray<bool> compatibleSPNode(sptree,true);
    edge eMCS;
    node skelNode;
    adjEntry adj;
    forall_edges(eMCS,m_MCS)
        if (m_eMCStoWeight[eMCS] == WEIGHT_NOT_COMPATIBLE)
        {
            skelNode=spSPQRTree.skeletonOfReal(eMCS).treeNode();
            compatibleSPNode[skelNode]=false;
            if (spSPQRTree.typeOf(skelNode)==SPQRTree::NodeType::PNode)
            {
                forall_adj(adj,skelNode)
                    compatibleSPNode[adj->twinNode()]=false;
            }
        }

    // Determine for each node of the SP-Tree the connected component regarding compatibleSPNode. These are the 2-connected parts of the MCS
    int component=0,numNodes;
    NodeArray<int> connectedComponent(sptree,-1);
    std::queue<node> conQ;
    std::vector<node> componentsOrdered;
    node nSP,nQueue,twin;
    forall_nodes(nSP,sptree)
    {
        if (compatibleSPNode[nSP] && connectedComponent[nSP]==-1)
        {
            conQ.push(nSP);
            numNodes=1;
            while (!conQ.empty())
            {
                nQueue=conQ.front();
                connectedComponent[nQueue]=component;
                componentsOrdered.push_back(nQueue);
                forall_adj(adj,nQueue)
                {
                    twin=adj->twinNode();
                    if (compatibleSPNode[twin] && connectedComponent[twin]==-1)
                    {
                        conQ.push(twin);
                        ++numNodes;
                    }
                }
                conQ.pop();
            }
            if (numNodes==1 && spSPQRTree.typeOf(nSP)==SPQRTree::NodeType::PNode) // single P-Nodes are not 2-connected
            {
                connectedComponent[nSP]=-1;
                compatibleSPNode[nSP]=false;
                componentsOrdered.pop_back();
            }
            else
                ++component;
        }
    }

    if (componentsOrdered.empty()) // no 2-connected component remaining
        return;

    // For each component sum the edge and node weights. Also check if possible required vertex is part of the component. Then store total weight in m_D
    wType weight=0;
    bool m_vG_Contained=(m_vG==nullptr?true:false);
    NodeArray<int> nodeComponent(m_MCS,component);
    std::vector<edge> componentEdges;
    std::vector<node> componentNodes;
    --component;
    edge eSkel,eG,eH;
    node nMCS,nG,nH;

    while (!componentsOrdered.empty())
    {
        nSP=componentsOrdered.back();
        const Skeleton &skel=spSPQRTree.skeleton(nSP);
        const Graph &skelGraph=skel.getGraph();
        forall_edges(eSkel,skelGraph)
            if (!skel.isVirtual(eSkel))
            {
                eMCS=skel.realEdge(eSkel);
                componentEdges.push_back(eMCS);
                weight+=m_eMCStoWeight[eMCS];
                nMCS=eMCS->source();
                if (nodeComponent[nMCS]>component) // node weight not yet added
                {
                    nodeComponent[nMCS]=component;
                    componentNodes.push_back(nMCS);
                    weight+=m_vMCStoWeight[nMCS];
                    if (m_vMCStovG[nMCS]==m_vG)
                        m_vG_Contained=true;
                }
                nMCS=eMCS->target();
                if (nodeComponent[nMCS]>component)
                {
                    nodeComponent[nMCS]=component;
                    componentNodes.push_back(nMCS);
                    weight+=m_vMCStoWeight[nMCS];
                    if (m_vMCStovG[nMCS]==m_vG)
                        m_vG_Contained=true;
                }
            }
        componentsOrdered.pop_back();
        if (componentsOrdered.empty() || connectedComponent[componentsOrdered.back()]<component) // next node in Q is of another component or Q empty
        {
            if (m_vG_Contained) // store computed weight in m_D
            {
                while (!componentEdges.empty())
                {
                    eMCS=componentEdges.back();
                    eG = m_eMCStoeG[eMCS];
                    eH = m_eGtoeH[eG];
                    m_D[eG][eH][(size_t)m_detectCaseCache[eG][eH]] = weight;
                    componentEdges.pop_back();
                }

                // add new solution, if it is large enough
                if (weight > maxWeight || (enumerate && weight == maxWeight))
                {
                    if (weight > maxWeight) // delete possible previous smaller solution(s)
                    {
                        maximumMappingList.clear();
                        maxWeight=weight;
                    }
                    maximumMappingList.push_front(vector<Mapping>());
                    while (!componentNodes.empty())
                    {
                        nMCS=componentNodes.back();
                        nG = m_vMCStovG[nMCS];
                        nH = m_vGtovH[nG];
                        if (m_vG!=nG) // do not store given vertex
                            maximumMappingList.front().push_back(Mapping(m_bbp->getBCG().original(m_vGtoaG[nG]),m_vHtoaH[nH],m_vMCStoExpand[nMCS]));
                        componentNodes.pop_back();
                    }
                }
            }
            else
            {
                componentEdges.clear();
                componentNodes.clear();
            }
            --component;
            weight=0;
            m_vG_Contained=(m_vG==nullptr?true:false);
        }
    }
}

//void MaxComPart2Tree::c_Clique(vector<Productnode> &P, vector<Productnode> &Q, vector<Productnode> &X, vector<Productnode> &Y, bool biconnected)
void MaxComPart2Tree::c_Clique(vector<ProductnodeNoExpand> &P, vector<ProductnodeNoExpand> &Q, bool biconnected, wType curWeight)
{
    // maximum time check
    if (m_maxCliqueTime>0 && (clock()-m_cliqueTimeStart) / (  CLOCKS_PER_SEC / 1000 ) > m_maxCliqueTime)
        return;
++m_cliqueRecursions;
    // number of adjacent vertices; modify weight
    int numAdjacency=0;
    curWeight += productNodeWeight(m_R.back());
    for (size_t i=0; i< m_R.size(); ++i)
    {
        if (getEdgeType(m_R.back(),m_R[i])==C_EDGE)
        {
            ++numAdjacency;
            curWeight += productEdgeWeight(m_R.back(),m_R[i]);
        }
    }

    // biconnect status
    if (m_R.size()<3) // need at least 3 vertices
        biconnected=false;
    else
    {
        // count with how many vertices the new one is connected
        if (numAdjacency == 1) // biconnectivity lost as new vertex is single connected
            biconnected=false;
        else // connected to at least 2 vertices
        {
            if (!biconnected) // if it was biconnected before, this will not change. Else recomputation required
            {
                m_bi_count=0;
                for (size_t i=0; i<m_R.size(); ++i)
                {
                    m_bi_num[i]=-1;
                }
                biconnected=!hasAC(0,-1);
            }
        }
    }

    // report clique, if biconnected
    if (biconnected)
        ReportClique(curWeight);

    // Zeile 1 bis 2
/*	if (P.empty())
    {
        if (X.empty())
            ;//ReportClique(biconnected);
        return;
    }*/

    vector<ProductnodeNoExpand> Q_new, P_new;//, X_new, Y_new;
    // Zeile 4 bis 5, Index, ab dem P beginnt.
    for (size_t P_start=0; P_start<P.size();)
    {
        // ui auslesen
        ProductnodeNoExpand ui=P[P_start];

        // Zeile 6
        ++P_start;

        // Zeile 7
        m_R.push_back(ui);

        // remove vertices, which do not have at least two connections through c_edges
        stack<ProductnodeNoExpand> sInsufficientNeighbors;
        for (size_t i=P_start; i<P.size(); ++i)
            if (getEdgeType(ui,P[i])==NO_EDGE)
                removeNeighborship(P[i],&sInsufficientNeighbors);

        for (size_t i=0; i<Q.size(); ++i)
            if (getEdgeType(ui,Q[i])==NO_EDGE)
                removeNeighborship(Q[i],&sInsufficientNeighbors);

        // lines 8,9
        Q_new.clear();
        P_new.clear();
        for (size_t i=P_start; i<P.size(); ++i)
            if (!m_productNodeNeighbors[P[i].nodeG * m_sizeH + P[i].nodeH].getwasRemoved()) // vertex not removed?
                if (getEdgeType(ui,P[i])!=NO_EDGE)
                    P_new.push_back(P[i]);

        for (size_t i=0; i<Q.size(); ++i)
            if (!m_productNodeNeighbors[Q[i].nodeG * m_sizeH + Q[i].nodeH].getwasRemoved()) // vertex not removed?
            {
                if (getEdgeType(ui,Q[i])==D_EDGE)
                    Q_new.push_back(Q[i]);
                else if (getEdgeType(ui,Q[i])==C_EDGE)
                    P_new.push_back(Q[i]);
            }


    /*	// Zeile 8
        Q_new.clear();
        for (size_t i=0; i<Q.size(); ++i)
            if (getEdgeType(ui,Q[i])==D_EDGE)
                Q_new.push_back(Q[i]);

        // Zeile 9
        P_new.clear();
        for (size_t i=P_start; i<P.size(); ++i)
            if (getEdgeType(ui,P[i])!=NO_EDGE)
                P_new.push_back(P[i]);
        for (size_t i=0; i<Q.size(); ++i)
            if (getEdgeType(ui,Q[i])==C_EDGE)
                P_new.push_back(Q[i]);*/

        // Zeile 10
        //Y_new.clear();
        //for (size_t i=0; i<Y.size(); ++i)
        //	if (getEdgeType(ui,Y[i])==D_EDGE)
        //		Y_new.push_back(Y[i]);

        // Zeile 11
        /*X_new.clear();
        for (size_t i=0; i<X.size(); ++i)
            if (getEdgeType(ui,X[i])!=NO_EDGE)
                X_new.push_back(X[i]);
        for (size_t i=0; i<Y.size(); ++i)
            if (getEdgeType(ui,Y[i])==C_EDGE)
                X_new.push_back(Y[i]);*/

        // Zeile 12
        // c_Clique(P_new, Q_new, X_new, Y_new);
        c_Clique(P_new, Q_new,biconnected,curWeight);

        // restore vertices, which have been excluded through NO_EDGE
        while (!sInsufficientNeighbors.empty())
        {
            addNeighborship(sInsufficientNeighbors.top());
            sInsufficientNeighbors.pop();
        }

        // letztes Element aus R wieder entfernen, um Originalzustand herzustellen. Erspart Kopie von R.
        m_R.pop_back();

        // Zeile 13
    //	X.push_back(ui);
    }
}

wType MaxComPart2Tree::productNodeWeight(ProductnodeNoExpand &n)
{
    return m_productNodeWeight[n.nodeG * m_sizeH + n.nodeH];
    /*if (m_nodesG[n.nodeG] == m_vG) // along vertex of fixed mapping no expansion into other B-nodes allowed
    {
        return m_bbp->getWeight(m_vGtoaG[m_nodesG[n.nodeG]],m_vHtoaH[m_nodesH[n.nodeH]]);
    }
    else
    {
        return m_bbp->getWeightPlusCVMatching(m_vGtoaG[m_nodesG[n.nodeG]],m_vHtoaH[m_nodesH[n.nodeH]],n.expand);
    }*/
}

wType MaxComPart2Tree::productEdgeWeight(ProductnodeNoExpand &from, ProductnodeNoExpand &to)
{
    int indexG=max(from.nodeG,to.nodeG) * (max(from.nodeG,to.nodeG)+1) / 2 + min(from.nodeG,to.nodeG);
    int indexH=max(from.nodeH,to.nodeH) * (max(from.nodeH,to.nodeH)+1) / 2 + min(from.nodeH,to.nodeH);
    return m_bbp->getWeight(m_eGtoaeG[(*m_adjG)[indexG]],m_eHtoaeH[(*m_adjH)[indexH]]);
}

bool MaxComPart2Tree::MaxComPart2Tree::hasAC(int u, int parent)
{
    // Count of children in DFS Tree
    int firstChild = -1;

    // Initialize discovery time and low value
    m_bi_num[u] = m_bi_low[u] = m_bi_count++;

    for (int i=0; i<(int) m_R.size(); ++i)
    {
        if (getEdgeType(m_R[u],m_R[i])==C_EDGE)
        {
            //cout << "("<<u<<","<<i<<") low["<<u<<"]="<<m_bi_low[u]<<" low["<<i<<"]="<<m_bi_low[i]<<" \n";
            if (m_bi_num[i] == -1)
            {
                if (firstChild == -1)
                    firstChild = i;

                if (hasAC(i,u))
                    return true;

                if (m_bi_low[i] >= m_bi_num[u] && (parent != -1 ||  i != firstChild))
                    return true;

                m_bi_low[u]=min(m_bi_low[u],m_bi_low[i]);
            }
            else if (i != parent)
                m_bi_low[u] = min(m_bi_low[u],m_bi_num[i]);
        }
    }
    return false;
}

void MaxComPart2Tree::ReportClique(wType weight)
{
//	cout << "report clique\n";

    // check for 2-connectivity
    /*bool biconnected=false;
    if (m_R.size() > 2)
    {
        m_bi_count=0;
        for (size_t i=0; i<m_R.size(); ++i)
        {
            m_bi_num[i]=-1;
        }

        biconnected=!hasAC(0,-1);
    }

    // only keep biconnected cliques*/
    /*if (!biconn)
    {
        //cout << "not biconnected\n";
        return;
    }*/
    //cout << "biconnected\n";
    /*if (biconn != biconnected)
        cout << "Biconn != biconnected\n";*/


    // add egde and vertex weight and set pointers to maximumMappingList[aH] and maxWeight[aH];
    //wType weight=0;
    //bool vG_Contained=(m_vG==nullptr?true:false);
    for (size_t i=0; i<m_R.size(); ++i)
    {
        if (m_nodesG[m_R[i].nodeG]==m_vG)
        {
            //vG_Contained=true;
            node aH=m_vHtoaH[m_nodesH[m_R[i].nodeH]];
            m_maximumMappingListPtr=&(*m_maximumMappingArrayListPtr)[aH];
            m_maxWeightNonOuterPtr=&(*m_maxWeightArrayNonOuterPtr)[aH];
            break;
        }
        /*weight += productNodeWeight(m_R[i]);
        for (size_t j=i+1; j<m_R.size(); ++j)
        {
            if (getEdgeType(m_R[i],m_R[j])==C_EDGE)
            {
                weight += productEdgeWeight(m_R[i],m_R[j]);
            }
        }*/
    }
    //if (weight != cliqueWeight)
    //	cout << " weight differs!!!!\n";

    //if (vG_Contained)
    //{
        if (weight>*m_maxWeightNonOuterPtr || (m_enumerate && weight == *m_maxWeightNonOuterPtr))
        {
            if (weight > *m_maxWeightNonOuterPtr) // delete possible previous smaller solution(s)
            {
                m_maximumMappingListPtr->clear();
                *m_maxWeightNonOuterPtr=weight;
            }
            m_maximumMappingListPtr->push_front(vector<Mapping>());
            node nG,nH;
            for (size_t i=0; i<m_R.size(); ++i)
            {
                nG=m_nodesG[m_R[i].nodeG];
                nH=m_nodesH[m_R[i].nodeH];
                if (m_vG!=nG) // do not store given vertex
                    m_maximumMappingListPtr->front().push_back(Mapping(m_bbp->getBCG().original(m_vGtoaG[nG]),m_vHtoaH[nH],m_productNodeExpand[m_R[i].nodeG * m_sizeH + m_R[i].nodeH]));
            }
        }
/*		for (size_t i=0; i<R.size(); ++i)
            cout << R[i];
        cout << " Total weight: "<<weight;
        cout << endl;*/
    //}
    //else
        //cout << "vg not contained\n";
}

int MaxComPart2Tree::getEdgeType(ProductnodeNoExpand from, ProductnodeNoExpand to)
{
    if (from.nodeG==to.nodeG || from.nodeH==to.nodeH)
        return NO_EDGE;
    int indexG=max(from.nodeG,to.nodeG) * (max(from.nodeG,to.nodeG)+1) / 2 + min(from.nodeG,to.nodeG);
    int indexH=max(from.nodeH,to.nodeH) * (max(from.nodeH,to.nodeH)+1) / 2 + min(from.nodeH,to.nodeH);
    if ((*m_adjG)[indexG] != nullptr && (*m_adjH)[indexH] != nullptr)
    {
        //int indexG=max(from.nodeG,to.nodeG) * (max(from.nodeG,to.nodeG)+1) / 2 + min(from.nodeG,to.nodeG);
        //int indexH=max(from.nodeH,to.nodeH) * (max(from.nodeH,to.nodeH)+1) / 2 + min(from.nodeH,to.nodeH);
        //return m_bbp->getWeight(m_eGtoaeG[(*m_adjG)[indexG]],m_eHtoaeH[(*m_adjH)[indexH]]);

        if (m_bbp->getWeight(m_eGtoaeG[(*m_adjG)[indexG]],m_eHtoaeH[(*m_adjH)[indexH]]) != WEIGHT_NOT_COMPATIBLE)
        //if (productEdgeWeight(from, to) != WEIGHT_NOT_COMPATIBLE)
            return C_EDGE;
        else
            return NO_EDGE;
    }
    if ((*m_adjG)[indexG] == nullptr && (*m_adjH)[indexH] == nullptr)
        return D_EDGE;
//	if ((*adjG)[from.nodeG*sizeG+to.nodeG] && (*adjH)[from.nodeH*sizeH+to.nodeH])
//		return C_EDGE;
//	if (!(*adjG)[from.nodeG*sizeG+to.nodeG] && !(*adjH)[from.nodeH*sizeH+to.nodeH])
//		return D_EDGE;
    return NO_EDGE;
}

// initialize MCS computation for outerplanar and not outerplanar biconnected graphs
void MaxComPart2Tree::init()
{
    // non outerplanar case first
    if (!m_outerplanar)
    {
        m_sizeG=m_G.numberOfNodes();
        m_sizeH=m_H.numberOfNodes();
        m_nodesProductGraph=m_sizeG*m_sizeH;

        // Make nodes accessible via index
        node n;
        m_nodesG = new node[m_sizeG];
        m_nodesH = new node[m_sizeH];
        forall_nodes(n,m_G)
        {
            m_nodesG[n->index()]=n;
        }
        forall_nodes(n,m_H)
        {
            m_nodesH[n->index()]=n;
        }

        // compute weights
        m_productNodeWeight = new wType[m_nodesProductGraph];
        m_productNodeExpand.resize(m_nodesProductGraph, NO_EXPANSION);
        int index=0;
        bool expand;

        for (int nG = 0; nG < m_sizeG; ++nG)
        {
            for (int nH = 0; nH < m_sizeH; ++nH, ++index)
            {
                if (m_nodesG[nG] == m_vG)
                {
                    m_productNodeWeight[index] = m_bbp->getWeight(m_vGtoaG[m_nodesG[nG]],m_vHtoaH[m_nodesH[nH]]);
                }
                else if (m_nodesG[nG] == m_xG)
                {
                    m_productNodeWeight[index] = WEIGHT_NOT_COMPATIBLE;
                }
                else
                {
                    m_productNodeWeight[index] = m_bbp->getWeightPlusCVMatching(m_bbp->getBCG().original(m_vGtoaG[m_nodesG[nG]]),m_vHtoaH[m_nodesH[nH]],expand);
                    m_productNodeExpand[index]=expand;
                }
            }
        }

        // neighbors and occurences of them
        m_productNodeNeighbors = new ProductNodeNeighbors[m_nodesProductGraph];


        // Use lower left triangular matrix, as graphs are undirected
        //adjG= new vector<bool>(sizeG*sizeG,false);
        //adjH= new vector<bool>(sizeH*sizeH,false);
        //adjG= new vector<bool>(sizeG * (sizeG+1) / 2,false);
        //adjH= new vector<bool>(sizeH * (sizeH+1) / 2,false);
        m_adjG= new vector<edge>(m_sizeG * (m_sizeG+1) / 2, nullptr);
        m_adjH= new vector<edge>(m_sizeH * (m_sizeH+1) / 2, nullptr);
        edge e;

        forall_edges(e,m_G)
        {
            int si=e->source()->index();
            int ti=e->target()->index();
            (*m_adjG)[max(si,ti) * (max(si,ti)+1) / 2 + min(si,ti)]=e;
            //(*adjG)[e->target()->index()*sizeG+e->source()->index()] = (*adjG)[e->source()->index()*sizeG+e->target()->index()] = true;
        //	if (e->target()->index()>=sizeG || e->source()->index()>=sizeG)
        //		cerr << "index out of range in init\n";
        }
        forall_edges(e,m_H)
        {
            int si=e->source()->index();
            int ti=e->target()->index();
            (*m_adjH)[max(si,ti) * (max(si,ti)+1) / 2 + min(si,ti)]=e;
            //(*adjH)[e->target()->index()*sizeH+e->source()->index()] = (*adjH)[e->source()->index()*sizeH+e->target()->index()] = true;
        //	if (e->target()->index()>=sizeH || e->source()->index()>=sizeH)
        //		cerr << "index out of range in init\n";
        }
        //cout << "Non outerplanar block, Productgraphsize="<<m_nodesProductGraph<<"  sizeG*sizeG="<< sizeG*sizeG << "adjG.size()="<<adjG->size()<< "  sizeH*sizeH="<< sizeH*sizeH << "adjH.size()="<<adjH->size()<<endl;

        // checking biconnectivity
        m_bi_num = new int[m_nodesProductGraph];
        m_bi_low = new int[m_nodesProductGraph];

        // precompute vertex weights

        return;
    }
    // In case of 2-connected partial two, the SP-Tree equals the SPQR-Tree. (computed before)
    //m_SPG = new StaticPlanarSPQRTree(m_G);
    //m_SPH = new StaticPlanarSPQRTree(m_H);

    m_map_snode = new NodeBijection(m_SPG->tree(), m_SPH->tree());
    m_vGtovH.init(m_G);
    m_vGtoMCS.init(m_G);
    m_vMCStovG.init(m_MCS);

    m_eGtoeH.init(m_G,NULL);
    //m_eGtoMCS.init(m_G,NULL);
    m_eMCStoeG.init(m_MCS,NULL);

    m_vMCStoWeight.init(m_MCS);
    m_eMCStoWeight.init(m_MCS);
    m_vMCStoExpand.init(m_MCS);

    edge e,f;
    m_detectCaseCache.init(m_G);
    forall_edges(e, m_G)
        m_detectCaseCache[e].init(m_H,CASE_NOT_COMPUTED);

    // m_D.init(m_G,m_H);
    m_D.init(m_G);
    forall_edges(e, m_G) {
        m_D[e].init(m_H);
        forall_edges(f, m_H) {
            m_D[e][f] = new wType[4]{WEIGHT_UNDEF,WEIGHT_UNDEF,WEIGHT_UNDEF,WEIGHT_UNDEF};
        }
    }
}

int MaxComPart2Tree::detectCase(edge eG, edge eH) {
    bool lessIndexG = (eG->source()->index() < eG->target()->index());
    bool lessIndexH = (m_vGtovH[eG->source()]->index() < m_vGtovH[eG->target()]->index());
    char orientation = (lessIndexG == lessIndexH) ? 0 : 1;

    // at least one adjacent face must be mapped
    // we associated the non-mapped faces onto each other if required
    int idFaceG1 = -1, idFaceG2 = -1, idFaceH1 = -1, idFaceH2 = -1;

    const Skeleton &skelG = m_SPG->skeletonOfReal(eG);
    node nG = skelG.treeNode();
    if (m_SPG->typeOf(nG) == SPQRTree::NodeType::PNode) {
        edge f;
        for(adjEntry adj : nG->adjEntries)
        {
            f=adj->theEdge();
        //forall_adj_edges(f,nG) {
            node adjSNode = f->opposite(nG);
            if (idFaceG1 == -1 && m_map_snode->map(adjSNode) != NULL) {
                idFaceG1 = adjSNode->index();
            } else {
                idFaceG2 = adjSNode->index();
            }
        }
    } else {
        idFaceG1 = nG->index();
        idFaceG2 = -1;
    }
    const Skeleton &skelH = m_SPH->skeletonOfReal(eH);
    node nH = skelH.treeNode();
    if (m_SPH->typeOf(nH) == SPQRTree::NodeType::PNode) {
        edge f;
        for(adjEntry adj : nH->adjEntries)
        {
            f=adj->theEdge();
        //forall_adj_edges(f,nH) {
            node adjSNode = f->opposite(nH);
            if (m_map_snode->invMap(adjSNode) != NULL &&
                m_map_snode->invMap(adjSNode)->index() == idFaceG1) {
                idFaceH1 = adjSNode->index();
            } else {
                idFaceH2 = adjSNode->index();
            }
        }
    } else {
        idFaceH1 = nH->index();
        idFaceH2 = -1;
    }

    lessIndexG = idFaceG1 < idFaceG2;
    lessIndexH = idFaceH1 < idFaceH2;
    if (lessIndexG == lessIndexH) {
        return orientation;
    } else {
        return orientation + 2;
    }
}

void MaxComPart2Tree::deleteCurrentSolution()
{
    m_map_snode->clear();
    m_MCS.clear();
}

void MaxComPart2Tree::addNeighborship(ProductnodeNoExpand& source)
{
    node nodeG = m_nodesG[source.nodeG];
    node nodeH = m_nodesH[source.nodeH];
    adjEntry adjnodeG,adjnodeH;
    node nodeTargetG,nodeTargetH; // target nodes
    ProductnodeNoExpand target;

    forall_adj(adjnodeG,nodeG) // iterate over neighbors of nodeG
    {
        nodeTargetG = adjnodeG->twinNode();
        target.nodeG = nodeTargetG->index();

        //if (nodeTargetG != m_xG) // check if target node in G is not excluded
        forall_adj(adjnodeH,nodeH) // iterate over neighbors of nodeH
        {
            nodeTargetH = adjnodeH->twinNode();
            target.nodeH = nodeTargetH->index();
            if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE) // (source and) target compatible
            {
                if (getEdgeType(source,target)==C_EDGE) // connected edge between both product nodes, therefore increment neighborship
                {
                    //m_productNodeNeighbors[source.nodeG * m_sizeH + source.nodeH].addNeighbor(target);
                    m_productNodeNeighbors[target.nodeG * m_sizeH + target.nodeH].addNeighbor(source);
                }
            }
        }
    }
    m_productNodeNeighbors[source.nodeG * m_sizeH + source.nodeH].setwasRemoved(false);
}

void MaxComPart2Tree::removeNeighborship(ProductnodeNoExpand source, stack<ProductnodeNoExpand>* sInsufficientNeighbors)
{
    int index=source.nodeG * m_sizeH + source.nodeH;
    if (m_productNodeNeighbors[index].getwasRemoved())
        return;

    m_productNodeNeighbors[index].setwasRemoved(true);
    if (sInsufficientNeighbors != nullptr)
        sInsufficientNeighbors->push(source);

    node nodeG = m_nodesG[source.nodeG];
    node nodeH = m_nodesH[source.nodeH];
    adjEntry adjnodeG,adjnodeH;
    node nodeTargetG,nodeTargetH; // target nodes
    ProductnodeNoExpand target;
    // iterate over neighbours of source

    forall_adj(adjnodeG,nodeG) // iterate over neighbors of nodeG
    {
        nodeTargetG = adjnodeG->twinNode();
        target.nodeG = nodeTargetG->index();

        //if (nodeTargetG != m_xG) // check if target node in G is not excluded
        forall_adj(adjnodeH,nodeH) // iterate over neighbors of nodeH
        {
            nodeTargetH = adjnodeH->twinNode();
            target.nodeH = nodeTargetH->index();
            if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE) // (source and) target compatible
            {
                if (getEdgeType(source,target)==C_EDGE) // connected edge between both product nodes, therefore decrement neighborship
                {
                    if (m_productNodeNeighbors[target.nodeG * m_sizeH + target.nodeH].removeNeighbor(source)) // queue vertices, which got too less neighbors
                    {
                        removeNeighborship(target, sInsufficientNeighbors);
                    }
                }
            }
        }
    }
}

void MaxComPart2Tree::eliminateNodesOfDegreeOneAndZeroInProductGraph()
{
    ProductnodeNoExpand source,target;
    int &nG=source.nodeG;
    int &nH=source.nodeH;

    // iterate over all product nodes pn and determine their c-connected neighbors and quantities
    for (nG = 0; nG < m_sizeG; ++nG)
    {
        //if (m_nodesG[nG] != m_xG) // omit excluded vertex
        for (nH = 0; nH < m_sizeH; ++nH)
        {
            if (productNodeWeight(source) != WEIGHT_NOT_COMPATIBLE)
            {
                addNeighborship(source);
            }
        }
    }


    // add nodes with too less neighbors to queue
    for (nG = 0; nG < m_sizeG; ++nG)
    {
        //if (m_nodesG[nG] != m_xG) // omit excluded vertex
        for (nH = 0; nH < m_sizeH; ++nH)
        {
            if (productNodeWeight(source) != WEIGHT_NOT_COMPATIBLE && !m_productNodeNeighbors[nG * m_sizeH + nH].getwasRemoved()
              && !m_productNodeNeighbors[nG * m_sizeH + nH].hasSufficientAmountOfNeighbors())
            {
                    removeNeighborship(source);
            }
        }
    }

}

unsigned MaxComPart2Tree::computeBlock(forward_list<std::vector<Mapping>>& maximumMappingList, wType &maxWeight, const bool &enumerate)
{
    m_enumerate=enumerate;
    maxWeight=WEIGHT_NOT_COMPATIBLE;
    m_maxWeightNonOuterPtr=&maxWeight;

    // clique computation for non outerplanar blocks
    if (!m_outerplanar)
    {
        eliminateNodesOfDegreeOneAndZeroInProductGraph();
        m_cliqueTimeStart=clock();
        //return;
        //cout << "computeBlock non outerplanar\n";
        m_maximumMappingListPtr=&maximumMappingList;

        //if (m_xG != nullptr)
        //	cout << "Exclude vertex "<<m_xG<<endl;
        // C-Clique-Init (G)
        vector<ProductnodeNoExpand> Q,P;//,X,Y;
        stack<ProductnodeNoExpand> sInsufficientNeighbors;

         // Nodes 0..t are Part of set T
        for (int t=0, uin1=0, uin2=0; t<m_nodesProductGraph; ++t)
        {
            ProductnodeNoExpand ui(uin1,uin2);
            // additional code for excluded vertices
            if (m_nodesG[uin1] != m_xG)
            {
                // Zeile 4
                m_R.clear();
                if (productNodeWeight(ui) != WEIGHT_NOT_COMPATIBLE && m_productNodeNeighbors[t].hasSufficientAmountOfNeighbors())
                {
                    m_R.push_back(ui);

                    // lines 5 to 8
                    Q.clear();
                    P.clear();
                    //X.clear();
                    //Y.clear();
                    ProductnodeNoExpand target(uin1,0);
                    // Iterate through all neighbours which were no starting vertices before
                    int &targetnG=target.nodeG, &targetnH=target.nodeH;
                    //for (int targetnG=0, targetnH=0, counter=0;	counter<m_nodesProductGraph; ++counter)
                    //for (int targetnG=(t+1)/sizeH, targetnH=(t+1)%sizeH, counter=t+1; counter<m_nodesProductGraph; ++counter)

                    // remove vertices, which do not have at least two connections through c_edges
                    for (int counter=targetnG*m_sizeH; counter<m_nodesProductGraph; ++counter)
                    {
                        //if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE && !m_productNodeNeighbors[counter].getwasRemoved() && t != counter)
                        if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE && t != counter)
                        {
                            if (getEdgeType(ui,target)==NO_EDGE) // NO_EDGE, therefore eliminate target recursively from neighborship
                            {
                                removeNeighborship(target, &sInsufficientNeighbors);
                            }
                        }
                        ++targetnH;
                        if (targetnH==m_sizeH)
                        {
                            targetnH=0;
                            ++targetnG;
                        }
                    }

                    targetnG=uin1+1; targetnH=0;
                    for (int counter=targetnG*m_sizeH; counter<m_nodesProductGraph; ++counter)
                    {
                        // if (m_nodesG[targetnG] != m_xG && productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE)

                        if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE && !m_productNodeNeighbors[counter].getwasRemoved())
                        {
                            if (getEdgeType(ui,target)==C_EDGE)
                            {
                                //if (counter>t) // P
                                { P.push_back(target);} //cout << "P"<<target; }
                                //else // X
                                //{ X.push_back(target);}// cout << "X"<<target; }
                            }
                            else if (getEdgeType(ui,target)==D_EDGE)
                            {
                                //if (counter>t) // Q
                                { Q.push_back(target);}// cout << "Q"<<target; }
                                //else // Y
                                //{ Y.push_back(target);}// cout << "Y"<<target; }
                            }
                        }
                        // nächstes Element aus D[ui] bzw C[ui] bestimmen
                        ++targetnH;
                        if (targetnH==m_sizeH)
                        {
                            targetnH=0;
                            ++targetnG;
                        }
                    }

                    // Aufruf von C-Clique (Zeile 9)
                    c_Clique(P,Q);
                    if (m_maxCliqueTime>0 && (clock()-m_cliqueTimeStart) / ( CLOCKS_PER_SEC / 1000)> m_maxCliqueTime)
                    {
//cout << m_cliqueRecursions<< " ";
//return m_cliqueRecursions;
                        return 1;  // clique computation time limit exceeded
                    }

                    // restore vertices, which have been excluded through NO_EDGE
                    while (!sInsufficientNeighbors.empty())
                    {
                        addNeighborship(sInsufficientNeighbors.top());
                        sInsufficientNeighbors.pop();
                    }

                    // eliminate last used node recursively from neighborship
                    removeNeighborship(ui);
                }
            }
            else // skip node uin1 if it is excluded
            {
                t+=m_sizeH-1;
                uin2=m_sizeH-1;
            }

            ++uin2;
            if (uin2==m_sizeH)
            {
                uin2=0;
                ++uin1;
            }
        }
        //maxWeight=maxWeightNonOuter;
    //	cout << "Clique Recursions:"<<m_cliqueRecursions<<endl;
//cout << m_cliqueRecursions<< " ";
//return m_cliqueRecursions;
        return 0;
    }

    // quadratic time computation for outerplanar blocks
    List<node> sNodesG = m_SPG->nodesOfType(StaticSPQRTree::NodeType::SNode);
    List<node> sNodesH = m_SPH->nodesOfType(StaticSPQRTree::NodeType::SNode);
    for(ListConstIterator< node > itG = (sNodesG).begin(); itG.valid(); ++itG)	{
//	forall_listiterators(node, itG, sNodesG) {
        node sG = *itG;
        Skeleton &skelSG = m_SPG->skeleton(sG);
        for(ListConstIterator< node > itH = (sNodesH).begin(); itH.valid(); ++itH)	{
//		forall_listiterators(node, itH, sNodesH) {
            node sH = *itH;
            Skeleton &skelSH = m_SPH->skeleton(sH);
            // check if skeleton graphs are compatible
            if (skelSG.getGraph().numberOfNodes() != skelSH.getGraph().numberOfNodes())
                continue;
            edge eGSkel;
            forall_edges(eGSkel, skelSG.getGraph()) {
                edge eG = getRealEdge(eGSkel, skelSG, m_SPG);
                node uGSkel = eGSkel->source();
                node vGSkel = eGSkel->target();
                node uG = skelSG.original(uGSkel);
                node vG = skelSG.original(vGSkel);
                edge eHSkel;
                forall_edges(eHSkel, skelSH.getGraph()) {
                    edge eH = getRealEdge(eHSkel, skelSH, m_SPH);
                    node uHSkel = eHSkel->source();
                    node vHSkel = eHSkel->target();
                    node uH = skelSH.original(uHSkel);
                    node vH = skelSH.original(vHSkel);

                    m_map_snode->addMap(sG,sH);
                    // first orientation
                    addNodeMappingToMCS(uG, uH);
                    addNodeMappingToMCS(vG, vH);

                    if (m_D[eG][eH][detectCase(eG, eH)] == WEIGHT_UNDEF) { // not yet computed
                        matchEdge(uGSkel, vGSkel, uHSkel, vHSkel, sG, sH, eGSkel, eHSkel);
                        matchCycle(uGSkel, vGSkel, uHSkel, vHSkel, sG, sH);
                        applyWeightsToMCS();
                        computeSplitIsomorphisms(maximumMappingList,maxWeight);
                    }
                    deleteCurrentSolution();

                    m_map_snode->addMap(sG,sH);
                    // second orientation: uG -> vH, vG -> uH
                    addNodeMappingToMCS(uG, vH);
                    addNodeMappingToMCS(vG, uH);
                    if (m_D[eG][eH][detectCase(eG, eH)] == WEIGHT_UNDEF) { // not yet computed
                        matchEdge(uGSkel, vGSkel, vHSkel, uHSkel, sG, sH, eGSkel, eHSkel); // todo remove parameter sp tree, use owner (somewhere in code)
                        matchCycle(uGSkel, vGSkel, vHSkel, uHSkel, sG, sH);
                        applyWeightsToMCS();
                        computeSplitIsomorphisms(maximumMappingList,maxWeight);
                    }
                    deleteCurrentSolution();
                }
            }
        }
    }
    return 0;
}

unsigned MaxComPart2Tree::computeFixedMapping(NodeArray<std::forward_list<std::vector<Mapping>>> &aGmaximumMappingList, NodeArray<wType> &maxWeight, const bool &enumerate)
{
    m_enumerate=enumerate;
    node aH,vH;
    forall_nodes(vH,m_H)
        maxWeight[m_vHtoaH[vH]]=WEIGHT_NOT_COMPATIBLE;

    if (!m_outerplanar)
    {
        eliminateNodesOfDegreeOneAndZeroInProductGraph();
        m_cliqueTimeStart=clock();
        //return;
        //cout << "computeFixedMapping with given vertex "<<m_vG<<endl;
        m_maxWeightArrayNonOuterPtr=&maxWeight;
        m_maximumMappingArrayListPtr=&aGmaximumMappingList;
        //m_maximumMappingListPtr=&maximumMappingList;
        // C-Clique-Init (G)
        vector<ProductnodeNoExpand> Q,P;//,X,Y;
        stack<ProductnodeNoExpand> sInsufficientNeighbors;

         // Nodes 0..t are Part of set T
        for (int t=0, uin1=0, uin2=0; t<m_nodesProductGraph; ++t)
        {
            // additional code for fixed vertex
            if (m_nodesG[uin1] == m_vG)
            {
                // Zeile 4
                m_R.clear();
                ProductnodeNoExpand ui(uin1,uin2);
                if (productNodeWeight(ui) != WEIGHT_NOT_COMPATIBLE && m_productNodeNeighbors[t].hasSufficientAmountOfNeighbors())
                {
                    m_R.push_back(ui);

                    // Zeile 5 bis 8
                    Q.clear();
                    P.clear();
                    //X.clear();
                    //Y.clear();
                    ProductnodeNoExpand target(0,0);
                    // Iterate through all neighbours which were no starting vertices before
                    int &targetnG=target.nodeG, &targetnH=target.nodeH;

                    // remove vertices, which do not have at least two connections through c_edges
                    for (int counter=targetnG*m_sizeH; counter<m_nodesProductGraph; ++counter)
                    {
                        if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE && t != counter)
                        {
                            if (getEdgeType(ui,target)==NO_EDGE) // NO_EDGE, therefore eliminate target recursively from neighborship
                            {
                                removeNeighborship(target, &sInsufficientNeighbors);
                            }
                        }
                        ++targetnH;
                        if (targetnH==m_sizeH)
                        {
                            targetnH=0;
                            ++targetnG;
                        }
                    }

                    targetnG=0; targetnH=0;
                    for (int counter=0; counter<m_nodesProductGraph; ++counter)
                    //for (int targetnG=(t+1)/sizeH, targetnH=(t+1)%sizeH, counter=t+1; counter<m_nodesProductGraph; ++counter)
                    //for (int targetnG=uin1+1, targetnH=0, counter=targetnG*sizeH; counter<m_nodesProductGraph; ++counter)
                    {

                        //if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE)
                        //if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE && m_productNodeNeighbors[counter].hasSufficientAmountOfNeighbors())
                        if (productNodeWeight(target) != WEIGHT_NOT_COMPATIBLE && !m_productNodeNeighbors[counter].getwasRemoved())
                        {
                            if (getEdgeType(ui,target)==C_EDGE)
                            {
                                //if (counter>t) // P
                                { P.push_back(target);} //cout << "P"<<target; }
                                //else // X
                                //{ X.push_back(target);}// cout << "X"<<target; }
                            }
                            else if (getEdgeType(ui,target)==D_EDGE)
                            {
                                //if (counter>t) // Q
                                { Q.push_back(target);}// cout << "Q"<<target; }
                                //else // Y
                                //{ Y.push_back(target);}// cout << "Y"<<target; }
                            }
                        }
                        // nächstes Element aus D[ui] bzw C[ui] bestimmen
                        ++targetnH;
                        if (targetnH==m_sizeH)
                        {
                            targetnH=0;
                            ++targetnG;
                        }
                    }
                    // Aufruf von C-Clique (Zeile 9)
                    //c_Clique(P,Q,X,Y);
                    c_Clique(P,Q);
                    if (m_maxCliqueTime>0 && (clock()-m_cliqueTimeStart) /( CLOCKS_PER_SEC /1000)> m_maxCliqueTime)
                    {
//cout << m_cliqueRecursions<< " ";
//return m_cliqueRecursions;
                        return 1;
                    }

                    // restore vertices, which have been excluded through NO_EDGE
                    while (!sInsufficientNeighbors.empty())
                    {
                        addNeighborship(sInsufficientNeighbors.top());
                        sInsufficientNeighbors.pop();
                    }

                    // eliminate last used node recursively from neighborship
                    removeNeighborship(ui);
                }

            }
            else // skip node uin1 if it is not m_vG
            {
                t+=m_sizeH-1;
                uin2=m_sizeH-1;
            }
            // Nächstes ui bestimmen
            ++uin2;
            if (uin2==m_sizeH)
            {
                if (m_nodesG[uin1] == m_vG) // everything on fixed vertex m_vG finished
                    break;
                uin2=0;
                ++uin1;
            }
        }
//cout << m_cliqueRecursions<< " ";
//return m_cliqueRecursions;
        return 0;
    }

    List<node> sNodesG = m_SPG->nodesOfType(StaticSPQRTree::NodeType::SNode);
    List<node> sNodesH = m_SPH->nodesOfType(StaticSPQRTree::NodeType::SNode);
    for(ListConstIterator< node > itG = (sNodesG).begin(); itG.valid(); ++itG)	{
//	forall_listiterators(node, itG, sNodesG) {
        node sG = *itG;
        Skeleton &skelSG = m_SPG->skeleton(sG);
        for(ListConstIterator< node > itH = (sNodesH).begin(); itH.valid(); ++itH)	{
//		forall_listiterators(node, itH, sNodesH) {
            node sH = *itH;
            Skeleton &skelSH = m_SPH->skeleton(sH);
            // check if skeleton graphs are compatible
            if (skelSG.getGraph().numberOfNodes() != skelSH.getGraph().numberOfNodes())
                continue;
            edge eGSkel;
            forall_edges(eGSkel, skelSG.getGraph()) {
                edge eG = getRealEdge(eGSkel, skelSG, m_SPG);
                node uGSkel = eGSkel->source();
                node vGSkel = eGSkel->target();
                node uG = skelSG.original(uGSkel);
                node vG = skelSG.original(vGSkel);
                if (uG != m_vG && vG != m_vG) // at least one vertex of the edge must be of the fixed mapping
                    continue;
                edge eHSkel;
                forall_edges(eHSkel, skelSH.getGraph()) {
                    edge eH = getRealEdge(eHSkel, skelSH, m_SPH);
                    node uHSkel = eHSkel->source();
                    node vHSkel = eHSkel->target();
                    node uH = skelSH.original(uHSkel);
                    vH = skelSH.original(vHSkel);
                    if (vG == m_vG)
                        aH=m_vHtoaH[vH];
                    else
                        aH=m_vHtoaH[uH];
                    m_map_snode->addMap(sG,sH);
                    // first orientation
                    addNodeMappingToMCS(uG, uH);
                    addNodeMappingToMCS(vG, vH);

                    if (m_D[eG][eH][detectCase(eG, eH)] == WEIGHT_UNDEF) { // not yet computed
                        matchEdge(uGSkel, vGSkel, uHSkel, vHSkel, sG, sH, eGSkel, eHSkel);
                        matchCycle(uGSkel, vGSkel, uHSkel, vHSkel, sG, sH);
                        applyWeightsToMCS();
                        //computeSplitIsomorphisms(maximumMapping[aH],maxWeight[aH]);
                        computeSplitIsomorphisms(aGmaximumMappingList[aH],maxWeight[aH]);
                    }
                    deleteCurrentSolution();

                    if (vG == m_vG)
                        aH=m_vHtoaH[uH];
                    else
                        aH=m_vHtoaH[vH];
                    m_map_snode->addMap(sG,sH);
                    // second orientation: uG -> vH, vG -> uH
                    addNodeMappingToMCS(uG, vH);
                    addNodeMappingToMCS(vG, uH);
                    if (m_D[eG][eH][detectCase(eG, eH)] == WEIGHT_UNDEF) { // not yet computed
                        matchEdge(uGSkel, vGSkel, vHSkel, uHSkel, sG, sH, eGSkel, eHSkel); // todo remove parameter sp tree, use owner (somewhere in code)
                        matchCycle(uGSkel, vGSkel, vHSkel, uHSkel, sG, sH);
                        applyWeightsToMCS();
                        //computeSplitIsomorphisms(maximumMapping[aH],maxWeight[aH]);
                        computeSplitIsomorphisms(aGmaximumMappingList[aH],maxWeight[aH]);
                    }
                    deleteCurrentSolution();
                }
            }
        }
    }
    return 0;
}


// note: parameters are skeleton elements
void MaxComPart2Tree::matchCycle(const node uGSkel, const node vGSkel, const node uHSkel, const node vHSkel, const node sG, const node sH) {
  
    node endG = uGSkel;
    node lastG = uGSkel;
    node currentG = vGSkel;
    node lastH = uHSkel;
    node currentH = vHSkel;

    // initial vertices
    edge eG=nullptr,eH=nullptr;
    for(adjEntry adj : currentG->adjEntries)
    {
        eG=adj->theEdge();
    //forall_adj_edges(eG,currentG) {
        if (eG->opposite(currentG) != lastG) {
            lastG = currentG;
            currentG = eG->opposite(currentG);
            break;
        }
    }
    for(adjEntry adj : currentH->adjEntries)
    {
        eH=adj->theEdge();
    //forall_adj_edges(eH,currentH) {
        if (eH->opposite(currentH) != lastH) {
            lastH = currentH;
            currentH = eH->opposite(currentH);
            break;
        }
    }

    while (currentG != endG) {
        addNodeMappingToMCS(m_SPG->skeleton(sG).original(currentG), m_SPH->skeleton(sH).original(currentH));
        matchEdge(lastG, currentG, lastH, currentH, sG, sH, eG, eH);

        // update vertices
        for(adjEntry adj : currentG->adjEntries)
        {
            eG=adj->theEdge();
        //forall_adj_edges(eG,currentG) {
            if (eG->opposite(currentG) != lastG) {
                lastG = currentG;
                currentG = eG->opposite(currentG);
                break;
            }
        }
        for(adjEntry adj : currentH->adjEntries)
        {
            eH=adj->theEdge();
        //forall_adj_edges(eH,currentH) {
            if (eH->opposite(currentH) != lastH) {
                lastH = currentH;
                currentH = eH->opposite(currentH);
                break;
            }
        }
    }
    matchEdge(lastG, currentG, lastH, currentH, sG, sH, eG, eH);
}

edge MaxComPart2Tree::getRealEdge(const edge sNodeSkeletonEdge, const Skeleton &skel, const StaticSPQRTree *sprqtree) {
    if (skel.isVirtual(sNodeSkeletonEdge)) {
        node pG = skel.twinTreeNode(sNodeSkeletonEdge);
        edge eG;
        forall_edges(eG, sprqtree->skeleton(pG).getGraph()) {
            if (!sprqtree->skeleton(pG).isVirtual(eG)) break;
        }
        return sprqtree->skeleton(pG).realEdge(eG);
    } else {
        return skel.realEdge(sNodeSkeletonEdge);
    }
}

// note: parameters are skeleton elements
void MaxComPart2Tree::matchEdge(const node uGSkel, const node vGSkel, const node uHSkel, const node vHSkel, const node sG, const node sH, const edge eGSkel, const edge eHSkel) {
//	cout << "entered matching" << endl;
    bool eGVirtual = m_SPG->skeleton(sG).isVirtual(eGSkel);
    bool eHVirtual = m_SPH->skeleton(sH).isVirtual(eHSkel);
    edge eG = getRealEdge(eGSkel, m_SPG->skeleton(sG), m_SPG);
    edge eH = getRealEdge(eHSkel, m_SPH->skeleton(sH), m_SPH);

//	cout << "the edge: " << eGSkel << " " << eHSkel << endl;
//	cout << "the real edge: " << eG << " " << eH << endl;
    addEdgeMappingToMCS(eG, eH);

//	cout << "mapped " << eG << " ------> " << eH << endl;

//	cout << "determined virtual or real" << endl;
    if (eGVirtual && eHVirtual) {
//		cout << "complex matching required" << endl;
        node pG = m_SPG->skeleton(sG).twinTreeNode(eGSkel);
        node pH = m_SPH->skeleton(sH).twinTreeNode(eHSkel);

        node sGPrime=nullptr;
        edge e;
        for(adjEntry adj : pG->adjEntries)
        {
            e=adj->theEdge();
        //forall_adj_edges(e,pG) {
            node u = e->opposite(pG);
            if (u != sG) sGPrime = u;
        }

        node sHPrime=nullptr;
        for(adjEntry adj : pH->adjEntries)
        {
            e=adj->theEdge();
        //forall_adj_edges(e,pH) {
            node u = e->opposite(pH);
            if (u != sH) sHPrime = u;
        }

        // check compatibility
        if (m_SPG->skeleton(sGPrime).getGraph().numberOfNodes() ==
                m_SPH->skeleton(sHPrime).getGraph().numberOfNodes() ) {


            node uGPrime = getCopy(m_SPG->skeleton(sG).original(uGSkel),m_SPG->skeleton(sGPrime));
            node vGPrime = getCopy(m_SPG->skeleton(sG).original(vGSkel),m_SPG->skeleton(sGPrime));
            node uHPrime = getCopy(m_SPH->skeleton(sH).original(uHSkel),m_SPH->skeleton(sHPrime));
            node vHPrime = getCopy(m_SPH->skeleton(sH).original(vHSkel),m_SPH->skeleton(sHPrime));

            m_map_snode->addMap(sGPrime, sHPrime);
            matchCycle(uGPrime, vGPrime, uHPrime, vHPrime, sGPrime, sHPrime);
        }

    }
}

node MaxComPart2Tree::getCopy(const node v, const Skeleton &S)
{
    //node skelV;
    //forall_nodes(skelV, S.getGraph())
    for(node skelV: S.getGraph().nodes)
    {
        if(S.original(skelV) == v)
            return skelV;
    }
    return NULL;
}
