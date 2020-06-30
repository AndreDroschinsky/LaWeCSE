#include <bbpmcsi.h>
//#include <test.h>
#include <queue>
#include <list>

//#include <ogdf/planarity/PlanarizationLayout.h>

// Construct BC-Trees, select initial B-Node, computes embedding into the plane
BBP_MCSI::BBP_MCSI(LabelFunction &labelFunction, InputGraph &iG, InputGraph &iH, bool enumerate, bool graphical_display, const vector<string>* simpleLabelToString, clock_t maxCliqueTime, wType distancePenalty)
:m_labelFunction(labelFunction),m_ig_G(&iG),m_ig_H(&iH),m_G(&iG.graph),m_H(&iH.graph),
     m_SimpleLabelToString(simpleLabelToString), m_enumerate(enumerate),
     m_VbcG(iG.m_VbcG), m_VbcH(iH.m_VbcG), m_bGAG(iG.m_bGAG), m_bHAH(iH.m_bGAG), m_aGtoSPaG(iG.m_aGtoSPaG), m_aHtoSPaH(iH.m_aGtoSPaG),
     m_SPaGtoaG(iG.m_SPaGtoaG), m_SPaHtoaH(iH.m_SPaGtoaG), m_SPeGtoaeG(iG.m_SPeGtoaeG), m_SPeHtoaeH(iH.m_SPeGtoaeG), m_bGOuterplanar(iG.m_bGOuterplanar), m_bHOuterplanar(iH.m_bGOuterplanar),
     m_bGToSPQRTree(iG.m_bGToSPQRTree), m_bHToSPQRTree(iH.m_bGToSPQRTree),
     m_weightBBPIso(WEIGHT_UNDEF),m_maxCliqueTime(maxCliqueTime),m_distancePenalty(distancePenalty)
#ifdef GRAPHICS
, m_graphical_display(graphical_display)
#endif
{
    Graph &G = iG.graph;
    Graph &H = iH.graph;
    m_nodeLabelSimpleG=iG.nodeLabel;
    m_nodeLabelSimpleH=iH.nodeLabel;
    m_edgeLabelSimpleG=iG.edgeLabel;
    m_edgeLabelSimpleH=iH.edgeLabel;
    m_simpleLabelG=iG.simpleLabel;
    m_simpleLabelH=iH.simpleLabel;

    // Compute BC Trees and vertices V(b) for all B-Nodes of G and H and the single parts of their auxiliary graphs if not done yet
    m_BC_G=initVbcAndSingleParts(iG);
    m_BC_H=initVbcAndSingleParts(iH);

/*	if (iG.bcTree == nullptr)
    {
        iG.bcTree=m_BC_G=new BCTree(iG.graph);
        initVbcAndSingleParts(iG);
        //initVbcAndSingleParts(*m_BC_G,m_VbcG,m_bGAG,m_aGtoSPaG,m_SPaGtoaG,m_SPeGtoaeG,m_bGOuterplanar);
    }
    else
    {
        m_BC_G=iG.bcTree;
    }

    if (iH.bcTree == nullptr)
    {
        iH.bcTree=m_BC_H=new BCTree(iH.graph);
        initVbcAndSingleParts(iH);
        //initVbcAndSingleParts(*m_BC_H,m_VbcH,m_bHAH,m_aHtoSPaH,m_SPaHtoaH,m_SPeHtoaeH,m_bHOuterplanar);
    }
    else
    {
        m_BC_H=iH.bcTree;
    }*/

    m_BCG=&m_BC_G->bcTree();
    m_BCH=&m_BC_H->bcTree();
    m_AG=&m_BC_G->auxiliaryGraph();
    m_AH=&m_BC_H->auxiliaryGraph();

    // Select initial B-Node
    m_b_initial=m_BC_G->bcTree().firstNode();


    // Initialize local mappings between two blocks
    m_mcp2txG.init(*m_BCG);
    m_mcp2tvG.init(*m_BCG);

    node b;
    forall_nodes(b,*m_BCG)
        if (m_VbcG[b].size()>2)
        {
            m_mcp2txG[b].init(*m_BCH,nullptr);
            m_mcp2tvG[b].init(*m_BCH,nullptr);
        }

    // Initialize weights and local mappings for BBP_Edge, where vG, vH are nullptrs
    m_Weight_BBPE_bGTobH.init(*m_BCG);
    m_LocalIso_BBPE_bGTobH.init(*m_BCG);
    m_LocalIso_BBPE_bGTobH_first.init(*m_BCG);
    if (m_enumerate)
    {
        m_EnumBlockBridgeMappingIterator_bGTobH.init(*m_BCG);
    }
    node n;
    forall_nodes(n,*m_BCG)
    {
        m_Weight_BBPE_bGTobH[n].init(*m_BCH,WEIGHT_UNDEF);
        m_LocalIso_BBPE_bGTobH[n].init(*m_BCH,forward_list<vector<Mapping>>());
        m_LocalIso_BBPE_bGTobH_first[n].init(*m_BCH, false);
        if (m_enumerate)
        {
            m_EnumBlockBridgeMappingIterator_bGTobH[n].init(*m_BCH);
            /*node nH;
            forall_nodes(nH,*m_BCH)
            {
                m_EnumBlockBridgeMappingIterator_bGTobH[n][nH]=m_LocalIso_BBPE_bGTobH[n][nH].end();
            }*/
        }
    }
    // Initialize weights and local mappings for BBP_Edge, where vG, vH are set
    m_Weight_BBPE_aGToaH.init(*m_AG);
    m_LocalIso_BBPE_aGToaH.init(*m_AG);
    if (m_enumerate)
    {
        m_EnumBlockBridgeMappingIterator_aGToaH.init(*m_AG);
    }
    forall_nodes(n,*m_AG)
    {
        m_Weight_BBPE_aGToaH[n].init(*m_AH,WEIGHT_UNDEF);
        //m_LocalIso_BBPE_aGToaH[n].init(*m_AH,vector<Mapping>());
        m_LocalIso_BBPE_aGToaH[n].init(*m_AH,forward_list<vector<Mapping>>());
        if (m_enumerate)
        {
            m_EnumBlockBridgeMappingIterator_aGToaH[n].init(*m_AH);
            node nH;
            forall_nodes(nH,*m_AH)
            {
                m_EnumBlockBridgeMappingIterator_aGToaH[n][nH]=m_LocalIso_BBPE_aGToaH[n][nH].end();
            }
        }
    }

    // Initialize weights and local mapping pointers for BBP_SingleVertex
    m_Weight_BBPSV_vGTovH.init(G);
    m_LocalIso_BBPSV_vGTovH.init(G);
    /*if (m_enumerate)
    {
        m_EnumBlockBridgeSV_Iterator.init(G);
    }*/
    forall_nodes(n,G)
    {
        m_Weight_BBPSV_vGTovH[n].init(H,WEIGHT_UNDEF);
        m_LocalIso_BBPSV_vGTovH[n].init(H,nullptr);
        /*m_LocalIso_BBPSV_vGTovH[n].init(H,vector< forward_list<vector<Mapping>>* >());
        if (m_enumerate)
        {
            m_EnumBlockBridgeSV_Iterator[n].init(H);
            node nH;
            forall_nodes(nH,H)
            {
                m_EnumBlockBridgeSV_Iterator[n][nH]=m_LocalIso_BBPSV_vGTovH[n][nH].end();
            }
        }*/
    }

    // Initialize MWM for BBP_Edge and subinstances
    m_Matching.init(G);
    m_LaWeCSExpand.init(G);
    if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE)
    {
        m_SkippedVertex.init(G);
        m_LaWeCSMapping.init(G);
        m_LaWeCS_SkippedRootWeight.init(G);
        m_LaWeCS_SkippedRootMapping.init(G);
    }
    m_submatching_index.init(G);
    m_submatching_block.init(G);
    m_instance_ID.init(*m_AH,0);
    forall_nodes(n,G)
    {
        if (m_BC_G->typeOfGNode(n) == BCTree::GNodeType::CutVertex)
        {
            m_Matching[n].init(H, nullptr);

            m_submatching_index[n].init(H);
            m_submatching_block[n].init(H, nullptr);
        }
        m_LaWeCSExpand[n].init(*m_AH, MWM_EXPANSION);
        if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE)
        {
            m_SkippedVertex[n].init(H, nullptr);
            m_LaWeCSMapping[n].init(*m_AH,Mapping());
        }
    }

    // Store CVs of H
    m_CVH.reserve(m_BC_G->numberOfCComps());
    forall_nodes(n,H)
        if (m_BC_H->typeOfGNode(n) == BCTree::GNodeType::CutVertex)
            m_CVH.push_back(n);

    adjEntry adjBH;
    for (node cv_orig:m_CVH)
    {
        unsigned bHIndex = 0;
        node cv_BC=m_BC_H->bcproper(cv_orig);
        forall_adj(adjBH,cv_BC)
        {
            m_instance_ID[m_BC_H->repVertex(cv_orig,adjBH->twinNode())]=bHIndex++;
        }
        m_instance_ID[m_BC_H->rep(cv_orig)]=bHIndex;
    }

#ifndef LAWECSU_NDEBUG
    test::PrintBCTree(*m_BC_G);
    forall_nodes(n,*m_BCG)
    {
        cout << "Vertices in component "<<n<<". ";
        for (node v : m_VbcG[n])
            cout << m_BC_G->original(v)<<",";
        cout <<endl;
    }
    cout <<endl;
    test::PrintBCTree(*m_BC_H);
    forall_nodes(n,*m_BCH)
    {
        cout << "Vertices in component "<<n<<". ";
        for (node v : m_VbcH[n])
            cout << m_BC_H->original(v)<<",";
        cout <<endl;
    }
    cout <<endl;
#endif

#ifdef GRAPHICS
    // Initialize graphical display
    if (m_graphical_display)
    {
        m_graphVisualizationG=new GraphVisualization;
        m_graphVisualizationH=new GraphVisualization;
        m_graphVisualizationG->init(G,*m_BC_G,sf::Color(245,252,245,255),m_nodeLabelSimpleG,m_edgeLabelSimpleG,m_SimpleLabelToString,m_simpleLabelG,iG.fogNodeEdgeOffset);
        m_graphVisualizationH->init(H,*m_BC_H,sf::Color(245,250,252,255),m_nodeLabelSimpleH,m_edgeLabelSimpleH,m_SimpleLabelToString,m_simpleLabelH,iH.fogNodeEdgeOffset);
    }
#endif
}

BBP_MCSI::~BBP_MCSI()
{
    // Optimize this, to consider only relevant CVs as aG
    node cG;
    forall_nodes(cG, *m_G)
    {
        if (m_BC_G->typeOfGNode(cG)== BCTree::GNodeType::CutVertex)
            for (node cH : m_CVH)
            {
                if (m_Matching[cG][cH]!= nullptr)
                    delete m_Matching[cG][cH];
            }
        node vH;
        if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE)
        {
            forall_nodes(vH, *m_H)
            {
                if (m_SkippedVertex[cG][vH]!= nullptr)
                    delete m_SkippedVertex[cG][vH];
            }
        }
    }
    node bG,bH;
    forall_nodes(bG,*m_BCG)
        if (m_VbcG[bG].size()>2)
            forall_nodes(bH,*m_BCH)
            {
                delete m_mcp2txG[bG][bH];
                if (bG!=m_b_initial)
                    delete m_mcp2tvG[bG][bH];
            }
    //delete(m_BC_G);
    //delete(m_BC_H);
#ifdef GRAPHICS
    if (m_graphical_display)
    {
        delete m_graphVisualizationG;
        delete m_graphVisualizationH;
    }
#endif
}

void BBP_MCSI::compute_BC_G_parents(node v)
{
    adjEntry adjv;
    forall_adj(adjv, v)
    {
        node neighbour=adjv->twinNode();
        if (m_parent_BC_G[neighbour] == nullptr)
        {
            m_parent_BC_G[neighbour] = v;
            compute_BC_G_parents(neighbour);
        }
    }
}

// computes the weight of a BBPMCSI between G and H
wType BBP_MCSI::computeSize()
{
    if (m_weightBBPIso == WEIGHT_UNDEF)
    {
        m_parent_BC_G.init(*m_BCG,nullptr);
        m_parent_BC_G[m_b_initial] = m_b_initial;
        compute_BC_G_parents(m_b_initial);

        //computeMatchings(m_b_initial,nullptr);
        m_weightBBPIso=setSx();
    }
    return m_weightBBPIso;
}

SkippedVertex::SkippedVertex(unsigned neighbours)
{
    m_weight.resize(neighbours, WEIGHT_UNDEF);
}

wType SkippedVertex::getWeight(unsigned jobIndex) const
{
    //if (jobIndex >= 0 && jobIndex<m_weight.size())
        return m_weight[jobIndex];
    /*else
    {
        cerr << "SkippedVertex::getWeight out of range\n";
        exit(1);
    }*/
}

void SkippedVertex::setWeight(unsigned jobIndex, wType weight)
{
    //if (jobIndex >= 0 && jobIndex<m_weight.size())
    //{
        m_weight[jobIndex]=weight;
        //cout << "weight of "<<jobIndex<< " set to "<< weight<<endl;
    /*}
    else
    {
        cerr << "SkippedVertex::setWeight out of range\n";
        exit(1);
    }*/
}

void SkippedVertex::updateStoredBestSolutions(node aH, wType weight, const Mapping& mapping)
{
    if (aH == m_best_aH) // special case if current best solution is from same aH
    {
        if (weight > m_best_weight) // new best solution; a 2nd best solution may be ignored
        {
            m_best_weight = weight;
            m_best_mapping = mapping;
        }
    }
    else
    {
        if (weight > m_best_weight) // new best solution
        {
            m_2nd_best_weight = m_best_weight;
            m_2nd_best_aH = m_best_aH;
            m_2nd_best_mapping = m_best_mapping;

            m_best_weight = weight;
            m_best_aH = aH;
            m_best_mapping = mapping;
        }
        else if (weight > m_2nd_best_weight) // equal or worse than best, better than second best
        {
            m_2nd_best_weight = weight;
            m_2nd_best_aH = aH;
            m_2nd_best_mapping = mapping;
        }
    }
}

wType SkippedVertex::getBestWeightAndMapping(node excluded, Mapping &mapping) const
{
    if (excluded == m_best_aH)
    {
        mapping = m_2nd_best_mapping;
        return m_2nd_best_weight;
    }
    mapping = m_best_mapping;
    return m_best_weight;
}

wType BBP_MCSI::getSkippedVertexValue(const node vG, const node aH)
{
    // MCS, no skipped vertices allowed
    if (m_distancePenalty == WEIGHT_NOT_COMPATIBLE)
        return WEIGHT_NOT_COMPATIBLE;

    node vH=m_BC_H->original(aH);

    SkippedVertex* skippedVertex = m_SkippedVertex[vG][vH];
    if (skippedVertex == nullptr)
    {
        //cout << "children nodes are "<<vG<<"  "<<aH;//<<endl;
        //if (m_BC_H->typeOfGNode(vH)==m_BC_G->GNodeType::CutVertex)
        //cout <<" m_instance_ID:"<<m_instance_ID[aH]<<endl;
        int degree = m_BC_H->bcproper(vH)->degree();
        if (degree == 0) // special case of one isolated block only. Then we cannot skip vertices
        {
            degree++;
        }
        m_SkippedVertex[vG][vH] = skippedVertex = new SkippedVertex(degree);
    }
    if (skippedVertex ->getWeight(m_instance_ID[aH]) == WEIGHT_UNDEF)
    {
        //cout << "computing...vG="<<vG<<" d(vG)="<<vG->degree()<< " ";
        //cout << "rep vg="<<m_BC_G->rep(vG)<<" ";
        wType weight = WEIGHT_NOT_COMPATIBLE;
        adjEntry adj;

        // maximum solution from children of aH
        if (m_BC_H->typeOfGNode(vH)== BCTree::GNodeType::CutVertex) // only CVs have children
        {
            if (skippedVertex -> getComputedInstance() == -1) // first sub instance
            {
                node cH=m_BC_H->bcproper(vH);
                forall_adj(adj,cH)
                {
                    node childBlock=adj->twinNode();
                    if (m_VbcH[childBlock].size() == 2) // only bridges are allowed
                    {
                        node childaH = (m_BC_H->original(m_VbcH[childBlock][0]) == vH ? m_VbcH[childBlock][1]:m_VbcH[childBlock][0]);
                        if (m_BC_H->repVertex(vH,childBlock) != aH) // those bridges need to be children in H
                        {
                            // solution from skipped vertex
                            wType wSkipped = getSkippedVertexValue(vG,childaH) - m_distancePenalty;

                            // update best and 2nd best solution
                            if (wSkipped != WEIGHT_NOT_COMPATIBLE)
                                skippedVertex -> updateStoredBestSolutions(childaH, wSkipped, m_LaWeCSMapping[vG][childaH]);

                            if (wSkipped >= weight)
                            {
                                weight = wSkipped;
                                m_LaWeCSMapping[vG][aH] = m_LaWeCSMapping[vG][childaH];
                            }

                            // solution from matching
                            bool expand;
                            wType wMatching = getWeightPlusCVMatching(vG,childaH,expand) - m_distancePenalty;

                            // update best and 2nd best solution
                            if (wMatching != WEIGHT_NOT_COMPATIBLE)
                                skippedVertex -> updateStoredBestSolutions(childaH, wMatching, Mapping(vG,childaH,expand));

                            if (wMatching >= weight)
                            {
                                weight = wMatching;
                                m_LaWeCSMapping[vG][aH] = Mapping(vG,childaH,expand);
                            }
                        }
                        else // store child for later computation
                        {
                            skippedVertex -> setSkippedChild(childaH);
                            skippedVertex -> setComputedInstance(m_instance_ID[aH]);
                        }
                    }
                }
            }
            else // other sub instances are derived from main instance
            {
                if (skippedVertex -> getComputedInstance() > -1) // 2nd sub instance. Then check if previously skipped child yields a new best or 2nd best solution.
                {
                    skippedVertex -> setComputedInstance(-2);
                    node childSkippedaH = skippedVertex -> getSkippedChild();

                    // solution from skipped vertex
                    wType wSkipped = getSkippedVertexValue(vG, childSkippedaH) - m_distancePenalty;

                    // update best and 2nd best solution
                    if (wSkipped != WEIGHT_NOT_COMPATIBLE)
                        skippedVertex -> updateStoredBestSolutions(childSkippedaH, wSkipped, m_LaWeCSMapping[vG][childSkippedaH]);

                    // solution from matching
                    bool expand;
                    wType wMatching = getWeightPlusCVMatching(vG,childSkippedaH,expand) - m_distancePenalty;

                    // update best and 2nd best solution
                    if (wMatching != WEIGHT_NOT_COMPATIBLE)
                        skippedVertex -> updateStoredBestSolutions(childSkippedaH, wMatching, Mapping(vG,childSkippedaH,expand));
                }

                // all other solutions are best, if best is not excluded; otherwise 2nd best
                node childExcluded = nullptr;
                if (aH->degree() == 1) // if aH is in a bridge, we need to exclude that bridge
                {
                    childExcluded = aH->firstAdj()->twinNode();
                }
                Mapping mBest;
                wType wBest = skippedVertex -> getBestWeightAndMapping(childExcluded, mBest);
                if (wBest > weight)
                {
                    weight = wBest;
                    m_LaWeCSMapping[vG][aH] = mBest;
                }
            }
        }

        // maximum solution from children of vG
        if (m_BC_G->typeOfGNode(vG)== BCTree::GNodeType::CutVertex) // only CVs have children
        {
            node cG=m_BC_G->bcproper(vG);
            forall_adj(adj,cG)
            {
                node childBlock=adj->twinNode();
                if (m_parent_BC_G[childBlock] == cG && m_VbcG[childBlock].size() == 2) // only children in G that are bridges are allowed
                {
                    node childaG = (m_BC_G->original(m_VbcG[childBlock][0]) == vG ? m_VbcG[childBlock][1]:m_VbcG[childBlock][0]);
                    //node childvG = (m_BC_G->original(m_VbcG[childBlock][0]) == vG ? m_BC_G->original(m_VbcG[childBlock][1]):m_BC_G->original(m_VbcG[childBlock][0]));
                    node childvG = m_BC_G->original(childaG);
                    wType wSkipped = getSkippedVertexValue(childvG,aH) - m_distancePenalty;
                    if (wSkipped > weight)
                    {
                        weight = wSkipped;
                        m_LaWeCSMapping[vG][aH] = m_LaWeCSMapping[childvG][aH];
                    }
                    bool expand;
                    wType wMatching = getWeightPlusCVMatching(childvG,aH,expand) - m_distancePenalty;
                    if (wMatching >= weight)
                    {
                        weight = wMatching;
                        m_LaWeCSMapping[vG][aH] = Mapping(childvG,aH,expand);
                        //cout << "Added LawMap "<<childvG<<","<<m_BC_H->original(aH)<<" in "<<vG<<","<<m_BC_H->original(aH)<<"expand="<<expand<<endl;
                    }
                }
            }
        }
        skippedVertex->setWeight(m_instance_ID[aH], weight);
    }
    return skippedVertex->getWeight(m_instance_ID[aH]);
}

wType BBP_MCSI::computeMatchingEdgeWeight(node bG, node bH, node vG, node vH, wType &wMatching, wType &wSkippedVertex)
{
    wMatching = WEIGHT_NOT_COMPATIBLE;
    if (compatible(vG,vH))
        wMatching = bbpEdge(bG, bH, nullptr,vG,vH)-w(vG,vH);

    // weight with skipped vertex (topological path)
    wSkippedVertex = WEIGHT_NOT_COMPATIBLE;
    if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE)
    {
        if (m_VbcG[bG].size() == 2 && m_VbcH[bH].size() == 2) // skipped vertices only for bridges
        {

            node childvG = (m_BC_G->original(m_VbcG[bG][0]) == vG ? m_BC_G->original(m_VbcG[bG][1]):m_BC_G->original(m_VbcG[bG][0])); //
            node childaH = (m_BC_H->original(m_VbcH[bH][0]) == vH ? m_VbcH[bH][1]:m_VbcH[bH][0]);

            wSkippedVertex = getSkippedVertexValue(childvG, childaH);

            if (wSkippedVertex != WEIGHT_NOT_COMPATIBLE || wMatching != wMatching)
            {
                if (wSkippedVertex > wMatching)
                {
                    m_LaWeCSExpand[childvG][childaH] = SKIPPED_VERTEX;
                    //cout << "skipped:"<<childvG<<","<<childaH<<" "<<endl;
                }
                else if (wMatching > wSkippedVertex)
                {
                    m_LaWeCSExpand[childvG][childaH] = MWM_EXPANSION;
                }
                else
                {
                    m_LaWeCSExpand[childvG][childaH] = MWM_AND_SKIPPED_VERTEX;
                }
            }
        }
    }
    return max(wMatching, wSkippedVertex);
}

wType BBP_MCSI::getMatchingValue(const node vG, const node aH)
{
    node vH=m_BC_H->original(aH);
    MWM* matching = m_Matching[vG][vH];
    // compute matching values if not done yet
    if (matching == nullptr)
    {
        //cout << "matching not computed\n";
        adjEntry adjBG,adjBH;
        node bG_,cH,bH_;//cG,cvG,bG_,;
        unsigned bGIndex,bHIndex;
        wType weightMax, weightMatching, weightSkippedVertex;

        // Reference on new instance
        m_Matching[vG][vH] = matching = new MWM(m_enumerate);
        //m_SkippedVertex[cvG][cvH] = new SkippedVertex(m_enumerate);
        //caG=BC_G->repVertex(cvG,bG);
        //m_Matching[cvG][cvH]=matching;

        // Add Agents
        node cG=m_BC_G->bcproper(vG);
        node bG=m_parent_BC_G[cG];
        forall_adj(adjBG,cG)
        {
            bG_=adjBG->twinNode();
            if (bG_ != bG)
                matching->addAgent(bG_);
        }
        // Add Jobs
        cH=m_BC_H->bcproper(vH);
        forall_adj(adjBH,cH)
            matching->addJob(adjBH->twinNode());

        // Duplicate vertices
        //matching->addGraphCopy();

        // Add weights
        bGIndex=0;
        forall_adj(adjBG,cG)
        {
            bG_=adjBG->twinNode();
            if (bG_ != bG)
            {
                bHIndex=0;
                forall_adj(adjBH,cH)
                {

                    bH_=adjBH->twinNode();
                    //weight=bbpEdge(bG_,bH_,xG,cvG,cvH)-w(cvG,cvH);
                    if (m_BC_H->repVertex(vH,bH_) != aH) // omit skipped component
                    {
                        weightMax = computeMatchingEdgeWeight(bG_, bH_, vG, vH, weightMatching, weightSkippedVertex);
                        if (weightMax>0 || (weightMax == 0 && m_enumerate))
                            matching->addEdge(bGIndex,bHIndex,weightMax);
                    }
                    else
                    {
                        m_submatching_index[vG][m_BC_H->original(aH)]=bHIndex;
                        m_submatching_block[vG][m_BC_H->original(aH)]=bH_;
                    }
                    ++bHIndex;
                }
                ++bGIndex;
            }
        }
        //matching->addGraphCopy();
#ifndef LAWECSU_NDEBUG
        cout << "MatchingCVs ("<<cvG<<","<< cvH<<"). ";
#endif
        //matching->compute();

        // Compute references to matching bipartite graph instances
        /*bHIndex=0;
        forall_adj(adjBH,cH)
            m_MWM_ID[vG][m_BC_H->repVertex(vH,adjBH->twinNode())]=bHIndex++;
        m_MWM_ID[vG][m_BC_H->rep(vH)]=bHIndex;*/
#ifndef LAWECSU_NDEBUG
        node aH=m_BC_H->rep(cvH);
        getMatchingValue(cvG,aH);
        forall_adj(adjBH,cH)
            getMatchingValue(cvG,m_BC_H->repVertex(cvH,adjBH->twinNode()));
#endif
    }
    else if (!matching->computedAllMatchings()) // only one submatching computed yet
    {
        node bH_=m_submatching_block[vG][vH];
        if ((int) m_instance_ID[aH] != matching->computedMatchings() && bH_ != nullptr) // we need another matching now
        {
            // add missing edges of full matching instance
            node cG=m_BC_G->bcproper(vG);
            unsigned bGIndex=0;
            unsigned bHIndex=m_submatching_index[vG][vH];
            adjEntry adjBG;
            node bG=m_parent_BC_G[cG];
            MWM* matching=m_Matching[vG][vH];
            wType weightMax, weightMatching, weightSkippedVertex;
            forall_adj(adjBG,cG)
            {
                node bG_=adjBG->twinNode();
                if (bG_ != bG)
                {
                    weightMax = computeMatchingEdgeWeight(bG_, bH_, vG, vH, weightMatching, weightSkippedVertex);
                    if (weightMax>0 || (weightMax == 0 && m_enumerate))
                        matching->addEdge(bGIndex,bHIndex,weightMax);
                    ++bGIndex;
                }
            }

        }
    }
    return matching->getWeight(m_instance_ID[aH]);
}

const vector<MWM::AgentJobPair>& BBP_MCSI::getMatching(const node vG, const node aH) const
{
    /*if (m_Matching[vG][m_BC_H->original(aH)]==nullptr)
    {
        cerr << "Matching not computed\n";
        exit (1);
    }*/
    //cout << "mwmid="<<m_MWM_ID[vG][aH]<<endl;
    //return m_Matching[vG][m_BC_H->original(aH)]->getMatching(m_MWM_ID[vG][aH]);
    return m_Matching[vG][m_BC_H->original(aH)]->getMatching(m_instance_ID[aH]);

}

wType BBP_MCSI::w(const node vG, const node vH) const
{
    //return 0;
    if (m_nodeLabelSimpleG != nullptr)
    {
        auto it = m_labelFunction.weightTable.find(LABEL_MULTIPLYER * (*m_nodeLabelSimpleG)[vG] + (*m_nodeLabelSimpleH)[vH]);
        if (it == m_labelFunction.weightTable.end()) // individual entry not found
        {
            if ((*m_nodeLabelSimpleG)[vG]==(*m_nodeLabelSimpleH)[vH])
            {
//				if (1.0 == m_labelFunction.sameNodeLabel)
                    return m_labelFunction.sameNodeLabel;
/*				else
                {
                    cerr << "labelFunction sameNodeLabel not working\n";
                    exit(1);
                }*/
            }
            else
            {
                // if (WEIGHT_NOT_COMPATIBLE == m_labelFunction.differentNodeLabel)
                    return m_labelFunction.differentNodeLabel;
                /*else
                {
                    cerr << "labelFunction differentNodeLabel not working\n";
                    exit(1);
                }*/
            }
        }
        else
        {
            return it->second; // individual entry
        }
    }
    else
        return 1.0;
}

wType BBP_MCSI::w(const edge eG, const edge eH) const
{
    //return 0.0;
    if (m_edgeLabelSimpleG != nullptr)
    {
        auto it = m_labelFunction.weightTable.find(LABEL_MULTIPLYER * (*m_edgeLabelSimpleG)[eG] + (*m_edgeLabelSimpleH)[eH]);
        if (it == m_labelFunction.weightTable.end()) // individual entry not found
        {
            if ((*m_edgeLabelSimpleG)[eG]==LABEL_NOT_COMPATIBLE || (*m_edgeLabelSimpleH)[eH]==LABEL_NOT_COMPATIBLE)
                return WEIGHT_NOT_COMPATIBLE;
            else if ((*m_edgeLabelSimpleG)[eG]==(*m_edgeLabelSimpleH)[eH])
            {
                //if (1.0 == m_labelFunction.sameEdgeLabel)
                    return m_labelFunction.sameEdgeLabel;
                /*else
                {
                    cerr << "labelFunction sameEdgeLabel not working\n";
                    exit(1);
                }*/
            }
            else
            {
                //if (WEIGHT_NOT_COMPATIBLE == m_labelFunction.differentEdgeLabel)
                    return m_labelFunction.differentEdgeLabel;
                /*else
                {
                    cerr << "labelFunction differentEdgeLabel not working\n";
                    exit(1);
                }*/
            }
        }
        else
        {
            return it->second; // individual entry
        }
    }
    else
        return 1.0;
}

wType BBP_MCSI::getWeightPlusCVMatching(const node vG, const node aH, bool &expand)
{
    node vH=m_BC_H->original(aH);
    wType weight=w(vG,vH);
    expand = false;
    if (weight != WEIGHT_NOT_COMPATIBLE && m_BC_G->typeOfGNode(vG)== BCTree::GNodeType::CutVertex && m_BC_H->typeOfGNode(vH)== BCTree::GNodeType::CutVertex)
    {
        wType matchingWeight=getMatchingValue(vG,aH);
        if (matchingWeight>0)
        {
            weight += getMatchingValue(vG,aH);
            expand = true;
        }
    }
    return weight;
}

// computes an isomorphism of weight weightBBPIso
void BBP_MCSI::computeIsomorphism()
{
    bool nonEmptyIsomorphismPresent=false;
    computeSize(); // In case this did not happen yet
    //cout << "Computing Isomorphism\n";
    //cout << "Weight :"<<m_weightBBPIso<<endl;

    queue<BGxGPair> setSxQ;
    setSxQ.push(BGxGPair(m_b_initial,nullptr));
    node bG,xG,bH,vG,vH,c,xG_,bG_;
    adjEntry adjC,adjB;
    int numNodesbG,numNodesbH;
    while (!setSxQ.empty())
    {
        BGxGPair &bGxGPair=setSxQ.front();
        bG=bGxGPair.bG;
        xG=bGxGPair.xG;
        numNodesbG=m_BC_G->numberOfNodes(bG);
        // Part of Sx, where at least one edge is mapped
        forall_nodes(bH,*m_BCH)
        {
            numNodesbH=m_BC_H->numberOfNodes(bH);
            // beides Brücken oder beides Blöcke
            if ((numNodesbG == 2 && numNodesbH == 2) || (numNodesbG > 2 && numNodesbH > 2))
            {
                if (bbpEdge(bG, bH, xG)==m_weightBBPIso) // in this set there is a BBP-MCSI
                {
                    if (m_enumerate)
                    {
                        m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::BLOCK_BRIDGE_BG,bG,bH,1));
                        //cout << "block: prior to maximize\n";
                        maximizeIsomorphism();
                        //cout << "block: after maximize\n";
                        nonEmptyIsomorphismPresent=true;
                        outputIsomorphism();
                        //cout << "block: output 1 done\n";
                        enumerateIsomorphisms(); // enumerate further isomorphisms
                        //cout << "block: enum done\n";
                    }
                    else
                    {
                        //cout << "BBP-Edge: ";
                        outputIsomorphism(m_LocalIso_BBPE_bGTobH[bG][bH].front());
                        return;
                    }
                }
            }
        }

        // Part of Sx, a single vertex is mapped
        for (node v : m_VbcG[bG])
        {
            vG=m_BC_G->original(v);
            if (xG != vG) // vG not excluded
            {
                forall_nodes(vH,*m_H)
                {
                    if (compatible(vG,vH) && bbpSingleVertex(bG, vG, vH)==m_weightBBPIso)
                    {
                        if (m_enumerate)
                        {
                            // single local mapping, maybe expanded into other B-nodes
                            wType wSV = w(vG,vH);
                            if (m_BC_G->typeOfGNode(vG)== BCTree::GNodeType::CutVertex && m_BC_H->typeOfGNode(vH)== BCTree::GNodeType::CutVertex)
                                wSV += getMatchingValue(vG,m_BC_H->rep(vH));
                            if (wSV == m_weightBBPIso)
                            //if (m_LocalIso_BBPSV_vGTovH[vG][vH]==nullptr)
                            //if (m_LocalIso_BBPSV_vGTovH[vG][vH].empty()==true)
                            {
                                //cout << "sv: vertices "<<vG<<", "<<vH<<" m_EnumMappingSourceNode.size()="<<m_EnumMappingSourceNode.size()<< endl;
                                m_EnumMappingSourceNode.push_back(vG);
                                m_EnumMappingTargetNode.push_back(vH);
                                if (bbpSingleVertex(bG, vG, vH)>=w(vG,vH)) // expansion
                                {
                                    m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::MATCHING,vG,m_BC_H->rep(vH),1));
                                    maximizeIsomorphism();
                                    nonEmptyIsomorphismPresent=true;
                                    outputIsomorphism();
                                    enumerateIsomorphisms(); // enumerate further isomorphisms
                                }
                                /*else
                                {
                                    nonEmptyIsomorphismPresent=true;
                                    outputIsomorphism();
                                    removeEnumMapping(1);
                                }*/
                            }
                            //else
                            if (m_LocalIso_BBPSV_vGTovH[vG][vH] != nullptr) // BBP-Edge of some neighboring block
                            {
                                bH=m_BC_H->bcproper(vH);
                                forall_adj(adjB,m_BC_G->bcproper(vG))
                                {
                                    bG_=adjB->twinNode();
                                    if (bG_ != bG)
                                    {
                                        if (bbpEdge(bG_,bH,nullptr,vG,vH)==m_weightBBPIso)
                                        {
                                            m_EnumMappingSourceNode.push_back(vG);
                                            m_EnumMappingTargetNode.push_back(vH);
                                            m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::BLOCK_BRIDGE_AG,m_BC_G->repVertex(vG,bG_),m_BC_H->repVertex(vH,bH),1));
                                            maximizeIsomorphism();
                                            nonEmptyIsomorphismPresent=true;
                                            outputIsomorphism();
                                            enumerateIsomorphisms();
                                            // removeEnumMapping(1); // Mappings are cleared at the end of enumerateIsomorphisms();
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            // single local mapping, maybe expanded into other B-nodes
                            if (m_LocalIso_BBPSV_vGTovH[vG][vH]==nullptr)
                            //if (m_LocalIso_BBPSV_vGTovH[vG][vH].empty()==true)
                            {
                                vector<Mapping> mapping;
                                mapping.push_back(Mapping(vG,m_BC_H->rep(vH),bbpSingleVertex(bG, vG, vH)>w(vG,vH)?true:false));
                                //cout << "BBP-SV: ";
                                outputIsomorphism(mapping);
                            }
                            else // BBP-Edge of some neighboring block
                            {
                                //cout << "BBP-SV->E: ";
                                m_LocalIso_BBPSV_vGTovH[vG][vH]->emplace_back(Mapping(vG,m_BC_H->rep(vH),false));
                                outputIsomorphism(*m_LocalIso_BBPSV_vGTovH[vG][vH]);
                                m_LocalIso_BBPSV_vGTovH[vG][vH]->pop_back();
                            }
                            return;
                        }
                    }
                }
                // Possible LaWeCSu solution where vG is a skipped inner vertex
                if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE && m_BC_G->typeOfGNode(vG)== BCTree::GNodeType::CutVertex)
                {
                    if (m_LaWeCS_SkippedRootWeight[vG] == m_weightBBPIso)
                    {
                        outputIsomorphism(m_LaWeCS_SkippedRootMapping[vG]);
                        return;
                    }
                }
            }
        }

        // No vertex of bG is mapped, recursion implemented as queue
        forall_adj(adjC,bG)
        {
            c=adjC->twinNode(); // neighboring CV in BC-tree
            xG_=m_BC_G->original(m_VbcG[c][0]); // same node in G
            if (xG != xG_)
            {
                forall_adj(adjB,c)
                {
                    bG_=adjB->twinNode();
                    if (bG_ != bG)
                        setSxQ.push(BGxGPair(bG_,xG_));
                }
            }
        }
        setSxQ.pop();
    }
    // no isomorphism found, therefore display empty mapping
    if (!nonEmptyIsomorphismPresent)
    {
#ifdef GRAPHICS
        if (m_graphical_display)
            display();
#endif
    }
}

// outputs local isomorphism starting in given Mapping
void BBP_MCSI::outputIsomorphism(vector<Mapping>& localMapping)
{
    list<Mapping*> globalMapping;
    // transfer local mapping into a queue
    for (auto it=localMapping.begin();it != localMapping.end(); ++it)
        globalMapping.push_back(&*it);
    // output mapping and enqueue expansion, if indicated
#ifdef GRAPHICS
    int mcs_id=0;
#endif
    cout <<"Weight "<< m_weightBBPIso<<": ";
    while (!globalMapping.empty())
    {
        Mapping& mapping=*globalMapping.front();
        // Output mapping
        cout << "("<<mapping.vG->index()+m_ig_G->fogNodeEdgeOffset <<","<<m_BC_H->original(mapping.aH)->index()+m_ig_H->fogNodeEdgeOffset <<") ";
#ifdef GRAPHICS
        if (m_graphical_display)
        {
            m_graphVisualizationG->setVertexMcsID(mapping.vG,++mcs_id);
            m_graphVisualizationH->setVertexMcsID(m_BC_H->original(mapping.aH),mcs_id);
        }
#endif
        // Enqueue further mappings based on computed matching
        if (mapping.expand)
        {
            //cout << "M ";
            const vector<MWM::AgentJobPair>& matching=getMatching(mapping.vG,mapping.aH);
            //cout << matching.size();
            for (auto itM=matching.begin(); itM != matching.end(); ++itM)
            {
                //cout << "("<<itM->matchingAgent <<","<<itM->matchingJob <<"), ";
                node aG=m_BC_G->repVertex(mapping.vG,itM->matchingAgent);
                //node vG=m_BC_G->original(aG);
                node aH=m_BC_H->repVertex(m_BC_H->original(mapping.aH),itM->matchingJob);
                //cout << " M vG,aG,vH,aH="<<vG<<","<<aG<<","<<m_BC_H->original(aH)<<","<<aH<<" ";
                LaWeCSExpand laWeCSExpand = MWM_EXPANSION; // default to MWM

                //wType wMatching,wSkippedVertex;
                if (aG->degree()==1 && aH->degree()==1) // both chidren are bridges, then check if expansion through skipped vertices
                {
                    node childaH=aH->firstAdj()->twinNode();
                    //node childvH=m_BC_H->original(aH);
                    node childvG=m_BC_G->original(aG->firstAdj()->twinNode());
                    laWeCSExpand = m_LaWeCSExpand[childvG][childaH];
                    if (laWeCSExpand == SKIPPED_VERTEX)
                    {
                        //cout << "LawMapping vG,vH: "<<m_LaWeCSMapping[childvG][childaH].vG<<","<<m_BC_H->original(m_LaWeCSMapping[childvG][childaH].aH);
                        globalMapping.push_front(&m_LaWeCSMapping[childvG][childaH]);
                    }
                }
                //cout << "laW="<<laWeCSExpand<< " ";
                if (laWeCSExpand == MWM_EXPANSION || laWeCSExpand == MWM_AND_SKIPPED_VERTEX)
                {
                    for (auto it=m_LocalIso_BBPE_aGToaH[aG][aH].front().begin();it != m_LocalIso_BBPE_aGToaH[aG][aH].front().end(); ++it)
                    //	if (it->vG != mapping.vG) // do not store given vertex
                        globalMapping.push_front(&*it);
                }
            }
        }
        globalMapping.remove(&mapping);
    }
    cout << endl;
#ifdef GRAPHICS
    if (m_graphical_display)
        display();
#endif
}

void BBP_MCSI::maximizeIsomorphism()
{
    while (!m_EnumStackInProgressToMaximize.empty())
    {
        EnumStackElement &stackElement = m_EnumStackInProgressToMaximize.top();
        m_EnumStackToEnumerate.push(stackElement);
        expansionType expType=stackElement.type;
        node nG=stackElement.nG;
        node nH=stackElement.nH;
        m_EnumStackInProgressToMaximize.pop();
        if(expType==expansionType::BLOCK_BRIDGE_BG)
        {
            //cout << "Push on stack BLOCK_BRIDGE_BG bG="<<nG<<" ,bH="<<nH<<", numel="<<stackElement.num_elements_in_hierarchy<<endl;
            // only one method for both blocks and bridges
            enumerateBBPE_bGTobH(nG, nH);
        }
        else if(expType==expansionType::BLOCK_BRIDGE_AG)
        {
            //cout << "Push on stack BLOCK_BRIDGE_AG vG="<< m_BC_G->original(nG)<<" ,vH="<<m_BC_H->original(nH)<<", numel="<<stackElement.num_elements_in_hierarchy<<endl;
            if (nG->degree()==1) // node of edge has degree 1, any node of a block has degree >1
            {
                enumerateBbpEdgeaG(nG, nH);
            }
            else
            {
                enumerateBBPE_aGToaH(nG, nH);
                // m_LocalIso_BBPE_bGTobH[bG][bH]
                //cout << "Maximize: Transferring / enumerating data from m_LocalIso_BBPE_aGToaH[aG][aH] block not yet implemented\n";
                //exit(1);
            }
        }
        else if(expType==expansionType::MATCHING)
        {
            if (m_BC_G->typeOfGNode(nG) == BCTree::GNodeType::CutVertex && m_BC_H->typeOfGNode(m_BC_H->original(nH)) == BCTree::GNodeType::CutVertex)
            {
                //cout << "Push on stack MATCHING vG="<< nG<<" ,vH="<<m_BC_H->original(nH)<<", numel="<<stackElement.num_elements_in_hierarchy<<endl;
                //unsigned mwmID=m_MWM_ID[nG][nH];
                unsigned mwmID=m_instance_ID[nH];
                //cout << "MWM ID="<<mwmID<<endl;
                MWM& matching=*(m_Matching[nG][m_BC_H->original(nH)]);
                const NodeArray<node>& vMatching_to_bGH= matching.getvEnumTovOriginal(mwmID);
                vector<edge>* matchingEdges=matching.enumNextMatching(mwmID);
                for (edge e:*matchingEdges)
                {
                    node aG=m_BC_G->repVertex(nG,vMatching_to_bGH[e->source()]);
                    node aH=m_BC_H->repVertex(m_BC_H->original(nH),vMatching_to_bGH[e->target()]);
                    m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::BLOCK_BRIDGE_AG,aG,aH,0));
                }
                if (matchingEdges->size()>0)
                    m_EnumStackInProgressToMaximize.top().num_elements_in_hierarchy=static_cast<unsigned int>(matchingEdges->size());
            }
            else
                m_EnumStackToEnumerate.pop();
        }
    }
}

void BBP_MCSI::outputIsomorphism()
{
    const int outputEveryXth=1;
    const int outputAtLeast=10;
#ifdef GRAPHICS
    int mcs_id=0;
#endif
    ++m_mcs_number;
    if (m_mcs_number%outputEveryXth== 0 || m_mcs_number<=outputAtLeast)
        cout << m_mcs_number<<". BBP-MCS with "<<m_EnumMappingSourceNode.size() << " nodes: ";
    for (std::vector<node,allocator<node>>::iterator it_s=m_EnumMappingSourceNode.begin(), it_t=m_EnumMappingTargetNode.begin(); it_s != m_EnumMappingSourceNode.end(); ++it_s, ++it_t)
    {
        if (m_mcs_number%outputEveryXth== 0 || m_mcs_number<=outputAtLeast)
            cout << "("<< it_s.operator *()->index()+m_ig_G->fogNodeEdgeOffset << "," << it_t.operator *()->index()+m_ig_H->fogNodeEdgeOffset<<") ";
#ifdef GRAPHICS
        if (m_graphical_display)
        {
            m_graphVisualizationG->setVertexMcsID(*it_s,++mcs_id);
            m_graphVisualizationH->setVertexMcsID(*it_t,mcs_id);
        }
#endif
    }
    if (m_mcs_number%outputEveryXth== 0 || m_mcs_number<=outputAtLeast)
        cout << endl;
#ifdef GRAPHICS
    if (m_graphical_display)
        display();
#endif
}

void BBP_MCSI::enumerateIsomorphisms()
{
    while(!m_EnumStackToEnumerate.empty())
    {
        EnumStackElement &stackElement = m_EnumStackToEnumerate.top();
        expansionType expType=stackElement.type;
        node nG=stackElement.nG;
        node nH=stackElement.nH;
        unsigned num_elements_in_hierarchy=stackElement.num_elements_in_hierarchy;
        if(expType==expansionType::BLOCK_BRIDGE_BG)
        {
            //cout << "Enum BLOCK_BRIDGE_BG\n";
            // Only one method for both blocks and bridges
            if (enumerateBBPE_bGTobH(nG,nH)) // other local mapping available?
            {
                maximizeIsomorphism();
                outputIsomorphism();
            }
            else
            {
                m_EnumStackToEnumerate.pop(); // BLOCK_BRIDGE_BG is only possible at the beginning. Therefore enumeration with this starting block or bridge complete
            }
        }
        else if(expType==expansionType::BLOCK_BRIDGE_AG)
        {
            //cout << "Get from stack BLOCK_BRIDGE_AG vG="<< m_BC_G->original(nG)<<" ,vH="<<m_BC_H->original(nH)<<", numel="<<stackElement.num_elements_in_hierarchy<<endl;
            if (nG->degree()==1) // bridge (blocks have node degree > 1 for any node)
            {
                // there is always exactly one mapping in each bridge_aG, therefore we cannot enumerate another
                removeEnumMapping(1);
                if (num_elements_in_hierarchy>0) // first block or bridge of matching. Then remove other blocks and bridges of current hierarchy
                {
                    while (num_elements_in_hierarchy>1)
                    {
                        m_EnumStackInProgressToMaximize.pop();
                        --num_elements_in_hierarchy;
                    }
                }
                else // add to maximize stack
                {
                    m_EnumStackInProgressToMaximize.push(stackElement);
                }
                m_EnumStackToEnumerate.pop();
            }
            else
            {
                //cout << "Enum BLOCK_BRIDGE_AG\n";
                //exit(1);
                if (enumerateBBPE_aGToaH(nG,nH)) // other local mapping available?
                {
                    maximizeIsomorphism();
                    outputIsomorphism();
                }
                else
                {
                    if (num_elements_in_hierarchy>0) // first block or bridge of matching. Then remove other blocks and bridges of current hierarchy
                    {
                        while (num_elements_in_hierarchy>1)
                        {
                            m_EnumStackInProgressToMaximize.pop();
                            --num_elements_in_hierarchy;
                        }
                    }
                    else // add to maximize stack
                    {
                        m_EnumStackInProgressToMaximize.push(stackElement);
                    }
                    m_EnumStackToEnumerate.pop(); // enumeration finished -- TODO: check if correct
                }
            }

        }
        else if(expType==expansionType::MATCHING)
        {
            //cout << "Get from stack MATCHING vG="<< nG<<" ,vH="<<m_BC_H->original(nH)<<", numel="<<stackElement.num_elements_in_hierarchy<<endl;
            unsigned mwmID=m_instance_ID[nH];
            MWM& matching=*(m_Matching[nG][m_BC_H->original(nH)]);
            const NodeArray<node>& vMatching_to_bGH= matching.getvEnumTovOriginal(mwmID);
            vector<edge> *matchingEdges=matching.enumNextMatching(mwmID);
            if (matchingEdges==nullptr) // last matching computed
            {
                if (num_elements_in_hierarchy>0) // first matching of block or bridge. Then remove other matchings of current hierarchy
                {
                    while (num_elements_in_hierarchy>1)
                    {
                        m_EnumStackInProgressToMaximize.pop();
                        --num_elements_in_hierarchy;
                    }
                }
                else // add to maximize stack
                {
                    m_EnumStackInProgressToMaximize.push(stackElement);
                }
                m_EnumStackToEnumerate.pop();
            }
            else // next matching computed
            {
                for (edge e:*matchingEdges)
                {
                    node aG=m_BC_G->repVertex(nG,vMatching_to_bGH[e->source()]);
                    node aH=m_BC_H->repVertex(m_BC_H->original(nH),vMatching_to_bGH[e->target()]);
                    m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::BLOCK_BRIDGE_AG,aG,aH,0));
                }
                if (matchingEdges->size()>0)
                    m_EnumStackInProgressToMaximize.top().num_elements_in_hierarchy=static_cast<unsigned>(matchingEdges->size());
                maximizeIsomorphism();
                outputIsomorphism();
            }
        }
    }


#ifdef GRAPHICS
    // clear current mapping
    if (m_graphical_display)
    {
        for (node n:m_EnumMappingSourceNode)
            m_graphVisualizationG->setVertexMcsID(n,0);
        for (node n:m_EnumMappingTargetNode)
            m_graphVisualizationH->setVertexMcsID(n,0);
    }
#endif
    m_EnumMappingSourceNode.clear();
    m_EnumMappingTargetNode.clear();
}

bool BBP_MCSI::enumerateBBPE_bGTobH(const node bG, const node bH)
{
    // NodeArray<NodeArray<forward_list<vector<Mapping>>::iterator>> m_EnumBlockBridgeMappingIterator_bGTobH
    // NodeArray<NodeArray<forward_list<vector<Mapping>>>> m_LocalIso_BBPE_bGTobH;
    // m_LocalIso_BBPE_bGTobH[n][nH].end();
    forward_list<vector<Mapping>>::iterator &mappingIt=m_EnumBlockBridgeMappingIterator_bGTobH[bG][bH];
    // New enumeration from start
    if (!m_LocalIso_BBPE_bGTobH_first[bG][bH])
    {
        m_LocalIso_BBPE_bGTobH_first[bG][bH] = true;
        mappingIt = m_LocalIso_BBPE_bGTobH[bG][bH].end();
    }
    //mappingIt = m_LocalIso_BBPE_bGTobH[bG][bH].begin();
    if (mappingIt == m_LocalIso_BBPE_bGTobH[bG][bH].end())
    {
        mappingIt = m_LocalIso_BBPE_bGTobH[bG][bH].begin();
    }
    // iterate to next stored mapping
    else
    {
        removeEnumMapping(static_cast<int>(mappingIt->size())); // remove previously stored mapping
        ++mappingIt;

        // All possible mappings have been enumerated
        if (mappingIt == m_LocalIso_BBPE_bGTobH[bG][bH].end())
            return false;
    }

    // copy stored mapping into current solution vector
    // emplace matchings on m_EnumStackInProgressToMaximize
    unsigned numMatchingVertices=0;
    for (vector<Mapping>::iterator vertexIt = mappingIt->begin(); vertexIt != mappingIt->end(); ++ vertexIt)
    {
        node vG=vertexIt->vG;
        node aH=vertexIt->aH;
        node vH=m_BC_H->original(aH);
        m_EnumMappingSourceNode.push_back(vG);
        m_EnumMappingTargetNode.push_back(vH);
        if (m_BC_G->typeOfGNode(vG)==BCTree::GNodeType::CutVertex && m_BC_H->typeOfGNode(vH)== BCTree::GNodeType::CutVertex)
        //if (vertexIt->expand)
        {
            m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::MATCHING,vG,aH,0));
            ++numMatchingVertices;
        }
    }
    if (numMatchingVertices > 0)
        m_EnumStackInProgressToMaximize.top().num_elements_in_hierarchy = numMatchingVertices; // total number of matching vertices.

    return true;
}

bool BBP_MCSI::enumerateBBPE_aGToaH(const node aG, const node aH)
{
    forward_list<vector<Mapping>>::iterator &mappingIt=m_EnumBlockBridgeMappingIterator_aGToaH[aG][aH];

    // New enumeration from start
    if (mappingIt == m_LocalIso_BBPE_aGToaH[aG][aH].end())
    {
        mappingIt = m_LocalIso_BBPE_aGToaH[aG][aH].begin();
    }
    // iterate to next stored mapping
    else
    {
        removeEnumMapping(static_cast<int>(mappingIt->size())); // remove previously stored mapping
        ++mappingIt;

        // All possible mappings have been enumerated
        if (mappingIt == m_LocalIso_BBPE_aGToaH[aG][aH].end())
            return false;
    }

    // copy stored mapping into current solution vector
    // emplace matchings on m_EnumStackInProgressToMaximize
    unsigned numMatchingVertices=0;
    for (vector<Mapping>::iterator vertexIt = mappingIt->begin(); vertexIt != mappingIt->end(); ++ vertexIt)
    {
        node vG=vertexIt->vG;
        node aH=vertexIt->aH;
        node vH=m_BC_H->original(aH);
        m_EnumMappingSourceNode.push_back(vG);
        m_EnumMappingTargetNode.push_back(vH);
        if (vertexIt->expand)
        {
            m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::MATCHING,vG,aH,0));
            ++numMatchingVertices;
        }
    }
    if (numMatchingVertices > 0)
        m_EnumStackInProgressToMaximize.top().num_elements_in_hierarchy = numMatchingVertices; // total number of matching vertices.

    return true;
}

void BBP_MCSI::enumerateBbpEdgeaG(const node aG, const node aH)
{
    Mapping &mapping=m_LocalIso_BBPE_aGToaH[aG][aH].front()[0]; // mapping of the two other vertices
    node vG=mapping.vG;
    //cout << " " <<m_EnumMappingSourceNode.size()<< " ";
    m_EnumMappingSourceNode.push_back(vG);
    m_EnumMappingTargetNode.push_back(m_BC_H->original(mapping.aH));
    //cout << " " <<m_EnumMappingSourceNode.size()<< " vG="<<vG<<" vH="<<m_BC_H->original(mapping.aH)<< " ";
    //if (mapping.expand)
        m_EnumStackInProgressToMaximize.emplace(EnumStackElement(expansionType::MATCHING,vG,mapping.aH,1));
}

void BBP_MCSI::enumerateBbpSingleVertex(const node vG, const node vH)
{
    cout << "enum single vertex\n";
}

void BBP_MCSI::removeEnumMapping(int numNodes)
{
    for (int i=0; i<numNodes; ++i)
    {
#ifdef GRAPHICS
        if (m_graphical_display)
        {
            m_graphVisualizationG->setVertexMcsID(m_EnumMappingSourceNode.back(),0);
            m_graphVisualizationH->setVertexMcsID(m_EnumMappingTargetNode.back(),0);
        }
#endif
#ifndef LAWECSU_NDEBUG
        if (m_EnumMappingSourceNode.empty())
        {
            cout << "Tried to pop from empty m_EnumMappingSourceNode\n";
            exit(1);
        }
#endif
        m_EnumMappingSourceNode.pop_back();
        m_EnumMappingTargetNode.pop_back();
/*#ifndef LAWECSU_NDEBUG
        if (m_EnumMappingSourceNode.empty())
        {
            cout << "m_EnumMappingSourceNode now empty\n";
        }
#endif*/
    }
}
// Computes weight of maximum isomorphism on CC(V(G)\{xG}, V(bG))
wType BBP_MCSI::setSx()
{
    queue<BGxGPair> setSxQ;
    setSxQ.push(BGxGPair(m_b_initial,nullptr));
    wType maxWeight=0;
    node bH,bG,xG,vG,vH,c,xG_,bG_;
    int numNodesbG,numNodesbH;
    adjEntry adjC,adjB;
    while (!setSxQ.empty())
    {
        BGxGPair &bGxGPair=setSxQ.front();
        bG=bGxGPair.bG;
        xG=bGxGPair.xG;
        numNodesbG=m_BC_G->numberOfNodes(bG);

        // Part of Sx, where at least one edge is mapped
        forall_nodes(bH,*m_BCH)
        {
            numNodesbH=m_BC_H->numberOfNodes(bH);
            // beides Brücken oder beides Blöcke
            if ((numNodesbG==2 && numNodesbH==2) || (numNodesbG>2 && numNodesbH > 2))
                maxWeight=max(maxWeight,bbpEdge(bG, bH, xG));
        }

        // Part of Sx, a single vertex is mapped
        for (node &v : m_VbcG[bG])
        {
            vG=m_BC_G->original(v);
            if (xG != vG) // vG not excluded
            {
                // MCS
                forall_nodes(vH,*m_H)
                    if (compatible(vG,vH))
                        maxWeight=max(maxWeight,bbpSingleVertex(bG, vG, vH));


                // LaWeCSE; Lemma 10 M_2
                // v is possibly an inner vertex of a topological path
                if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE && m_BC_G->typeOfGNode(vG)== BCTree::GNodeType::CutVertex)
                {
                    wType& laWeCS_SkippedRootWeight = m_LaWeCS_SkippedRootWeight[vG];
                    vector<Mapping>& laWeCS_SkippedRootMapping = m_LaWeCS_SkippedRootMapping[vG];
                    laWeCS_SkippedRootWeight = WEIGHT_NOT_COMPATIBLE;
                    // the children of H are determined through bridges
                    forall_nodes(bH,*m_BCH)
                    {
                        if (m_VbcH[bH].size()==2) // iterate over all bridges of H
                        {
                            // find all bridges connected to v. There must be at least two.
                            node cG=m_BC_G->bcproper(vG);
                            // best and 2nd best weight (children of v) for first and second vertex of bridge bH
                            wType wChildvG0_first = WEIGHT_NOT_COMPATIBLE, wChildvG0_second = WEIGHT_NOT_COMPATIBLE;
                            wType wChildvG1_first = WEIGHT_NOT_COMPATIBLE, wChildvG1_second = WEIGHT_NOT_COMPATIBLE;
                            // corresponding children v
                            node nChildvG0_first = nullptr, nChildvG0_second = nullptr;
                            node nChildvG1_first = nullptr, nChildvG1_second = nullptr;

                            forall_adj(adjB,cG)
                            {
                                bG_=adjB->twinNode();
                                if (bG_ != bG && m_VbcG[bG_].size() == 2) // child of v, which is a bridge
                                {
                                    node childvG = (m_BC_G->original(m_VbcG[bG_][0]) == vG ? m_BC_G->original(m_VbcG[bG_][1]):m_BC_G->original(m_VbcG[bG_][0]));

                                    // first orientation of edge
                                    node childaH0 = m_VbcH[bH][0];
                                    //node childvH0 = m_BC_H->original(childaH0);
                                    wType weightSkipped = getSkippedVertexValue(childvG, childaH0);
                                    bool expand;
                                    wType weightMatching = getWeightPlusCVMatching(childvG, childaH0, expand);
                                    //cout << "vG="<< vG->index()+1<<" cvG="<<childvG->index()+1<<" cvH="<<childvH0->index()+1<<" cvHE="<<m_BC_H->original(m_VbcH[bH][1])->index()+1
                                //			<< " WS="<<weightSkipped<<" WM="<<weightMatching<<endl;
                                    wType wMax = max(weightSkipped,weightMatching);
                                    if (wMax > wChildvG0_first) // new best weight
                                    {
                                        nChildvG0_second = nChildvG0_first;
                                        wChildvG0_second = wChildvG0_first;
                                        nChildvG0_first = childvG;
                                        wChildvG0_first = wMax;
                                    }
                                    else if (wMax > wChildvG0_second) // equal or worse than best, better than second
                                    {
                                        nChildvG0_second = childvG;
                                        wChildvG0_second = wMax;
                                    }

                                    // second orientation of edge
                                    node childaH1 = m_VbcH[bH][1];
                                    //node childvH1 = m_BC_H->original(childaH1);
                                    weightSkipped = getSkippedVertexValue(childvG, childaH1);
                                    weightMatching = getWeightPlusCVMatching(childvG, childaH1, expand);
                                    //cout << "vG="<< vG->index()+1<<" cvG="<<childvG->index()+1<<" cvH="<<childvH1->index()+1<<" cvHE="<<m_BC_H->original(m_VbcH[bH][0])->index()+1
                                    //		<< " WS="<<weightSkipped<<" WM="<<weightMatching<<endl;
                                    wMax = max(weightSkipped,weightMatching);
                                    if (wMax > wChildvG1_first) // new best weight
                                    {
                                        nChildvG1_second = nChildvG1_first;
                                        wChildvG1_second = wChildvG1_first;
                                        nChildvG1_first = childvG;
                                        wChildvG1_first = wMax;
                                    }
                                    else if (wMax > wChildvG1_second) // equal or worse than best, better than second
                                    {
                                        nChildvG1_second = childvG;
                                        wChildvG1_second = wMax;
                                    }

                                    // find best combination of different children
                                    node nBestChildvG0 = nullptr, nBestChildvG1 = nullptr;
                                    wType wLaWeCSbH = WEIGHT_NOT_COMPATIBLE;
                                    if (nChildvG0_first != nChildvG1_first) // top children are different, than this is the best solution
                                    {
                                        nBestChildvG0 = nChildvG0_first;
                                        nBestChildvG1 = nChildvG1_first;
                                        wLaWeCSbH = wChildvG0_first+wChildvG1_first;
                                    }
                                    else // top children are the same. Then best solution is second best of one + best of the other
                                    {
                                        if (wChildvG0_first + wChildvG1_second > wChildvG0_second + wChildvG1_first)
                                        {
                                            nBestChildvG0 = nChildvG0_first;
                                            nBestChildvG1 = nChildvG1_second;
                                            wLaWeCSbH = wChildvG0_first+wChildvG1_second;
                                        }
                                        else
                                        {
                                            nBestChildvG0 = nChildvG0_second;
                                            nBestChildvG1 = nChildvG1_first;
                                            wLaWeCSbH = wChildvG0_second+wChildvG1_first;
                                        }
                                    }
                                    wLaWeCSbH -= m_distancePenalty; // distance penalty for v

                                    // Solution with skipped vertex has greater weight than other solutions
                                    if (wLaWeCSbH > laWeCS_SkippedRootWeight && wLaWeCSbH > 0 && wLaWeCSbH > maxWeight)
                                    {
                                        //cout << "New best inner skipped found: "<<wLaWeCSbH<<" vG="<< vG->index()+1<<" nB0="<<nBestChildvG0->index()+1<<" nB1="<<nBestChildvG1->index()+1<<" edge:"<< m_BC_H->original(m_VbcH[bH][0])->index()+1<<","<<m_BC_H->original(m_VbcH[bH][1])->index()+1<<endl;
                                        //maxWeight = max(maxWeight, wLaWeCSbH);
                                        // laWeCS_SkippedRootWeight = wLaWeCSbH;
                                        laWeCS_SkippedRootWeight = maxWeight = wLaWeCSbH;

                                        // determine from which solution the new best weight came
                                        // first Mapping
                                        weightSkipped = getSkippedVertexValue(nBestChildvG0, childaH0);
                                        weightMatching = getWeightPlusCVMatching(nBestChildvG0, childaH0, expand);
                                        laWeCS_SkippedRootMapping.clear();
                                        if (weightMatching > weightSkipped)
                                        {
                                            laWeCS_SkippedRootMapping.push_back(Mapping(nBestChildvG0, childaH0, expand));
                                        }
                                        else
                                        {
                                            laWeCS_SkippedRootMapping.push_back(m_LaWeCSMapping[nBestChildvG0][childaH0]);
                                        }
                                        // second Mapping
                                        weightSkipped = getSkippedVertexValue(nBestChildvG1, childaH1);
                                        weightMatching = getWeightPlusCVMatching(nBestChildvG1, childaH1, expand);
                                        if (weightMatching > weightSkipped)
                                        {
                                            laWeCS_SkippedRootMapping.push_back(Mapping(nBestChildvG1, childaH1, expand));
                                        }
                                        else
                                        {
                                            laWeCS_SkippedRootMapping.push_back(m_LaWeCSMapping[nBestChildvG1][childaH1]);
                                        }
                                        //m_LaWeCS_SkippedRootWeight[vG][bH] = wLaWeCSbH;
                                    }
                                }
                            }
                        }
                    }
                    // debug out
                    //if (laWeCS_SkippedRootWeight > 0)
                    //	cout << "Best inner skipped vertex for vertex vG="<<vG->index()+1<<": "<<laWeCS_SkippedRootMapping.vG->index()+1<<","<<m_BC_H->original(laWeCS_SkippedRootMapping.aH)->index()+1<<","<<laWeCS_SkippedRootMapping.expand
                    //	<< " / "<< laWeCS_SkippedRootMapping.vG2->index()+1<<","<<m_BC_H->original(laWeCS_SkippedRootMapping.aH2)->index()+1<<","<<laWeCS_SkippedRootMapping.expand2<< endl;
                }
            }
        }

        // No vertex of bG is mapped, recursion (implemented as queue)
        forall_adj(adjC,bG)
        {
            c=adjC->twinNode(); // neighboring CV in BC-tree
            xG_=m_BC_G->original(m_VbcG[c][0]); // same node in G
            if (xG != xG_)
            {
                forall_adj(adjB,c)
                {
                    bG_=adjB->twinNode();
                    if (bG_ != bG)
                        setSxQ.push(BGxGPair(bG_,xG_));
                }
            }
        }

        setSxQ.pop();
    }
    //cout << "Maxweight = "<<maxWeight<<endl;
    return maxWeight;
}

void BBP_MCSI::bbpEdgeBridgebGbH(const node bG, const node bH, const edge eG, const edge eH, bool flipEdges, bool skipvG, bool skipTwinG)
{
    node vG = m_BC_G->original(eG->source());
    node twinG = m_BC_G->original(eG->target());
    node aH = eH->source();
    node twinaH = eH->target();
    if (flipEdges)
    {
//		swap(vH, twinH);
        swap(aH, twinaH);
    }
    node vH = m_BC_H->original(aH);
    node twinH = m_BC_H->original(twinaH);

    wType weight = WEIGHT_NOT_COMPATIBLE;
    Mapping mappingvG(vG, aH, false);
    Mapping mappingTwinG(twinG, twinaH, false);
    bool skippedVertexPossible = true;

    // vertices need to be compatible or skipped
    if ((compatible(vG,vH) || skipvG) && (compatible(twinG,twinH) || skipTwinG))
    {
        weight = 0;

        if (!skipvG && !skipTwinG) // add edge's weight, if no vertex is skipped
        {
            weight = w(m_BC_G->original(eG),m_BC_H->original(eH));
        }

        if (!skipvG) // vG is mapped
        {
            weight += getWeightPlusCVMatching(vG, aH, mappingvG.expand);
        }
        else // skipped vertex
        {
            wType weightSkipped = getSkippedVertexValue(vG, aH);
            if (weightSkipped > 0)
            {
                mappingvG = m_LaWeCSMapping[vG][aH];
                weight += weightSkipped;
            }
            else
            {
                skippedVertexPossible = false;
            }
        }

        if (skippedVertexPossible)
        {
            if (!skipTwinG) // twinG is mapped
            {
                weight += getWeightPlusCVMatching(twinG, twinaH, mappingTwinG.expand);
            }
            else // skipped vertex
            {
                wType weightSkipped = getSkippedVertexValue(twinG, twinaH);
                if (weightSkipped > 0)
                {
                    mappingTwinG = m_LaWeCSMapping[twinG][twinaH];
                    weight += weightSkipped;
                }
                else
                {
                    skippedVertexPossible = false;
                }
            }
        }

        // add solution, if possible
        if (skippedVertexPossible)
        {
            forward_list<vector<Mapping>> &mappingList=m_LocalIso_BBPE_bGTobH[bG][bH];
            wType &storedWeight = m_Weight_BBPE_bGTobH[bG][bH];

            // erase previous solutions if new weight is larger
            if (weight > storedWeight)
            {
                mappingList.clear();
            }

            // add new solution if it is larger; or if it is same weight with enumeration and true MCS
            if (weight > storedWeight || (weight == storedWeight && m_enumerate))
            {
                storedWeight = weight;
                mappingList.push_front(vector<Mapping>());
                mappingList.front().push_back(mappingvG);
                mappingList.front().push_back(mappingTwinG);
            }
        }
    }
}

// Weight of maximum isomorphism phi, where at least one edge of B-Node bG is mapped to B-Node bH,
//  restricted to CC(V(G)\{xG},V(bG)) and phi(vG)=vH
// bG, bH are B-Nodes; xG is a possible excluded vertex in G; vG, vH are vertices of G,H
wType BBP_MCSI::bbpEdge(const node bG, const node bH, const node xG, const node vG, const node vH)
{
    int numNodesbG=m_BC_G->numberOfNodes(bG);
    int numNodesbH=m_BC_H->numberOfNodes(bH);
#ifndef LAWECSU_NDEBUG
    if ((numNodesbG==2 && numNodesbH == 2) || (numNodesbG!=2 && numNodesbH != 2))
    {
        if ((vG == nullptr && m_weightBBPIso==m_Weight_BBPE_bGTobH[bG][bH]) || (vG != nullptr && m_weightBBPIso==m_Weight_BBPE_aGToaH[m_BC_G->repVertex(vG,bG)][m_BC_H->repVertex(vH,bH)]))
        {
            cout << "BBP-Edge between ";
            if (m_BC_H->numberOfNodes(bH)==2)
                cout << "Bridges";
            else
                cout << "Blocks";
            cout << " "<<bG<<","<<bH;
            if (xG != nullptr)
                cout << ", exclude: "<<xG;
            if (vG != nullptr)
                cout << ", mapping: "<<vG<<"->"<<vH;
        }
    }
#endif

    wType *maxWeight;
    node aG=nullptr,aH=nullptr;

    // Check, if weight has already been computed
    if (vG == nullptr) // Lookup in Weight_BBPE_bGTobH
    {
        if ((m_Weight_BBPE_bGTobH[bG][bH]) != WEIGHT_UNDEF)
        {
            return m_Weight_BBPE_bGTobH[bG][bH];
        }
        else
            maxWeight = &m_Weight_BBPE_bGTobH[bG][bH];
    }
    else // Lookup in Weight_BBPE_aGToaH
    {
        aG=m_BC_G->repVertex(vG,bG);
        aH=m_BC_H->repVertex(vH,bH);
        if ((m_Weight_BBPE_aGToaH[aG][aH]) != WEIGHT_UNDEF)
        {
            return m_Weight_BBPE_aGToaH[aG][aH];
        }
        else
            maxWeight = &m_Weight_BBPE_aGToaH[aG][aH];
    }
    // Default to not compatible
    *maxWeight = WEIGHT_NOT_COMPATIBLE;

    // Mapping between two bridges
    if (numNodesbG==2 && numNodesbH==2)
    {
        // cout << "two bridges ";
        edge eG=*m_BC_G->hEdges(bG).begin();
        edge eH=*m_BC_H->hEdges(bH).begin();
        // cout << "eGS:"<<m_BC_G->original(eG)->source()<<" ";
        // cout << "eGT:"<<m_BC_G->original(eG)->target()<<" ";
        // Check, if none of the 2 vertices in bG is excluded and edges are compatible
        if (xG!=m_BC_G->original(m_VbcG[bG][0]) && xG!=m_BC_G->original(m_VbcG[bG][1]))
        {
            // cout << "not excluded ";
            // Mapping of vertices is given
            if (vG!=nullptr)
            {
                if (compatible(m_BC_G->original(eG),m_BC_H->original(eH)))
                {
                    // cout << "vG no nullptr ";
                    // Other vertices of the edge
                    const node twinG=m_BC_G->original(aG->firstAdj()->twinNode());
                    const node twinaH=aH->firstAdj()->twinNode();
                    const node twinH=m_BC_H->original(twinaH);

                    // Check if compatible
                    if (compatible(twinG,twinH) && compatible(vG,vH))
                    {
                        // cout << "compatible ";
                        // local weight
                        *maxWeight = w(m_BC_G->original(eG),m_BC_H->original(eH))+w(twinG,twinH)+w(vG,vH);
                        // m_LocalIso_BBPE_aGToaH[aG][aH].push_back(Mapping(vG,aH,false)); // do not store given vertex
                        // Expansion along twinG,twinH
                        wType weightMatching=0;
                        bool expand=false;
                        if (m_BC_G->typeOfGNode(twinG)== BCTree::GNodeType::CutVertex && m_BC_H->typeOfGNode(twinH)== BCTree::GNodeType::CutVertex)
                        {
                            *maxWeight += (weightMatching=getMatchingValue(twinG,twinaH));
                            if (weightMatching > 0)
                                expand=true;
                        }
                        m_LocalIso_BBPE_aGToaH[aG][aH].push_front(vector<Mapping>());
                        m_LocalIso_BBPE_aGToaH[aG][aH].front().push_back(Mapping(twinG,twinaH,expand));
                    }
                }
            }
            // no mapping given, so try both different mappings and expansions
            else
            {
                // MCS weights
                if (compatible(m_BC_G->original(eG),m_BC_H->original(eH)))
                {
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, false, false, false);
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, true, false, false);
                }
                // skipped vertex weights
                if (m_distancePenalty != WEIGHT_NOT_COMPATIBLE)
                {
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, false, true, false);
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, false, true, true);
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, false, false, true);
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, true, true, false);
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, true, true, true);
                    bbpEdgeBridgebGbH(bG, bH, eG, eH, true, false, true);
                }
            }
        }
    }
    // Mapping between two blocks
    else if (numNodesbG>2 && numNodesbH>2)
    {
        if (vG==nullptr) // no fixed mapping given, therefore compute maximum isomorphism where xG is excluded
        {
            node xSPaG=nullptr;
            if (xG!=nullptr)
                xSPaG=m_aGtoSPaG[m_BC_G->repVertex(xG,bG)];
            //m_mcp2txG[bG][bH] = new MaxComPart2Tree(m_bGAG[bG], m_bHAH[bH], m_SPaGtoaG[bG], m_SPaHtoaH[bH], m_SPeGtoaeG[bG], m_SPeHtoaeH[bH], m_bGToSPQRTree[bG], m_bHToSPQRTree[bH], this, xSPaG, nullptr);
            m_mcp2txG[bG][bH] = new MaxComPart2Tree (*m_ig_G, *m_ig_H, bG, bH, this, xSPaG, nullptr,m_maxCliqueTime);
            //m_mcp2txG[bG][bH] = new MaxComPart2Tree (*m_ig_G, *m_ig_H, bG, bH, this, xSPaG, nullptr, m_bGAG[bG], m_bHAH[bH], m_SPaGtoaG[bG], m_SPaHtoaH[bH], m_SPeGtoaeG[bG], m_SPeHtoaeH[bH], m_bGToSPQRTree[bG], m_bHToSPQRTree[bH]);
            m_unfinishedCliqueComputations += m_mcp2txG[bG][bH]->computeBlock(m_LocalIso_BBPE_bGTobH[bG][bH],*maxWeight,m_enumerate);
        }
        else // fixed mapping vG to vH
        {
            if (compatible(vG,vH))
            {
                //m_mcp2tvG[bG][bH] = new MaxComPart2Tree(m_bGAG[bG], m_bHAH[bH], m_SPaGtoaG[bG], m_SPaHtoaH[bH], m_SPeGtoaeG[bG], m_SPeHtoaeH[bH], m_bGToSPQRTree[bG], m_bHToSPQRTree[bH], this, nullptr, m_aGtoSPaG[m_BC_G->repVertex(vG,bG)]);
                m_mcp2tvG[bG][bH] = new MaxComPart2Tree(*m_ig_G, *m_ig_H, bG, bH, this, nullptr, m_aGtoSPaG[m_BC_G->repVertex(vG,bG)],m_maxCliqueTime);
                //m_mcp2tvG[bG][bH] = new MaxComPart2Tree(*m_ig_G, *m_ig_H, bG, bH, this, nullptr, m_aGtoSPaG[m_BC_G->repVertex(vG,bG)],m_bGAG[bG], m_bHAH[bH], m_SPaGtoaG[bG], m_SPaHtoaH[bH], m_SPeGtoaeG[bG], m_SPeHtoaeH[bH], m_bGToSPQRTree[bG], m_bHToSPQRTree[bH]);
                m_unfinishedCliqueComputations += m_mcp2tvG[bG][bH]->computeFixedMapping(m_LocalIso_BBPE_aGToaH[aG],m_Weight_BBPE_aGToaH[aG],m_enumerate); // compute isomorphisms with mapping from vG to all vertices of bH
                *maxWeight = m_Weight_BBPE_aGToaH[aG][aH];
            }
        }
    }

    // Store computed weight and return it
    /*if (vG == nullptr)
        m_Weight_BBPE_bGTobH[bG][bH] = maxWeight;
    else // Lookup in Weight_BBPE_aGToaH
        m_Weight_BBPE_aGToaH[aG][aH] = maxWeight;*/

    return *maxWeight;
}


// Weight of maximum isomorphism phi, where exactly one vertex of B-Node bG is mapped, phi(vG)=vH
// bG is a B-Node; vG, vH are vertices of G,H
wType BBP_MCSI::bbpSingleVertex(const node bG, const node vG, const node vH)
{
#ifndef LAWECSU_NDEBUG
    if (m_weightBBPIso==m_Weight_BBPSV_vGTovH[vG][vH])
    {
        cout << "BBP-SingleVertex in block "<<bG;
        cout << ", mapping: "<<vG<<"->"<<vH;
    }
#endif

    wType &maxWeight = m_Weight_BBPSV_vGTovH[vG][vH];
    wType bbpeWeight;
    adjEntry adjB;
    node bH,bG_;

    // Check, if weight has already been computed
    if (maxWeight != WEIGHT_UNDEF)
    {
#ifndef LAWECSU_NDEBUG
        if (maxWeight == m_weightBBPIso)
            cout << ", Return value: "<<maxWeight<<endl;
#endif
        return maxWeight;
    }
#ifndef LAWECSU_NDEBUG
    else
    {
        if (maxWeight == m_weightBBPIso)
            cout << endl;
    }
#endif

    // Map vG to vH
    maxWeight=w(vG,vH);

    // Additional vertices only possible, if vG is a CV, and vG,vH are compatible
    if (m_BC_G->typeOfGNode(vG)== BCTree::GNodeType::CutVertex && compatible(vG,vH))
    {
        // vH is a CV, then expand into all blocks except bG
        if (m_BC_H->typeOfGNode(vH)== BCTree::GNodeType::CutVertex)
        {
            maxWeight = w(vG,vH) + getMatchingValue(vG,m_BC_H->rep(vH));
        }
        // vH is no CV, then compute BBP_Edge for all blocks where vG is in, and the single one, where vH is in
        else
        {
            bH=m_BC_H->bcproper(vH);
            forall_adj(adjB,m_BC_G->bcproper(vG))
            {
                bG_=adjB->twinNode();
                if (bG_ != bG)
                {
                    bbpeWeight=bbpEdge(bG_,bH,nullptr,vG,vH);
                    if (bbpeWeight >= maxWeight)
                    {
                        maxWeight=bbpeWeight;
                        m_LocalIso_BBPSV_vGTovH[vG][vH]=&m_LocalIso_BBPE_aGToaH[m_BC_G->repVertex(vG,bG_)][m_BC_H->repVertex(vH,bH)].front();
                    }
                }
            }
        }
    }

    // Store computed weight and return it
    //m_Weight_BBPSV_vGTovH[vG][vH] = maxWeight;

#ifndef LAWECSU_NDEBUG
    cout << "BBP-SingleVertex in block "<<bG;
    cout << ", mapping: "<<vG<<"->"<<vH;
    cout << ", Return value: "<<maxWeight<<endl;
#endif
    return maxWeight;
}

// display the graphs
void BBP_MCSI::display()
{
#ifdef GRAPHICS
    int imageSizex=(m_graphVisualizationG->getSizex()+m_graphVisualizationH->getSizex()+24); // space of 20 between graphs, 2 to left border, 2 to right border
    int imageSizey=(max(m_graphVisualizationG->getSizey(),m_graphVisualizationH->getSizey())+4); // top 2, bottom 2

    //double scaling=m_lastscaling==-1?min(2000.0/imageSizex,1000.0/imageSizey):m_lastscaling;
    double scaling=m_lastscaling==-1?min(800.0/imageSizex,1000.0/imageSizey):m_lastscaling;

    int windowSizex=imageSizex*scaling;
    int windowSizey=imageSizey*scaling;

    m_graphVisualizationG->computeTextureAndCoords(scaling);
    m_graphVisualizationH->computeTextureAndCoords(scaling);
    //computeIsomorphism();
    m_graphVisualizationG->applyColors();
    m_graphVisualizationH->applyColors();

    if (!m_window.isOpen())
    {
        string title="LaWeCSE -  N: toggle node index,  E: toggle edge index,  B: toggle block index,  C: toggle block colors,  L: toggle labels,  SPACE: close view";
        m_window.create(sf::VideoMode(windowSizex, windowSizey), title.c_str(), sf::Style::Titlebar | sf::Style::Close | sf::Style::Resize);
        //m_window.setPosition(sf::Vector2<int>(900,0));
        m_window.setPosition(sf::Vector2<int>(20,200));
        m_window.setVerticalSyncEnabled(true);
    }
    //RenderWindow window(VideoMode(windowSizex, windowSizey), "BBP-MCS -  N: toggle node index,  E: toggle edge index,  B: toggle block index,  C: toggle block colors,  L: toggle labels,  SPACE: close view", sf::Style::Titlebar | sf::Style::Close | sf::Style::Resize);
    //window.setVerticalSyncEnabled(true);
    //window.setPosition(sf::Vector2<int>(900,0));

    m_graphVisualizationG->setTexturePosition(2*scaling,(2+(m_graphVisualizationG->getSizey()<m_graphVisualizationH->getSizey()?(m_graphVisualizationH->getSizey()-m_graphVisualizationG->getSizey())/2:0))*scaling);
    m_graphVisualizationH->setTexturePosition((22+m_graphVisualizationG->getSizex())*scaling,(2+(m_graphVisualizationG->getSizey()>m_graphVisualizationH->getSizey()?(m_graphVisualizationG->getSizey()-m_graphVisualizationH->getSizey())/2:0))*scaling);

    m_graphVisualizationG->drawAll();
    m_graphVisualizationH->drawAll();

    bool update=false;
    m_lastscaling=scaling;
    sf::Clock recomputeClock;
    while (m_window.isOpen())
    {
        windowSizex=m_window.getSize().x;
        windowSizey=m_window.getSize().y;
        Event event;
        while (m_window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                m_lastscaling = -1; // compute new size  on next call
                m_window.close();
            }
            else if (event.type == sf::Event::Resized)
            {
                windowSizex=event.size.width;
                windowSizey=event.size.height;
                scaling=min(1.0*event.size.width/imageSizex,1.0*event.size.height/imageSizey);
            }
            else if (event.type ==  sf::Event::KeyPressed)
            {
                if (event.key.code == sf::Keyboard::N)
                {
                    m_graphVisualizationG->toogleShowVertices();
                    m_graphVisualizationH->toogleShowVertices();
                    update=true;
                }
                else if (event.key.code == sf::Keyboard::B)
                {
                    m_graphVisualizationG->toogleShowBlocks();
                    m_graphVisualizationH->toogleShowBlocks();
                    update=true;
                }
                else if (event.key.code == sf::Keyboard::C)
                {
                    m_graphVisualizationG->toogleBlockColors();
                    m_graphVisualizationH->toogleBlockColors();
                    update=true;
                }
                else if (event.key.code == sf::Keyboard::E)
                {
                    m_graphVisualizationG->toogleShowEdges();
                    m_graphVisualizationH->toogleShowEdges();
                    update=true;
                }
                else if (event.key.code == sf::Keyboard::L)
                {
                    m_graphVisualizationG->toogleLabels();
                    m_graphVisualizationH->toogleLabels();
                    update=true;
                }
                else if (event.key.code == sf::Keyboard::Space)
                {
                    return;
                }
            }
        }
        if (scaling == m_lastscaling)
            recomputeClock.restart();
        if (recomputeClock.getElapsedTime().asMilliseconds()>500 || update) // update window at most every 500 ms
        {
            update=false;
            m_lastscaling=scaling;
            m_graphVisualizationG->computeTextureAndCoords(scaling);
            m_graphVisualizationH->computeTextureAndCoords(scaling);
            m_graphVisualizationG->setTexturePosition(2*scaling,(2+(m_graphVisualizationG->getSizey()<m_graphVisualizationH->getSizey()?(m_graphVisualizationH->getSizey()-m_graphVisualizationG->getSizey())/2:0))*scaling);
            m_graphVisualizationH->setTexturePosition((22+m_graphVisualizationG->getSizex())*scaling,(2+(m_graphVisualizationG->getSizey()>m_graphVisualizationH->getSizey()?(m_graphVisualizationG->getSizey()-m_graphVisualizationH->getSizey())/2:0))*scaling);
            m_graphVisualizationG->applyColors();
            m_graphVisualizationH->applyColors();
            m_graphVisualizationG->drawAll();
            m_graphVisualizationH->drawAll();
            recomputeClock.restart();
        }
        m_window.setView(View(FloatRect(0,0, windowSizex, windowSizey)));
        m_window.clear(sf::Color::White);
        m_window.draw(m_graphVisualizationG->getTexture());
        m_window.draw(m_graphVisualizationH->getTexture());
        m_window.display();
    }
#endif
}
