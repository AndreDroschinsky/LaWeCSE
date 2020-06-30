// todo: alle Vergleich mit >, < unter Einbeziehung von EQDIST
/*
 * main.cpp
 *
 *  Created on: 22.06.2015
 *      Author: Andre Droschinsky
 */
#define MAIN
#include <global.h>
#include <TreeGen.h>
#pragma warning(push, 0)
#include <ogdf/decomposition/BCTree.h>
#pragma warning(pop)
//#include <test.h>
#include <random>
#include <bbpmcsi.h>
#include <auxiliary.h>

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

//#include <ogdf/planarity/PlanarizationGridLayout.h>

#ifdef _DEBUG
    #include <assert.h>
#endif

using namespace ogdf;
using namespace std;



int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = static_cast<int>(strlen(line));
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


void foo(int *& a)
{
    a=new int[10];
    for (int i=0; i<10; ++i)
        a[i]=i;
}

int main(int argc, char *argv[]) {
    cout.precision(15);
/*cout << "thread started\n";
    thread nothing(donothing);
nothing.join();
cout << "thread joined\n";
return 0;
    int *a;
    foo(a);
    for (int i=0; i<10; ++i)
        cout << a[i]<< " ";
    return 0;




    cpu_set_t mask;
    int status;

    CPU_ZERO(&mask);
    CPU_SET(0, &mask);
    status = sched_setaffinity(0, sizeof(mask), &mask);
    if (status != 0)
    {
        perror("sched_setaffinity");
    }*/
//createRandomTreeFile(10,5,5,1,1000,false);
//createFogStarGraphFile(1,201,10);
    cout.setf(ios::unitbuf);
    clock_t rndstart=clock();
    //rndstart=1000;
    //setSeed(rndstart);
    srand((unsigned int)rndstart);

    /*// previously used
    bool showComputation=true;
    bool showSummary=false;
    size_t totalRuns=1;
    bool random=false;
    bool r_half_split=false;
    bool r_same_size=false;
    bool checkOuterplanarity=false;
    int checkPlanarity=0;
    bool fullCompare=false;
    bool firstToAll=false;
    vector<InputGraph*> inputGraph;
    vector<vector<size_t>*> graphBySize;
    bool fogOutputBC=false;
    */

    // configuration, some parameters can be specified with command line
    string dbFilename="";
    string db2Filename="";
    ComputationType computationType = ComputationType::COMP_NONE;
    int argument1=0;
    int argument2=0;
    bool argument3=false; // for enumeration with -c
    bool argument4=false; // for graphical display with -c
    string seqFilename="";
    string comparedGraphsFilename="";
    string fogOutputFilename="";
    size_t maxMCSInMemory=1;
    size_t readQueueSize=200;
    int maxCliqueTime=-1;
    bool appendCliqueComputationTimeouts=false;
    int computationThreads=0;
    bool verbose=false; // true for console output, if output into file is active
    bool removeMultiEdgesAndSelfLoops=true;
    LabelFunction labelFunction;
    string labelFunctionFile;
    SimCoefficient simCoefficient= SimCoefficient::SC_WallisEtAl;
    wType distancePenalty=WEIGHT_NOT_COMPATIBLE;

    if (argc==1)
    {
        cout << "No parameters specified. Type "<<argv[0]<<" --help for a list of parameters.\n";
        return 0;
    }
    for (int i=1; i<argc; ++i)
    {
        string argument=argv[i];
        if (argument=="-crtf" && i+6<argc) // database file containing the graphs
        {
            int param1=stoi(argv[++i]); // numTrees
            int param2=stoi(argv[++i]); // minsize
            int param3=stoi(argv[++i]); // maxsize
            int param4=stoi(argv[++i]); // numlabels
            int param5=stoi(argv[++i]); // maxdegree
            int param6=stoi(argv[++i]); // twoEvenSizedSets
            createRandomTreeFile(param1,param2,param3, param4,param5,param6);
        }
        else if (argument=="-crof" && i+9<argc) // database file containing the graphs
        {
            int param1=stoi(argv[++i]);
            int param2=stoi(argv[++i]);
            int param3=stoi(argv[++i]);
            int param4=stoi(argv[++i]);
            int param5=stoi(argv[++i]);
            int param6=stoi(argv[++i]);
            int param7=stoi(argv[++i]);
            int param8=stoi(argv[++i]);
            int param9=stoi(argv[++i]);
            createRandomOuterplanarFile(param1,param2,param3, param4,param5,param6, param7,param8,param9);
        }
        else if (argument=="--help") // database file containing the graphs
        {
            cout << "This program computes block and bridge preserving maximum common subgraphs (BBP-MCS) between input graphs.\n"
                 << "To compute the size of an BBPMCS pairwise between all graphs from two input files, specify...\n"
                 << "-i <PatternFile> <MolecularDatabase>  Search pattern(s) and molecular database. Accepts .GML and .fog inputs.\n"
                 << "-t <N> (optional)                     Computation threads. N=0 for single thread (default); N=k for k computation and one additional IO thread.\n"
                 << "-v (optional)                         Console output if output file is specified.\n"
                 << "Example: "<<argv[0]<<" -i coffein.gml allthegraphs.gml -o results -t 3\n"
                 << "WARNING: <PatternFile> is completely loaded and processed into the main memory. Make sure it is of reasonable size.\n"<<endl

                 << "To compute a single BBPMCS isomorphism, specify...\n"
                 << "-c <GraphFile> <ID1> <ID2>  Graph database (.gml or .fog) and positions in the file, starting from 1.\n"
                 << "-e (optional)               Enumerate all BBP-MCSs.\n"
#ifdef GRAPHICS
                 << "-d (optional)               Graphical output.\n"
#endif
                 << "Example: "<<argv[0]<<" -c somegraphs.gml 1 2\n"<<endl

                 << "To compare repeatedly one graph to the following N graphs, specify...\n"
                 << "-oneToX <GraphFile> <N>  Graph database (.gml or .fog) and parameter N > 0.\n"
                 << "Example: "<<argv[0]<<" -oneToX somegraphs.gml 3\n"
                 << " (this compares 1st to 2nd, 1st to 3rd, 1st to 4th; 5th to 6th,..., 5th to 8th; 9th to 10th,...)\n"<<endl

                 << "OPTIONAL parameters.\n"
                 << "-o <outputFile>  Raw output into file. Columns of this file as follows.\n"
                 << "       Position and label of pattern, position and label of molecular graph, |pattern|, |Mol. Gr.|, similarity.\n"
                 << "-l <LabelFile>   Specify the weight function on the labels.\n"
                 << "-ct <T>          Maximum time per clique mcs computation in ms.\n"
                 << "-cto             Append number of clique computations exceeding time limit in output.\n"
                 << "-sc <1--8>       Similarity Coefficient (Table 2 of Effectiveness of graph-based and fingerprint-based similarity measures for virtual screening of 2D chemical structure databases).\n"
                 << "-p <D>           Distance penalty <double> for skipped vertices. Default: infinity, i.e., no skipped vertices. Incompatible with -e.\n"
                // << "-rms             Remove parallel (multi) edges and selfloops.\n"
                 << "\nWARNING: This software is beta and work in progress. Use at your own risk.\n";
            return 0;
        }
#ifdef GRAPHICS
        else if (argument=="-d")
            argument4=true;
#else
        else if (argument=="-d")
        {
            cerr << "Graphical output has not been compiled. Add #define GRAPHICS in global.h\n";
        }
#endif

        else if (argument=="-t" && i+1<argc) // database file containing the graphs
        {
            try
            {
                int t=stoi(argv[++i]);
                if (t>256 || t<0)
                {
                    cerr << "Warning: Number of computation threads must be between 0 and 256\n";
                }
                else
                {
                    computationThreads = t;
                }
            }
            catch(...)
            {
                cerr << "Please provide the a number between 0 and 63 when using -t.\n";
            }
        }
        else if (argument=="-mb" && i+1<argc) // maximum number of bbpmcs computations in memory before delete occurs
        {
            maxMCSInMemory=stoi(argv[++i]);
            if (maxMCSInMemory<1)
            {
                maxMCSInMemory=1;
            }
            //maxMCSInMemory_specified = true;
        }
        else if (argument=="-rq" && i+1<argc) // maximum number of graphs in read ahead queue
        {
            try
            {
                int rQS=stoi(argv[++i]);
                if (rQS > 0)
                {
                    readQueueSize = rQS;
                }
            }
            catch(...)
            {

            }
        }
        else if (argument=="-v") // database file containing the graphs
        {
            verbose=true;
        }
        else if (argument=="-c" && i+3<argc) // database and ids to compare
        {
            dbFilename=argv[++i];
            try
            {
                int arg=stoi(argv[++i])-1;
                if (arg < 0)
                {
                    cerr << "Graph ID must be at least 1.\n";
                }
                else
                {
                    argument1 = arg;
                }
            }
            catch (...)
            {
                cerr << "ID must be a natural number.\n";
            }

            try
            {
                int arg = stoi(argv[++i])-1;
                if (arg < 0)
                {
                    cerr << "Graph ID must be at least 1.\n";
                }
                else
                {
                    argument2 = arg;
                }
            }
            catch (...)
            {
                cerr << "ID must be a natural number.\n";
            }

            computationType = ComputationType::COMP_C;
        }
        else if (argument=="-i" && i+2<argc) // 2 input files to be compared pairwise
        {
            dbFilename=argv[++i];
            db2Filename=argv[++i];
            computationType = ComputationType:: COMP_I1;
        }
        else if (argument=="-e") // enumerate solutions
        {
            argument3=true;
            if (distancePenalty != WEIGHT_NOT_COMPATIBLE)
            {
                cerr << "Enumeration is not compatible with distance penalty. Switched to MCS computation.\n";
                distancePenalty = WEIGHT_NOT_COMPATIBLE;
            }
        }
        else if (argument=="-l" && i+1<argc) // label function file
        {
            labelFunctionFile=argv[++i];
        }
        else if (argument=="-ct" && i+1<argc) // label function file
        {
            try
            {
                int ct=stoi(argv[++i]);
                if (ct < 0)
                {
                    cerr << "Maximum clique computation time must be nonnegative.\n";
                }
                else
                {
                    maxCliqueTime=ct;
                }
            }
            catch(...)
            {
                cerr << "Maximum clique computation time must be a natural number.\n";
            }
        }
        else if (argument=="-cto") // label function file
        {
            appendCliqueComputationTimeouts=true;
        }
        else if (argument=="-sc" && i+1<argc) // label function file
        {
            try
            {
                int sc=stoi(argv[++i]);
                if (sc < 1 || sc > 8)
                {
                    cerr << "Similarity Coefficient must be a natural number between 1 and 8.\n";
                }
                else
                    simCoefficient=(SimCoefficient) sc;
            }
            catch(...)
            {
                cerr << "Similarity Coefficient must be a natural number between 1 and 8.\n";
            }
        }
        /*else if (argument=="-rms") // remove parallel edges
        {
            removeMultiEdgesAndSelfLoops=true;
        }*/
        else if (argument=="-o" && i+1<argc) // output the compared graphs into specified file
        {
            comparedGraphsFilename=argv[++i];
        }
        else if (argument=="-oneToX" && i+2<argc) // compare blockwise one graph to X others in the same file from top to bottom
        {
            dbFilename=argv[++i];
            try
            {
                int arg=stoi(argv[++i]);
                if (arg<1)
                {
                    cerr << "Argument of oneToX must be at least 1"<<endl;
                }
                else
                {
                    argument1 = arg;
                    computationType = ComputationType::COMP_OneToX;
                }
            }
            catch(...)
            {
                cerr << "Argument of oneToX must be a natural number\n";
            }
        }
        else if (argument=="-p" && i+1<argc) // distance penalty
        {
            try
            {
                wType p = stod(argv[++i]);
                if (p<0)
                {
                    cerr << "Distance penalty must be non negative."<<endl;
                }
                else
                {
                    distancePenalty=p;
                    if (argument3)
                    {
                        cerr << "Enumeration is not compatible with distance penalty. Deactivated Enumeration.\n";
                        argument3 = false;
                    }
                }
            }
            catch(...)
            {
                cerr << "Distance penalty must be non negative value.\n";
            }
        }
        /*else if (argument=="-fogoutput" && i+1<argc) // fog format output file
        {
            fogOutputFilename=argv[++i];
        }
        else if (argument=="-fogoutputBC" && i+1<argc) // fog format output file
        {
            fogOutputFilename=argv[++i];
            fogOutputBC=true;
        }
        else if (argument=="-inputsequence" && i+1<argc) // read in sequence of graph numbers to be compared
        {
            seqFilename=argv[++i];
        }
        else if (argument=="-transform" && i+3<argc) // read in sequence of graph numbers to be compared
        {
            string input1=argv[++i];
            string input2=argv[++i];
            string output=argv[++i];
            transform(input1.c_str(),input2.c_str(),output.c_str());
        }*/
        else
        {
            cerr << "Unrecognized parameter "<<argument<<" or too few arguments.\n";
        }
    }

    if (computationType== ComputationType::COMP_NONE)
    {
        cerr << "Please give parameter -i, -c, or -oneToX."<<endl;
        return 1;
    }


    map<string,labelType> stringLabelToSimpleLabel;
    vector<string> SimpleLabelToString;

    // Read weight function file
    if (labelFunctionFile != "")
    {
        readLableFile(labelFunctionFile,labelFunction,stringLabelToSimpleLabel,SimpleLabelToString);
    }

    // Initialize graph reader class
    ReadGraphDB readGraphDB(labelFunction,stringLabelToSimpleLabel,SimpleLabelToString,static_cast<int>(readQueueSize),
            fogOutputFilename.c_str(),comparedGraphsFilename==""?nullptr:comparedGraphsFilename.c_str(),
            computationThreads,maxMCSInMemory, maxCliqueTime, distancePenalty, removeMultiEdgesAndSelfLoops,appendCliqueComputationTimeouts,verbose,
            computationType,argument1,argument2,argument3,argument4,simCoefficient);

    if (computationType == ComputationType::COMP_OneToX || computationType == ComputationType::COMP_C)
    {
        return readGraphDB.readFile(dbFilename.c_str());
    }

    if (computationType == ComputationType::COMP_I1)
    {
        if (readGraphDB.readFile(dbFilename.c_str()) == 1)
            return 1;
        return readGraphDB.readFile(db2Filename.c_str());
    }
    return 0;

    // no more code in main!






    /*InputGraph minIGraph;
    Graph &minGraph=minIGraph.graph;
    node mGn1=minGraph.newNode();
    node mGn2=minGraph.newNode();
    node mGn3=minGraph.newNode();
    minGraph.newEdge(mGn1,mGn2);
    minGraph.newEdge(mGn1,mGn3);
    minGraph.newEdge(mGn2,mGn3);*/


    // Test for Enumeration of Matchings
        /*node a0=testGraph.newNode();
        node a1=testGraph.newNode();
        node a2=testGraph.newNode();
        //node a3=testGraph.newNode();
        //node a4=testGraph.newNode();
        node j0=testGraph.newNode();
        node j1=testGraph.newNode();
        node j2=testGraph.newNode();
        //node j3=testGraph.newNode();
        //node j4=testGraph.newNode();
        testMatching.addAgent(a0);
        testMatching.addAgent(a1);
        testMatching.addAgent(a2);
        //testMatching.addAgent(a3);
        //testMatching.addAgent(a4);
        testMatching.addJob(j0);
        testMatching.addJob(j1);
        testMatching.addJob(j2);
        //testMatching.addJob(j3);
        //testMatching.addJob(j4);
        testMatching.addGraphCopy();
        testMatching.addEdge(0,0,1);
        testMatching.addEdge(0,1,1);
        testMatching.addEdge(1,0,1);
        testMatching.addEdge(1,1,1);
        testMatching.addEdge(2,2,1);
        testMatching.addEdge(2,1,1);
        testMatching.addEdge(1,2,1);
        testMatching.addEdge(0,2,1);
        testMatching.addEdge(2,0,1);
        testMatching.addEdge(0,3,1);
        testMatching.addEdge(1,3,1);
        testMatching.addEdge(2,3,1);
        testMatching.addEdge(3,3,1);
        testMatching.addEdge(3,2,1);
        testMatching.addEdge(3,1,1);
        testMatching.addEdge(3,0,1);

        testMatching.addEdge(1,2,1);
        testMatching.addEdge(2,3,1);
        testMatching.addEdge(3,2,1);
        testMatching.addEdge(3,3,1);
        testMatching.addEdge(1,2,3);
        testMatching.addEdge(3,4,1);
        testMatching.addEdge(4,4,1);
        testMatching.addEdge(2,4,2);
        int ctr=0;
        if (false && withEnumeration) for (unsigned memoryleakt=0; memoryleakt<totalRuns; ++memoryleakt)
        {
            MWM testMatching(withEnumeration);
            Graph testGraph;

            int sizeA=rand()%10+3, sizeJ=rand()%10+3;
            //sizeA=5, sizeJ=5;
            vector<node> agents;
            vector<node> jobs;
            for (int i=0; i<sizeA; ++i)
            {
                agents.push_back(testGraph.newNode());
                testMatching.addAgent(agents[i]);
            }
            for (int j=0; j<sizeJ; ++j)
            {
                jobs.push_back(testGraph.newNode());
                testMatching.addJob(jobs[j]);
            }
            testMatching.addGraphCopy();
            for (int i=0; i<sizeA; ++i)
                for (int j=0; j<sizeJ; ++j)
                    testMatching.addEdge(i,j,rand()%5-1);

            testMatching.compute();
            cout << "computed Matchings"<<endl;

            clock_t enumMatchingTime=clock();
            for (int instanz=sizeJ; instanz<sizeJ+1; ++instanz)
            {
                cout << "Starting instance "<<instanz<<endl;
                vector<edge>* enumMatching=testMatching.enumNextMatching(instanz);
                while (enumMatching != nullptr)
                {
                    ++ctr;
                    enumMatching=testMatching.enumNextMatching(instanz);
                }

            }
            showtime(clock()-enumMatchingTime);
            if ((memoryleakt+1)%500==0)
            {
                cout << "Mem used:"<<getValue()<<" kB";
                cout << "  Total matchings enumerated:"<<ctr<<endl;
                ctr=0;
            }
            cout << "Enumeration finished.\n";
            return 0;
        }*/
        /*Graph G,H;
        G.newNode();
        H.newNode();
        for (int i=0; i<5; ++i) // 65
        {
            node g=G.chooseNode();
            G.newEdge(g,G.newNode());
            node h=H.chooseNode();
            H.newEdge(h,H.newNode());
        }
        node g=G.newNode();
        for (int i=0; i<60; ++i)
            G.newEdge(g,G.newNode());
        node h=H.newNode();
        for (int i=0; i<60; ++i)
            H.newEdge(h,H.newNode());*/


/*
    ofstream comparedGraphs;
    if (comparedGraphsFilename!="")
    {
        comparedGraphs.open(comparedGraphsFilename.c_str(),ios_base::out);
        if (!comparedGraphs.is_open())
            cout << "Could not open compared graphs file "<<comparedGraphsFilename<<" for output."<<endl;
        else
        {
            comparedGraphs << "BBPMCS: Input file "<<dbFilename<<" - Next line: number of computations. Other lines: id1, id2, |G|, |H|, |MCS|, time in microseconds";
            if (firstToAll)
                comparedGraphs << ", similarity";
            comparedGraphs << endl<< totalRuns<<endl;
        }
    }
    // input sequence
    ifstream inputSeq;
    stringstream sstr;
    char line[256];
    if (seqFilename!="")
    {
        inputSeq.open(seqFilename.c_str(),ios_base::in);
        if (!inputSeq.is_open())
            cout << "Could not open sequence file "<<seqFilename<<" for input."<<endl;
        else
        {
            inputSeq.getline(line,255); // some information, not needed here
            inputSeq.getline(line,255); // number of lines
            sstr<< line<<endl;
            sstr>> totalRuns;
        }
    }

    wType weight;
    size_t graphIndex=inputGraph.size();

    if (fullCompare) // compare all graphs
    {
        totalRuns=graphIndex*graphIndex;
    }

    if (firstToAll)
    {
        totalRuns *= graphIndex;
    }





    vector<BBP_MCSI*> bbp;
    //BBP_MCSI *bbpmalloc = new BBP_MCSI(labelFunction,minIGraph,minIGraph);
    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
//	int seed=2334;
//	rng.seed(seed);
    uniform_int_distribution<> dis(0,graphIndex-1); // guaranteed unbiased
    uniform_int_distribution<> dis1stHalf(0,(graphIndex-1)/2); // first half of graphs
    uniform_int_distribution<> dis2ndHalf((graphIndex+1)/2,graphIndex-1); // second half

    cout<<endl;

    //if (!showComputation)
    //	cout << "Computing isomorphisms .. ";
    unsigned summaryindex;
    clock_t totaltimestart=clock(),timestart,time,deletetime=0,totaltime;
    map<unsigned,ComputationTime> computationTime; //unordered_
    // actual computation
    for (unsigned run=0; run<totalRuns; ++run)
    {
        // Initialisierung der BC-Bäume
        if (inputSeq.is_open())
        {
            inputSeq.getline(line,255);
            sstr << line;
            sstr>> graphIndex1;
            sstr>> graphIndex2;
            sstr.str("");
        }
        else if (random)
        {
            if (r_half_split)
            {
                graphIndex1=dis1stHalf(rng);
                graphIndex2=dis2ndHalf(rng);
            }
            else
            {
                graphIndex1=dis(rng);
                if (r_same_size) // second graph randomly chosen among graphs of same size
                    graphIndex2=graphBySize[inputGraph[graphIndex1]->graph.numberOfNodes()]->at(uniform_int_distribution<>(0,graphBySize[inputGraph[graphIndex1]->graph.numberOfNodes()]->size()-1)(rng));
                else
                    graphIndex2=dis(rng);
            }
        }
        if (graphIndex1 >= inputGraph.size() || graphIndex2 >= inputGraph.size() || graphIndex1 < 0 || graphIndex2 < 0)
        {
            cout << "GraphID out of range. ID1="<<graphIndex1+1<<" ID2="<<graphIndex2+1<<" Database size: "<<inputGraph.size()<<" Terminating program."<<endl;
            totalRuns=run;
            break;
        }
        Graph *graph1=&inputGraph[graphIndex1]->graph;
        Graph *graph2=&inputGraph[graphIndex2]->graph;
        summaryindex=(graph1->numberOfNodes()+graph2->numberOfNodes())/2;
        if (showComputation)
            //cout << "Starting computation of BBP-MCS between graph G="<< graphIndex1<< " ("<<graphIndexToDBIndex[graphIndex1]<< ". entry in DB) and graph H="<<graphIndex2<<" ("<<graphIndexToDBIndex[graphIndex2] << ". entry in DB)"
            cout << "Starting computation of "<<run+1<<". BBP-MCS between graph G="<< inputGraph[graphIndex1]->graphLabel<<"("<<inputGraph[graphIndex1]->graphIndexToDBIndex+1 << ")" << " and graph H="
              << inputGraph[graphIndex2]->graphLabel<<"("<<inputGraph[graphIndex2]->graphIndexToDBIndex+1 << ")" << " with |G|="<< graph1->numberOfNodes()<< ", |H|="<< graph2->numberOfNodes()<<" .. ";
        timestart=clock();

        bbp.push_back(new BBP_MCSI(labelFunction,*inputGraph[graphIndex1],*inputGraph[graphIndex2],withEnumeration, graphical_display,&SimpleLabelToString));
        weight=bbp[run]->computeSize();
        time=clock()-timestart;

        if (comparedGraphs.is_open())
        {
            if (!fullCompare)
            {
                comparedGraphs << inputGraph[graphIndex1]->graphIndexToDBIndex << " "<< inputGraph[graphIndex2]->graphIndexToDBIndex<<" "<<graph1->numberOfNodes()+graph1->numberOfEdges()
                        <<" "<<graph2->numberOfNodes()+graph2->numberOfEdges()<<" "<<weight<<" "<<time;
                if (firstToAll)
                {
                    comparedGraphs << " " <<weight/(graph1->numberOfNodes()+graph2->numberOfNodes()+graph1->numberOfEdges()+graph2->numberOfEdges()-weight);
                    comparedGraphs << " " <<inputGraph[graphIndex2]->graphLabel;
                }
                comparedGraphs << endl;
            }
            else // for full compare build quadratic matrix with sizes only
                comparedGraphs <<weight;
        }

        auto it= computationTime.find(summaryindex);
        if (it == computationTime.end())
            computationTime.insert(pair<unsigned,ComputationTime> (summaryindex,ComputationTime(time)));
        else
            it->second.addTime(time);
        if (showComputation)
        {
            cout << "done.\nMaximum BBP-MCSI has a weight of " << weight<< ". Runtime: ";
            showtime(time);
        }
        //if (graphical_display || withEnumeration)
            bbp[run]->computeIsomorphism();
        if (run%maxMCSInMemory==maxMCSInMemory-1)
        {
            timestart=clock();
            for (size_t r=run-maxMCSInMemory+1; r<=run; ++r)
                delete bbp[r];
            //delete bbpmalloc;
            //bbpmalloc = new BBP_MCSI(labelFunction,minIGraph,minIGraph);
            deletetime+=clock()-timestart;
            cout << (run+1)/maxMCSInMemory << ") Memory currently allocated: "<<getValue()<< " kB\n";
        }
        if (fullCompare) // Compare
        {
            ++graphIndex2;
            if (graphIndex2 == graphIndex)
            {
                ++graphIndex1;
                graphIndex2=0;
                if (comparedGraphs.is_open())
                    comparedGraphs<<endl;
            }
            else if (comparedGraphs.is_open())
                comparedGraphs<<", ";
        }
        if (firstToAll) // first to all
        {
            ++graphIndex2;
            if (graphIndex2==graphIndex)
                graphIndex2=0;
        }
    }
    totaltime=clock()-totaltimestart;
    if (!showComputation)
        cout << "done.\n";
    cout << "Total computation time for "<< totalRuns << " BBP-MCSI computations: ";
    showtime(totaltime);
    //delete bbpmalloc;
    cout << "Deleting BBP-Instances .. ";
    timestart=clock();
    for (size_t r=totalRuns/maxMCSInMemory*maxMCSInMemory; r<totalRuns; ++r)
    {
        delete bbp[r];
    }
    bbp.clear();
    deletetime+=clock()-timestart;
    cout << "done. Time: ";
    showtime(deletetime);
    if (showSummary)
    {
        for (auto it = computationTime.begin(); it != computationTime.end(); ++it )
        {
            cout << "Average computation time for "<< it->second.instances<<" instance(s) of average size "<<it->first<< ": ";
            showtime(it->second.average());
        }
    }

    for (size_t i=0; i<graphIndex; ++i)
    {
        delete inputGraph[i];
    }
    if (graphBySizeptr!=nullptr)
        for (size_t i=0; i<graphBySize.size(); ++i)
            delete graphBySize[i];
    if (comparedGraphs.is_open())
        comparedGraphs.close();
    cout << "Terminated"<<endl;
    //	cout << "Done. RND-Init was: "<<rndstart<<endl;
    return 0;*/
}
