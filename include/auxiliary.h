#ifndef AUXILIARY_H
#define AUXILIARY_H
#include <global.h>
#include <bbpmcsi.h>
#pragma warning(push, 0)
#include <ogdf/fileformats/GraphIO.h>
#pragma warning(pop)
#include <thread>
#include <mutex>
#include <condition_variable>
#include <list>
#pragma warning(push, 0)
#include <ogdf/basic/Thread.h>
#pragma warning(pop)
//#include <ogdf/decomposition/BCTree.h>
//#include <ogdf/planarity/PlanarizationGridLayout.h>


BCTree* initVbcAndSingleParts(InputGraph& iG);
void showtime(time_t ticks);
void deletenode(Graph &G, int nodeIndex);
void connectnodes(Graph &G, int nodeIndexS, int nodeIndexT);
bool readLableFile(string& labelFile, LabelFunction &labelFunction, map<string,labelType>& stringLabelToSimpleLabel, vector<string>& simpleLabelToString);

enum class ComputationType { COMP_NONE, COMP_C, COMP_I1, COMP_I2, COMP_OneToX };
enum class SimCoefficient { SC_WallisEtAl=1, SC_BunkeShearer, SC_Asymmetric, SC_NormalizedJohnson, SC_Johnson, SC_SokalSneath, SC_Kulczynski, SC_McConnaughey };

class ReadGraphDB
{
private:
	LabelFunction& m_labelFunction;
	vector<InputGraph*> m_inputGraph;
	map<string,labelType>& stringLabelToSimpleLabel;
	vector<string>& simpleLabelToString;
	size_t readQueueSize;
	const char* fogOutputFilename;
	const char* compGraphsFileName;
	unsigned computationThreads;
	size_t maxMCSInMemory;
	int m_maxCliqueTime;
	wType m_distance_penalty;
	bool m_removeMultiEdgesAndSelfLoops;
	bool m_appendCliqueComputationTimeouts;
	bool verbose;
	ComputationType m_computationType;
	int m_argument1,m_argument2;
	bool m_argument3,m_argument4;
	list<BBP_MCSI*> bbpl;
	ofstream comparedGraphs;
	SimCoefficient m_simCoefficient;

	
	/*
	worker wait until notified
	if queue not empty:
		pop graph, notify main thread if queue was full, work, check for next
	if queue empty
	    if read finished stop
		else wait again

	main thread: read graphs, notify workers cv_work_avail, wait cv_work_processed if queue is full

	
	*/
	static void computeBBPMCS(int cpu=-1);
	static void outputAndDelete(list<BBP_MCSI*>& bbpl, size_t& graph_removal); // output list of BBPMCS into file or console and delete list and BBP_MCSIs

	static wType computeSimilarity(wType wMCS, wType wG, wType wH, SimCoefficient simCoefficient);

public:
	ReadGraphDB(LabelFunction& labelFunction, map<string,labelType>& stringLabelToSimpleLabel, vector<string>& simpleLabelToString,
	int readQueueSize=1, const char* fogOutputFilename=nullptr, const char* compGraphsFileName=nullptr, unsigned computationThreads=0,
	size_t maxMCSInMemory=200, int maxCliqueTime=0, wType distance_penalty=WEIGHT_NOT_COMPATIBLE, bool removeMultiEdgesAndSelfLoops=false, bool appendCliqueComputationTimeouts=false, bool verbose=false,
	ComputationType computationType= ComputationType::COMP_NONE, int argument1=0, int argument2=0, bool argument3=false, bool argument4=false, SimCoefficient simCoefficient= SimCoefficient::SC_WallisEtAl);
	~ReadGraphDB()
	{
		if (comparedGraphs.is_open())
			comparedGraphs.close();
	}
	int readFile(const char* filename);
};
void createFogStarGraphFile(int start, int end, int step);
void createRandomTreeFile(unsigned numTrees, unsigned minsize, unsigned maxsize, unsigned numlabels, unsigned maxdegree=0, bool twoEvenSizedSets=false);
void createRandomOuterplanarFile(unsigned numGraphs, unsigned minsize, unsigned maxsize, unsigned numlabels, bool twoEvenSizedSets, int maxBridgeLength, int maxComponentSize, int bridgeProbabilityPercent, int edgeExtendPercent);
void transform(const char* filenameinput1, const char* filenameinput2, const char* filenameoutput);
#endif
