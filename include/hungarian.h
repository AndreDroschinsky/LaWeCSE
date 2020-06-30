//
//  hungarian.h
//  mco
//
//  Created by Fritz Bökler on 03.04.14.
//
//

#ifndef __mco__lex__hungarian__
#define __mco__lex__hungarian__

//#include <functional>
#include <list>
#include <deque>
#include <vector>
#include <global.h>

#pragma warning(push, 0)
#include <ogdf/basic/Graph.h>
#pragma warning(pop)

//#include "point.h"
using namespace ogdf;
using namespace std;

namespace mco {
class HungarianMethod {
public:
    HungarianMethod() : value_(0) { }

    wType Solve(const Graph& graph,
                const EdgeArray<wType>& egdeweight,
                const vector<node>& agents,
                const vector<node>& non_agents
                );

    wType Solve_Max(const Graph& graph,
                const EdgeArray<wType>& egdeweight,
                const vector<node>& agents,
                const vector<node>& non_agents,
                NodeArray<node>& mate_,
                NodeArray<wType>& dual_variables_,
                unsigned iterations
                );
    
    wType value() {
            return value_;//_/2.0;
    }
    
private:
    wType value_;
    
};
}
#endif /* defined(__mco__hungarian__) */
