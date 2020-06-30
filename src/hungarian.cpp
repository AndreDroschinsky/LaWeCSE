//
//  hungarian.cpp
//  mco
//
//  Created by Fritz Bökler on 03.04.14.
//
//

#include <global.h>

#include <hungarian.h>

//using std::function;
/*using std::list;
using std::deque;
using std::abs;*/

#pragma warning(push, 0)
#include <ogdf/basic/Graph.h>
#pragma warning(pop)

/*using ogdf::Graph;
using ogdf::edge;
using ogdf::node;
using ogdf::NodeArray;
using ogdf::EdgeArray;*/


namespace mco {

void augment(node,
                    NodeArray<node>& exposed_,
                    NodeArray<node>& mate_,
                    NodeArray<node>& label_);


wType HungarianMethod::
Solve(const Graph& graph,
    const EdgeArray<wType>& egdeweight,
    const vector<node>& agents,
    const vector<node>& non_agents) {

    // EQDIST
    wType infinity_=numeric_limits<wType>::infinity();

    NodeArray<list<node>>   A_(graph);
    NodeArray<node>         mate_(graph);
    NodeArray<node>         exposed_(graph);
    NodeArray<node>         neighbour_(graph);
    NodeArray<node>         label_(graph);

    NodeArray<wType>        dual_variables_(graph, 0);
    NodeArray<wType>        slack_(graph, infinity_);
    NodeArray<unsigned>     count_(graph, 0);
    deque<node>             queue_;

    //list<node> non_agents;

    edge e;

    for(node n : agents)
        mate_[n]=n;

    for(node n : non_agents){
    //forall_nodes(n, graph) {
        mate_[n] = n;
        edge minimum = n->firstAdj()->theEdge();
        //forall_adj_edges(e, n)
        for(adjEntry adj : n->adjEntries)
        {
            e=adj->theEdge();
            if (egdeweight[e]<egdeweight[minimum])
            {
                minimum = e;
            }
        }
        dual_variables_[n] = egdeweight[minimum];
    }

    bool endstage;
    for(unsigned int s = 0; s < agents.size(); ++s) {
        endstage = false;

        // Initialize agents
        for(auto agent : agents) {
            A_[agent].clear();
            exposed_[agent] = agent;
            label_[agent] = agent;
            count_[agent] = 0;
        }

        // Initialize jobs
        for(auto job : non_agents) {
            slack_[job] = infinity_;
            neighbour_[job] = job;
        }

        // Initialize equality subgraph
        forall_edges(e, graph)
        if(abs(egdeweight[e] - (dual_variables_[e->source()] + dual_variables_[e->target()])) <= EQDIST) {
            auto job = e->target();
            auto agent = e->source();
            if(mate_[job] == job)
            exposed_[agent] = job;
            else if(agent != mate_[job])
            A_[agent].push_back(mate_[job]);
        }

        queue_.clear();
        for(auto agent : agents)
        if(mate_[agent] == agent) {
            if(exposed_[agent] != agent) {
                endstage = true;
                augment(agent, exposed_, mate_, label_);
                break;
            }
            queue_.push_back(agent);
            label_[agent] = agent;
            for(adjEntry adj : agent->adjEntries)
            {
                e=adj->theEdge();
            //forall_adj_edges(e, agent) {
                if(egdeweight[e]-dual_variables_[e->source()]-dual_variables_[e->target()] >= -EQDIST &&
                 slack_[e->target()] - (egdeweight[e]-dual_variables_[e->source()]-dual_variables_[e->target()]) >= EQDIST ) {

                       slack_[e->target()] = egdeweight[e] - dual_variables_[e->source()] - dual_variables_[e->target()];
                       if(neighbour_[e->target()] != e->target())
                       count_[neighbour_[e->target()]] -= 1;
                       count_[agent] += 1;
                       neighbour_[e->target()] = agent;
                   }
            }
        }

        if(endstage) {
           continue;
        }

        while(true) {

            while(!queue_.empty()) {
                auto agent = queue_.front();
                queue_.pop_front();
                for(auto next_agent : A_[agent])
                if(label_[next_agent] == next_agent) {
                    label_[next_agent] = agent;

                    if(exposed_[next_agent] != next_agent) {
                        augment(next_agent, exposed_, mate_, label_);
                        endstage = true;
                        break;
                    }

                    queue_.push_back(next_agent);

                    for(adjEntry adj : next_agent->adjEntries)
                    {
                        e=adj->theEdge();
                    //forall_adj_edges(e, next_agent) {
                        if(egdeweight[e]-dual_variables_[e->source()]-dual_variables_[e->target()] >= -EQDIST &&
                         slack_[e->target()] - (egdeweight[e]-dual_variables_[e->source()]-dual_variables_[e->target()]) >= EQDIST ) {

                               slack_[e->target()] = egdeweight[e] - dual_variables_[e->source()] - dual_variables_[e->target()];
                               if(neighbour_[e->target()] != e->target())
                               count_[neighbour_[e->target()]] -= 1;
                               count_[next_agent] += 1;
                               neighbour_[e->target()] = next_agent;
                           }
                    }
                }

                if(endstage) {
                    break;
                }
            }


            if(endstage)
            break;

            wType minimum_slack = infinity_;

            for(auto job : non_agents) {
                if(slack_[job] > EQDIST && slack_[job]-minimum_slack < -EQDIST) {
                    minimum_slack = slack_[job];
                }
            }

            wType theta = minimum_slack;
            theta *= 0.5;

            for(auto agent : agents) {
                if(label_[agent] != agent || count_[agent] > 0) {
                    dual_variables_[agent] += theta;
                } else {
                    dual_variables_[agent] -= theta;
                }
            }

            for(auto job : non_agents) {
                if(abs(slack_[job]) <= EQDIST) {
                    dual_variables_[job] -= theta;
                } else {
                    dual_variables_[job] += theta;
                }
            }

            theta = minimum_slack;

            for(auto job : non_agents)
            if(slack_[job] > EQDIST) {
                slack_[job] -= theta;
                if(abs(slack_[job]) <= EQDIST) {
                    if(mate_[job] == job) {
                        exposed_[neighbour_[job]] = job;
                        augment(neighbour_[job], exposed_, mate_, label_);
                        endstage = true;
                        break;
                    } else {
                        queue_.push_back(neighbour_[job]);
                        A_[neighbour_[job]].push_back(mate_[job]);
                    }
                }
            }

            if(endstage) {
                break;
            }
        }
    }

    value_ = 0;

    for(auto n : graph.nodes)
        value_ += dual_variables_[n];




#ifndef LAWECSU_NDEBUG
    for(node n : agents)
        std::cout << "(" << n << ", " << mate_[n] << "), ";
    //for(auto n : graph.nodes) {
        //if(n->index() < (signed) agents.size()) {
          //  std::cout << "(" << n << ", " << mate_[n] << "), ";
        //}
    //}
    std::cout << "Duals ";
    for(auto n : graph.nodes) {
        cout << "("<<n<<":"<<dual_variables_[n]<<"), ";
    }
    std::cout << "-- ";

#endif


#ifndef LAWECSU_NDEBUG
//    std::cout << value_ << std::endl;
#endif

    return value_;
}




wType HungarianMethod::
Solve_Max(const Graph& graph,
    const EdgeArray<wType>& egdeweight,
    const vector<node>& agents,
    const vector<node>& non_agents,
    NodeArray<node>& mate_,
    NodeArray<wType>& dual_variables_,
    unsigned iterations) {

    // EQDIST
    wType infinity_=numeric_limits<wType>::infinity();

    NodeArray<list<node>>   A_(graph);
    NodeArray<node>         exposed_(graph);
    NodeArray<node>         neighbour_(graph);
    NodeArray<node>         label_(graph);

    NodeArray<wType>        slack_(graph, infinity_);
    NodeArray<unsigned>     count_(graph, 0);
    deque<node>             queue_;

    edge e;

    bool endstage;
    for(unsigned s = 0; s < iterations; ++s) {
        endstage = false;

        // Initialize agents
        for(node agent : agents) {
            A_[agent].clear();
            exposed_[agent] = agent;
            label_[agent] = agent;
            count_[agent] = 0;
        }

        // Initialize jobs
        for(node job : non_agents) {
            slack_[job] = infinity_;
            neighbour_[job] = job;
        }

        // Initialize equality subgraph
        forall_edges(e, graph)
        if(abs(-egdeweight[e] - (dual_variables_[e->source()] + dual_variables_[e->target()])) <= EQDIST) {
            node job = e->target();
            node agent = e->source();
            if(mate_[job] == job)
            exposed_[agent] = job;
            else if(agent != mate_[job])
            A_[agent].push_back(mate_[job]);
        }

        queue_.clear();
        for(node agent : agents)
        if(mate_[agent] == agent) {
            if(exposed_[agent] != agent) {
                endstage = true;
                augment(agent, exposed_, mate_, label_);
                break;
            }
            queue_.push_back(agent);
            label_[agent] = agent;
            for(adjEntry adj : agent->adjEntries)
            {
                e=adj->theEdge();
            //forall_adj_edges(e, agent) {
                wType slack=-egdeweight[e]-dual_variables_[e->source()]-dual_variables_[e->target()];
                if(slack >= -EQDIST &&
                 slack_[e->target()] - slack >= EQDIST ) {

                       slack_[e->target()] = -egdeweight[e] - dual_variables_[e->source()] - dual_variables_[e->target()];
                       if(neighbour_[e->target()] != e->target())
                           count_[neighbour_[e->target()]] -= 1;
                       count_[agent] += 1;
                       neighbour_[e->target()] = agent;
                   }
            }
        }

        if(endstage) {
           continue;
        }

        while(true) {

            while(!queue_.empty()) {
                node agent = queue_.front();
                queue_.pop_front();
                for(node next_agent : A_[agent])
                if(label_[next_agent] == next_agent) {
                    label_[next_agent] = agent;

                    if(exposed_[next_agent] != next_agent) {
                        augment(next_agent, exposed_, mate_, label_);
                        endstage = true;
                        break;
                    }

                    queue_.push_back(next_agent);

                    for(adjEntry adj : next_agent->adjEntries)
                    {
                        e=adj->theEdge();
                    //forall_adj_edges(e, next_agent) {
                        wType slack=-egdeweight[e]-dual_variables_[e->source()]-dual_variables_[e->target()];
                        if(slack >= -EQDIST &&
                         slack_[e->target()] - slack >= EQDIST ) {

                               slack_[e->target()] = -egdeweight[e] - dual_variables_[e->source()] - dual_variables_[e->target()];
                               if(neighbour_[e->target()] != e->target())
                                   count_[neighbour_[e->target()]] -= 1;
                               count_[next_agent] += 1;
                               neighbour_[e->target()] = next_agent;
                           }
                    }
                }

                if(endstage) {
                    break;
                }
            }

            if(endstage)
            break;

            wType minimum_slack = infinity_;

            for(auto job : non_agents) {
                if(slack_[job] > EQDIST && slack_[job]-minimum_slack < -EQDIST) {
                    minimum_slack = slack_[job];
                }
            }

            wType theta = minimum_slack;
            theta *= 0.5;

            for(node agent : agents) {
                if(label_[agent] != agent || count_[agent] > 0) {
                    dual_variables_[agent] += theta;
                } else {
                    dual_variables_[agent] -= theta;
                }
            }

            for(node job : non_agents) {
                if(abs(slack_[job]) <= EQDIST) {
                    dual_variables_[job] -= theta;
                } else {
                    dual_variables_[job] += theta;
                }
            }

            theta = minimum_slack;

            for(node job : non_agents)
            if(slack_[job] > EQDIST) {
                slack_[job] -= theta;
                if(abs(slack_[job]) <= EQDIST) {
                    if(mate_[job] == job) {
                        exposed_[neighbour_[job]] = job;
                        augment(neighbour_[job], exposed_, mate_, label_);
                        endstage = true;
                        break;
                    } else {
                        queue_.push_back(neighbour_[job]);
                        A_[neighbour_[job]].push_back(mate_[job]);
                    }
                }
            }

            if(endstage) {
                break;
            }
        }
    }

    value_ = 0;

    for(node n : agents)
        value_ += dual_variables_[n];
    for(node n : non_agents)
        value_ += dual_variables_[n];



#ifndef LAWECSU_NDEBUG


 /*   std::cout << "Duals ";
    for(auto n : graph.nodes) {
        cout << "("<<n<<":"<<dual_variables_[n]<<"), ";
    }
    std::cout << "-- ";*/

#endif


#ifndef LAWECSU_NDEBUG
//    std::cout << value_ << std::endl;
#endif
    if (value_<0)
        value_=-value_;
    value_ /= 2.0;
    return value_;
}








void augment(node v,
                    NodeArray<node>& exposed_,
                    NodeArray<node>& mate_,
                    NodeArray<node>& label_) {

    while(label_[v] != v) {
        exposed_[label_[v]] = mate_[v];
        mate_[v] = exposed_[v];
        mate_[exposed_[v]] = v;
        v = label_[v];
    }
    mate_[v] = exposed_[v];
    mate_[exposed_[v]] = v;
}

}
