#pragma once
#include <istream>
#include "graph.hpp"
#include "graph_features.hpp"

Graphs::GraphFeatureEdgesTimestamp init(std::istream& in, char feature_set)
{
    long long n, m;
    in >> n >> m;
    Graphs::Graph::ListEdges LE(m);
    std::vector<double> TS(m);

    for (long long i = 0; i < m; ++i)
    {
        in >> LE[i].first >> LE[i].second;
        --LE[i].first;
        --LE[i].second;
        in >> TS[i];
    }

    std::vector<std::function<double(Graphs::Graph &, Graphs::Graph::Edge)>> features;

    if (feature_set == 'f')
    {
        features = Graphs::GraphFeatures::getAllFeatures();
    }
    else if (feature_set == 'l')
    {
        features = Graphs::GraphFeatures::getSimpleFeatures();
    }
    else if (feature_set == 'p')
    {
        features = Graphs::GraphFeatures::getPredictFeatures();
    }

    return Graphs::GraphFeatureEdgesTimestamp(n, LE, TS, features);
}