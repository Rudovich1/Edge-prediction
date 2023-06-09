#pragma once

#include <cmath>
#include <vector>
#include "graph.hpp"

namespace Graphs
{
    namespace GraphFeatures
    {
        double distance_V1_V2(Graph& graph, Graph::Edge edge)
        {
            double dist = graph.getDistance(edge.first, edge.second);
            if (dist == -1) 
            {
                return 1.;
            } 
            return dist / (double)(graph.getNumberVertices());
        }

        double numNeighbors(Graph& graph, Graph::Edge edge)
        {
            return (double)graph.getNumdAdjacentVertices(edge) / (graph.getNumberVertices() - 2);
        }

        double V1_Degree(Graph& graph, Graph::Edge edge)
        {
            return (double)graph.getNumberNeighbors(edge.first) / (graph.getNumberVertices() - 2);
        }

        double V2_Degree(Graph& graph, Graph::Edge edge)
        {
            return (double)graph.getNumberNeighbors(edge.second) / (graph.getNumberVertices() - 2);
        }

        double numSharedNeighbors(Graph& graph, Graph::Edge edge)
        {
            return (double)graph.getNumSharedAdjacentVertices(edge) / (graph.getNumberVertices() - 2);
        }

        double numNonAdjacentEdges(Graph& graph, Graph::Edge edge)
        {
            double max_number_edges = (graph.getNumberVertices()) * (graph.getNumberVertices() - 1) / 2.;
            double vC = graph.getNumdAdjacentVertices(edge);
            double num = graph.getNumberEdges() - vC;
            return num / (max_number_edges - vC);
        }

        double numXORNeighbors(Graph& graph, Graph::Edge edge)
        {
            double comb = graph.getNumdAdjacentVertices(edge);
            double inter = graph.getNumSharedAdjacentVertices(edge);
            double graph_size_tmp = graph.getNumberVertices();
            return (comb - inter) / (graph_size_tmp - 2);
        }

        double averageDistanceFrom_V1_ToNeighbors_V2(Graph& graph, Graph::Edge edge)
        {
            double mean_dist = 0;
            long long num = 0;

            for (auto& i: graph.getDistances(edge.first))
            {
                if (graph.inGraph({i.first, edge.second}))
                {
                    mean_dist += i.second;
                    ++num;
                }
            }

            if (num == 0)
            {
                return 1.;
            }

            mean_dist /= num;
            double graph_size_tmp = graph.getNumberVertices(); 
            return mean_dist / (graph_size_tmp);
        }

        double averageDistanceFrom_V2_ToNeighbors_V1(Graph& graph, Graph::Edge edge)
        {
            double mean_dist = 0;
            long long num = 0;

            for (auto& i: graph.getDistances(edge.second))
            {
                if (graph.inGraph({i.first, edge.first}))
                {
                    mean_dist += i.second;
                    ++num;
                }
            }

            if (num == 0)
            {
                return 1.;
            }

            mean_dist /= num;
            double graph_size_tmp = graph.getNumberVertices(); 
            return mean_dist / (graph_size_tmp);
        }

        double R1(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getR1();
            return num / (graph_size * (graph_size - 1));
        }

        double R2(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getR2();
            return num / (graph_size * std::pow(graph_size - 1, 2));
        }

        double R3(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getR3();
            return num / (graph_size * std::pow(graph_size - 1, 3));
        }

        double Re(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getRe();
            return num / (graph_size * std::pow(graph_size - 1, 3));
        }

        double pearsonCorrelationCoefficient(Graph& graph, Graph::Edge edge)
        {
            return (graph.getPearsonCorrelationCoefficient() + 1) / 2.;
        }

        double clusteringCoefficient(Graph& graph, Graph::Edge edge)
        {
            return graph.getClusteringCoefficient();
        }

        double graphDensity(Graph& graph, Graph::Edge edge)
        {
            return graph.getGraphDensity();
        }

        double eccentricityV1(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getEccentricity(edge.first);
            return num / (graph_size - 1);
        }

        double eccentricityV2(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getEccentricity(edge.second);
            return num / (graph_size - 1);
        }

        double radius(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getRadius();
            return num / (graph_size - 1);
        }

        double diameter(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getDiameter();
            return num / (graph_size - 1);
        }

        double percentile90(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getPercentile(90);
            return num / (graph_size - 1);
        }

        long double percentile50(Graph& graph, Graph::Edge edge)
        {
            double graph_size = graph.getNumberVertices();
            double num = graph.getPercentile(50);
            return num / (graph_size - 1);
        }

        double localSubgraph_R1(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return R1(local_subgraph, new_edge);
        }

        double localSubgraph_R2(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return R2(local_subgraph, new_edge);
        }

        double localSubgraph_R3(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return R3(local_subgraph, new_edge);
        }

        double localSubgraph_Re(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return Re(local_subgraph, new_edge);
        }

        double localSubgraph_PearsonCorrelationCoefficient(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return pearsonCorrelationCoefficient(local_subgraph, new_edge);
        }

        double localSubgraph_ClusteringCoefficient(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return clusteringCoefficient(local_subgraph, new_edge);
        }

        double localSubgraph_GraphDensity(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return graphDensity(local_subgraph, new_edge);
        }

        double localSubgraph_EccentricityV1(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return eccentricityV1(local_subgraph, new_edge);
        }

        double localSubgraph_EccentricityV2(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return eccentricityV2(local_subgraph, new_edge);
        }

        double localSubgraph_Radius(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return radius(local_subgraph, new_edge);
        }

        double localSubgraph_Diameter(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return diameter(local_subgraph, new_edge);
        }

        double localSubgraph_90Percentile(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return percentile90(local_subgraph, new_edge);
        }

        double localSubgraph_50Percentile(Graph& graph, Graph::Edge edge)
        {
            auto[local_subgraph, new_edge] = graph.getLocalSubgraph(edge);
            return percentile50(local_subgraph, new_edge);
        }

        double productNeighbors(Graph& graph, Graph::Edge edge)
        {
            double tmp1 = (double)graph.getNumberNeighbors(edge.first) / (graph.getNumberVertices() - 1);
            double tmp2 = (double)graph.getNumberNeighbors(edge.second) / (graph.getNumberVertices() - 1);
            return tmp1 * tmp2; 
        }

        double jacquardCoefficient(Graph& graph, Graph::Edge edge)
        {
            return graph.getJacquardCoefficient(edge);
        }

        double adarCoefficient(Graph& graph, Graph::Edge edge)
        {
            return graph.getAdarCoefficient(edge) * std::log(graph.getNumberVertices()) / (graph.getNumberVertices() - 2);
        }

        std::vector<std::function<double(Graph&, Graph::Edge)>> getAllFeatures()
        {
            std::vector<std::function<double(Graph&, Graph::Edge)>> features = 
            {
                distance_V1_V2,
                numNeighbors,
                V1_Degree,
                V2_Degree,
                numSharedNeighbors,
                numNonAdjacentEdges,
                numXORNeighbors,
                averageDistanceFrom_V1_ToNeighbors_V2,
                averageDistanceFrom_V2_ToNeighbors_V1,
                R1,
                R2,
                R3,
                Re,
                pearsonCorrelationCoefficient,
                clusteringCoefficient,
                graphDensity,
                eccentricityV1,
                eccentricityV2,
                radius,
                diameter,
                percentile90,
                percentile50,
                productNeighbors,
                jacquardCoefficient,
                adarCoefficient,
                localSubgraph_R1,
                localSubgraph_R2,
                localSubgraph_R3,
                localSubgraph_Re,
                localSubgraph_PearsonCorrelationCoefficient,
                localSubgraph_ClusteringCoefficient,
                localSubgraph_GraphDensity,
                localSubgraph_EccentricityV1,
                localSubgraph_EccentricityV2,
                localSubgraph_Radius,
                localSubgraph_Diameter,
                localSubgraph_90Percentile,
                localSubgraph_50Percentile
            };

            return features;
        }

        std::vector<std::function<double(Graph&, Graph::Edge)>> getSimpleFeatures()
        {
            std::vector<std::function<double(Graph&, Graph::Edge)>> features = 
            {
                distance_V1_V2,
                numNeighbors,
                V1_Degree,
                V2_Degree,
                numSharedNeighbors,
                numNonAdjacentEdges,
                numXORNeighbors,
                averageDistanceFrom_V1_ToNeighbors_V2,
                averageDistanceFrom_V2_ToNeighbors_V1,
                R1,
                R2,
                R3,
                Re,
                pearsonCorrelationCoefficient,
                clusteringCoefficient,
                graphDensity,
                productNeighbors,
                jacquardCoefficient,
                adarCoefficient
            };
            return features;
        }
    
        std::vector<std::function<double(Graph&, Graph::Edge)>> getPredictFeatures()
        {
            std::vector<std::function<double(Graph&, Graph::Edge)>> features = 
            {
                distance_V1_V2,
                numNeighbors,
                V1_Degree,
                V2_Degree,
                numSharedNeighbors,
                numNonAdjacentEdges,
                numXORNeighbors,
                averageDistanceFrom_V1_ToNeighbors_V2,
                averageDistanceFrom_V2_ToNeighbors_V1,
                R1,
                R2,
                R3,
                Re,
                pearsonCorrelationCoefficient,
                clusteringCoefficient,
                graphDensity,
                productNeighbors,
                jacquardCoefficient,
                adarCoefficient,
                localSubgraph_R1,
                localSubgraph_R2,
                localSubgraph_R3,
                localSubgraph_Re,
                localSubgraph_PearsonCorrelationCoefficient,
                localSubgraph_ClusteringCoefficient,
                localSubgraph_GraphDensity
            };

            return features;
        }
    }
}