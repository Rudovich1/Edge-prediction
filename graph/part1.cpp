// #define DEBUG

#include "graph.hpp"
#include "graph-reader.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


int main(int argc, char* argv[])
{  
    std::string str_graph_path;
    std::string str_out_path;
    std::ostream* out = &std::cout;
    #ifdef DEBUG
        argc = 3;
        str_graph_path = "../data/radoslaw_email_email.txt";
    #else
        str_graph_path = argv[1];
    #endif
    
    if (argc < 2 || argc > 3)
    {
        std::cout << "error 1";
        return 1;
    }
    if (argc == 3)
    {
        #ifdef DEBUG
            str_out_path = "log.txt";
        #else
            str_out_path = argv[2];
        #endif

        out = new std::ofstream(str_out_path);
    }

    std::filesystem::path graph_path(str_graph_path);

    Graph::Graph graph(Graph::GraphReader::readGraphData(graph_path));
    graph.setTime(ULLONG_MAX);

    *out << "Number of vertices: " << graph.getNumVertices() << '\n'
        << "Number of edges: " << graph.getNumEdges() << '\n'
        << "Graph density: " << graph.getDensity() << '\n';

    std::vector<std::vector<Graph::Vertex>> components = graph.getWeakComponents();
    size_t max_graph_id = 0;
    for (size_t i = 1; i < components.size(); ++i)
    {
        if (components[i].size() > components[max_graph_id].size())
        {
            max_graph_id = i;
        }
    }

    Graph::Graph max_component = graph.extractSubgraph(components[max_graph_id]);
    
    *out  << "Number components of weak connectivity: " << components.size() << '\n'
        << "Share of vertices in the maximum component: " << (double)max_component.getNumVertices() / graph.getNumVertices() << '\n'
        << "Number of vertices: " << max_component.getNumVertices() << '\n'
        << "Number of edges: " << max_component.getNumEdges() << '\n'
        << "Graph density: " << max_component.getDensity() << '\n' << '\n';

    
    std::vector<Graph::Vertex> random_subgraph = max_component.getRandomSubgraph(500);
    std::vector<Graph::Vertex> snowball_subgraph = max_component.getSnowballSubgraph(500);

    std::vector<size_t> static_data_rand = max_component.getStaticData(random_subgraph);
    std::vector<size_t> static_data_snow = max_component.getStaticData(snowball_subgraph);

    *out << "Radius: ";
    
    if (max_component.getNumVertices() < 4000 && max_component.getNumEdges() * max_component.getNumVertices() < 10'000'000)
    {
        *out << max_component.getRadius() << " / ";
    }
    else
    {
        *out << "__ / ";
    }
    *out << static_data_rand[0] << " / " << static_data_snow[0] << '\n';

    *out << "Diameter: ";
    
    if (max_component.getNumVertices() < 4000 && max_component.getNumEdges() * max_component.getNumVertices() < 10'000'000)
    {
        *out << max_component.getDiameter() << " / ";
    }
    else
    {
        *out << "__ / ";
    }
    *out << static_data_rand[1] << " / " << static_data_snow[1] << '\n';

    *out << "90th percentile: ";
    
    if (max_component.getNumVertices() < 4000 && max_component.getNumEdges() * max_component.getNumVertices() < 10'000'000)
    {
        *out << max_component.getPercentile(90) << " / ";
    }
    else
    {
        *out << "__ / ";
    }
    *out << static_data_rand[2] << " / " << static_data_snow[2] << '\n';
    
    *out << "Average clustering coefficient: ";
    if (max_component.getR2() < 10'000'000)
    {
        *out << max_component.getClusteringCoefficient() << " / ";
    }
    else
    {
        *out << "__ / ";
    }
    *out << max_component.getClusteringCoefficient(random_subgraph) << " / " << max_component.getClusteringCoefficient(snowball_subgraph) << '\n';

    *out << "Pearson correlation coefficient: " << max_component.getPearsonsCoefficient();

    out->flush();
}