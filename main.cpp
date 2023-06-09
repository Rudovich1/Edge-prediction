#include "graph.h"
#include "graph_features.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <set>

using namespace Graphs;

Graph initGraph(std::ifstream& in)
{
    size_t listSize, graphSize;
    in >> listSize >> graphSize;
    std::set<Graph::Edge> st;
    size_t x, y;
    for (size_t i = 0; i < listSize; ++i)
    {
        
        in >> x >> y;
        --x;
        --y;
        if (x != y)
        {
            st.insert({std::min(x, y), std::max(x, y)});
        }
    }
    Graph::ListEdges edges;
    for (auto& i: st)
    {
        edges.emplace_back(i);
    }
    return Graph(edges, graphSize);
}

void graphInfo(Graph& graph, std::ostream& out, bool full = false)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    graph.init();
    out << "Число вершин: " << graph.getNumberVertices() << '\n';
    out << "Число рёбер: " << graph.getNumberEdges() << '\n';
    out << "Плотность: " << graph.getGraphDensity() << '\n';
    out << "Число компонент слабой связности: " << graph.getNumWeakComponents() << '\n';
    out << "Доля вершин в максимальной по мощности компоненте слабой связности: " << (long double)graph.getMaxComponentSize() / ((long double) graph.getNumberVertices()) << '\n';
    
    Graph maxComponent(graph.getMaxComponent());
    maxComponent.init();

    out << "R1 наибольшей компоненты слабой связности: " << maxComponent.getR1() << '\n';
    out << "R2 наибольшей компоненты слабой связности: " << maxComponent.getR2() << '\n';
    out << "R3 наибольшей компоненты слабой связности: " << maxComponent.getR3() << '\n';
    out << "Re наибольшей компоненты слабой связности: " << maxComponent.getRe() << '\n';

    // out << "L наибольшей компоненты слабой связности: [ ";
    // std::vector<size_t> tmp = maxComponent.getL();
    // for (size_t i = 0; i < tmp.size() - 1; ++i)
    // {
    //     out << tmp[i] << ", ";
    // }

    // out << tmp[tmp.size() - 1] << " ]" << '\n';

    out << "Радиус наибольшей компоненты слабой связности (реальное \\ оценка 1 \\ оценка 2): " ;
    if (full)
    {
        out << maxComponent.getRadius() << "\\";
    }
    else
    {
        out << "X\\";
    }
    out << maxComponent.getRadius(maxComponent.getRandomSubgraph(std::min(maxComponent.getNumberVertices(), 1000ULL))) << "\\"
        << maxComponent.getRadius(maxComponent.getSnowballSampleSubgraph(
            std::min(maxComponent.getNumberVertices(), 1000ULL), std::uniform_int_distribution<size_t>(0, maxComponent.getNumberVertices() - 1)(gen))) << '\n';

    out << "Диаметр наибольшей компоненты слабой связности (реальное \\ оценка 1 \\ оценка 2): "; 
    if (full)
    {    
        out << maxComponent.getDiameter() << "\\";
    }
    else
    {
        out << "X\\";
    }
    out << maxComponent.getDiameter(maxComponent.getRandomSubgraph(std::min(maxComponent.getNumberVertices(), 1000ULL))) << "\\"
        << maxComponent.getDiameter(maxComponent.getSnowballSampleSubgraph(
            std::min(maxComponent.getNumberVertices(), 1000ULL), std::uniform_int_distribution<size_t>(0, maxComponent.getNumberVertices() - 1)(gen))) << '\n';

    out << "50 процентиль расстояния между вершинами наибольшей компоненты слабой связности (реальное \\ оценка 1 \\ оценка 2): " ;
    if (full)
    {
        out << maxComponent.getPercentile(50) << "\\";
    }
    else
    {
        out << "X\\";
    }
    out << maxComponent.getPercentile(50, maxComponent.getRandomSubgraph(std::min(maxComponent.getNumberVertices(), 1000ULL))) << "\\"
        << maxComponent.getPercentile(50, maxComponent.getSnowballSampleSubgraph(
           std::min(maxComponent.getNumberVertices(), 1000ULL), std::uniform_int_distribution<size_t>(0, maxComponent.getNumberVertices() - 1)(gen))) << '\n';
        
    out << "90 процентиль расстояния между вершинами наибольшей компоненты слабой связности (реальное \\ оценка 1 \\ оценка 2): " ;
    if (full)
    {
        out << maxComponent.getPercentile(90) << "\\";
    }
    else
    {
        out << "X\\";
    }
    out << maxComponent.getPercentile(90, maxComponent.getRandomSubgraph(std::min(maxComponent.getNumberVertices(), 1000ULL))) << "\\"
        << maxComponent.getPercentile(90, maxComponent.getSnowballSampleSubgraph(
            std::min(maxComponent.getNumberVertices(), 1000ULL), std::uniform_int_distribution<size_t>(0, maxComponent.getNumberVertices() - 1)(gen))) << '\n';

    out << "Средний кластерный коэффициент сети наибольшей компоненты слабой связности: " << maxComponent.getClusteringCoefficient() << '\n';
    out << "Коэффициент корреляции Пирсона: " << graph.getPearsonCorrelationCoefficient() << '\n' << '\n';
}

int main(int argc, char* argv[])
{
    bool full = false;
    int i = 1;
    if (argc <= i)
    {
        std::cout << "Err 1";
        return 1;
    }
    if (argv[i][0] == '-' && argv[i][1] == 'f')
    {
        full = true;
        ++i;
    }
    if (argc <= i)
    {
        std::cout << "Err 1";
        return 1;
    }
    std::ifstream in(argv[i]);
    if (!in.is_open())
    {
        std::cout << "Err 2";
        return 2;
    }
    ++i;
    if (argc <= i)
    {
        std::cout << "Err 1";
        return 1;
    }
    std::ofstream out(argv[i]);
    if (!out.is_open())
    {
        std::cout << "Err 3";
        return 3;
    }
    Graph graph(initGraph(in));
    graphInfo(graph, out, full);
}