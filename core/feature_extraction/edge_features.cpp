#include "graph.h"
#include "graph_features.h"
#include "graph_init.hpp"
#include "csv.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <exception>

class ArgParser
{
public:

    using Args = std::tuple<std::string, std::string, char, size_t, size_t, size_t, size_t>;

    static Args parse(int argc, char* argv[])
    {
        if (argc != 8)
        {
            throw std::logic_error("Invalid number of argument");
        }
        if (std::string(argv[3]).size() != 1)
        {
            throw std::logic_error("Invalid 3 argument");
        }
        Args res = {argv[1], argv[2], argv[3][0], std::stoull(argv[4]), std::stoull(argv[5]), std::stoull(argv[6]), std::stoull(argv[7])};
        if (std::get<2>(res) != 'f' && std::get<2>(res) != 'l')
        {
            throw std::logic_error("Invalid 3 argument");
        }
        if (std::get<5>(res) >= std::get<6>(res))
        {
            throw std::logic_error("Invalid time");
        }
        return res;
    }

    static void openInputFile(std::ifstream& ifs, const std::string& file_name)
    {
        ifs.open(file_name);
        if (!ifs.is_open())
        {
            throw std::runtime_error((std::string("File \"") + file_name + std::string("\" is not open")).c_str());
        }
    }

    static void openOutputFile(std::ofstream& ofs, const std::string& file_name)
    {
        ofs.open(file_name);
        if (!ofs.is_open())
        {
            throw std::runtime_error((std::string("File \"") + file_name + std::string("\" is not open")).c_str());
        }
    }
};

int main(int argc, char* argv[])
{
    std::ifstream in;
    std::ofstream out;
    char feature_set;
    Graphs::Graph::Edge edge;
    size_t t1;
    size_t t2;

    // in.open("../../data/download.tsv.edit-krwikiquote.txt");
    // out.open("res.txt");
    // feature_set = 'f';
    // edge = {1, 6};
    // t1 = 10;
    // t2 = 20;
    
    try
    {
        ArgParser::Args args = ArgParser::parse(argc, argv);
        ArgParser::openInputFile(in, std::get<0>(args));
        ArgParser::openOutputFile(out, std::get<1>(args));
        feature_set = std::get<2>(args);
        edge = {std::get<3>(args) - 1, std::get<4>(args) - 1};
        t1 = std::get<5>(args);
        t2 = std::get<6>(args);

        Graphs::GraphFeatureEdgesTimestamp graph = init(in, feature_set);
        in.close();
        graph.initTime(t1);
        if (graph.getEdgeTime(edge) <= t1)
        {
            throw std::logic_error((
                std::string("Edge {") + 
                std::to_string(edge.first + 1) + 
                std::string(", ") + 
                std::to_string(edge.second + 1) + 
                std::string("} already exists")
                ).c_str());
        }

        Graphs::GraphFeatureEdges::FeatureVector fv = graph.getData(t1, {edge})[0];
        std::vector<std::string> columns(fv.size() + 2);
        for (size_t i = 0; i < columns.size() - 4; ++i)
        {
            columns[i] = "fun " + std::to_string(i + 1);
        }
        columns[columns.size() - 4] = "time_calc";
        columns[columns.size() - 3] = "time_init";
        columns[columns.size() - 2] = "time_predict";
        columns[columns.size() - 1] = "predict";

        fv.emplace_back(graph.getNormTime(t2, t1));
        if (fv[fv.size() - 2] <= fv[fv.size() - 1] && fv[fv.size() - 2] != -1.L)
        {
            fv.emplace_back(1.L);
        }
        else
        {
            fv.emplace_back(0.L);
        }

        CSV::toCSV(columns, {fv}, out);
        out.close();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    
}