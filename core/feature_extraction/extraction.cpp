// #define DEBUG
#define BENCHMARK

#include "graph.hpp"
#include "graph_features.hpp"
#include "csv.hpp"
#include "data_preparation.hpp"
#include "graph_init.hpp"
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <chrono>

class ArgParser
{
public:

    using Args = std::tuple<std::string, std::string, long long, char, double>;

    static Args parse(int argc, char* argv[])
    {
        if (argc != 6)
        {
            throw std::logic_error("Invalid number of arguments");
        }
        if (std::string(argv[4]).size() != 1)
        {
            throw std::logic_error("Invalid 4 argument");
        }
        Args res = {argv[1], argv[2], std::stoll(argv[3]), argv[4][0], std::stod(argv[5])};
        if (std::get<3>(res) != 'f' && std::get<3>(res) != 'l' && std::get<3>(res) != 'p')
        {
            throw std::logic_error("Invalid 4 argument");
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
    long long size;
    char feature_set;
    double init_percent;

    // ArgParser::Args args = {"../../data/munmun_digg_0.3/edges_data.txt", "TMP_RES.txt", -1, 'p', 20};

    try
    {
        ArgParser::Args args = ArgParser::parse(argc, argv);
        ArgParser::openInputFile(in, std::get<0>(args));
        ArgParser::openOutputFile(out, std::get<1>(args));
        size = std::get<2>(args);
        feature_set = std::get<3>(args);
        init_percent = std::get<4>(args);

        Graphs::GraphFeatureEdgesTimestamp graph = init(in, feature_set);
        in.close();

        if (size == -1)
        {
            size = graph.getFullNumEdges() * init_percent / 100.;
        }

        std::vector<std::vector<double>> data = graph.getDataPartWithCalc(0., init_percent, size / 2, size / 4, size / 4);

        std::vector<std::string> columns(data[0].size());
        for (size_t i = 0; i < columns.size() - 4; ++i)
        {
            columns[i] = "fun " + std::to_string(i + 1);
        }
        columns[columns.size() - 4] = "time_calc";
        columns[columns.size() - 3] = "time_init";
        columns[columns.size() - 2] = "time_predict";
        columns[columns.size() - 1] = "predict";

        CSV::toCSV(columns, data, out);

        out.close();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    
}