#define BENCHMARK

#include "graph.hpp"
#include "graph_features.hpp"
#include "graph_init.hpp"
#include "data_preparation.hpp"
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

    using Args = std::tuple<std::string, std::string, char, long long, double, double>;

    static Args parse(int argc, char* argv[])
    {
        if (argc != 7)
        {
            throw std::logic_error("Invalid number of argument");
        }
        if (std::string(argv[3]).size() != 1)
        {
            throw std::logic_error("Invalid 3 argument");
        }
        Args res = {argv[1], argv[2], argv[3][0], std::stoll(argv[4]), std::stod(argv[5]), std::stod(argv[6])};
        if (std::get<2>(res) != 'f' && std::get<2>(res) != 'l' && std::get<2>(res) != 'p')
        {
            throw std::logic_error("Invalid 3 argument");
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
    long long size;
    double init_percent;
    double predict_percent;

    // ArgParser::Args args = {"../../data/munmun_digg_0.3/edges_data.txt", "tmp_res_123.txt", 'p', -1, 5, 10};
    
    try
    {
        ArgParser::Args args = ArgParser::parse(argc, argv);
        ArgParser::openInputFile(in, std::get<0>(args));
        ArgParser::openOutputFile(out, std::get<1>(args));
        feature_set = std::get<2>(args);
        size = std::get<3>(args);
        init_percent = std::get<4>(args);
        predict_percent = std::get<5>(args);

        Graphs::GraphFeatureEdgesTimestamp graph = init(in, feature_set);
        in.close();

        if (size == -1)
        {
            size = graph.getFullNumEdges() * (predict_percent - init_percent) / 100.;
        }

        std::vector<std::vector<double>> data = graph.getDataWithCalc(init_percent, predict_percent, size / 2, size / 4, size / 4);

        #ifdef BENCHMARK
            Benchmark::printRes();
        #endif

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