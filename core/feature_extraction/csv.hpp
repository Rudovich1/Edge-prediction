#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

namespace CSV
{
    void toCSV(const std::vector<std::string>& columns, const std::vector<std::vector<double>>& data, std::ofstream& out)
    {
        out << std::fixed << std::setprecision(9);
        if (!columns.empty())
        {
            for (size_t i = 0; i < columns.size() - 1; ++i)
            {
                out << columns[i] << ",";
            }
            out << columns[columns.size() - 1] << std::endl;
        }

        for (size_t i = 0; i < data.size(); ++i)
        {
            for (size_t j = 0; j < data[i].size() - 1; ++j)
            {
                out << data[i][j] << ',';
            }
            out << data[i][data[i].size() - 1] << std::endl;
        }
    }
}