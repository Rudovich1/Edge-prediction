#pragma once
#include <vector>
#include <random>

namespace DataPreproc
{
    using DataSet = std::vector<std::vector<double>>;
    void preproc(DataSet& data, double time)
    {
        for (auto& i: data)
        {
            i.emplace_back(time);
            if (i[i.size() - 2] <= time && i[i.size() - 2] != -1.)
            {
                i.emplace_back(1.);
            }
            else
            {
                i.emplace_back(0.);
            }
        }
    }

    DataSet norm(const DataSet& data)
    {
        DataSet neg, pos;
        for (auto& i: data)
        {
            if (i[i.size() - 1] == 0.L)
            {
                neg.emplace_back(i);
            }
            else
            {
                pos.emplace_back(i);
            }
        }
        size_t res_size = std::min(neg.size(), pos.size());

        std::random_device rd;
        std::mt19937 gen(rd());

        std::shuffle(neg.begin(), neg.end(), gen);
        std::shuffle(pos.begin(), pos.end(), gen);

        neg.resize(res_size);
        pos.resize(res_size);

        for (auto& i: pos)
        {
            neg.emplace_back(i);
        }

        std::shuffle(neg.begin(), neg.end(), gen);

        return neg;
    }
}