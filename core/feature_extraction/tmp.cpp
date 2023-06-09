#include <string>
#include <map>
#include <chrono>
#include <iostream>
#include <iomanip>

// #define BENCHMARK

// #ifdef BENCHMARK

struct Benchmark
{
    using TimePoint = std::chrono::steady_clock::time_point;

    static inline std::map<std::string, std::chrono::seconds> times_;
    std::string name_;
    TimePoint init_time_;

    Benchmark(const std::string& name): name_(name), init_time_(std::chrono::steady_clock::now()){}
    
    static long long getTime(const std::string& name)
    {
        if (times_.find(name) == times_.end())
        {
            return -1;
        }
        return times_[name].count();
    }

    static void printRes()
    {
        std::cout << "\nNum functions: " << times_.size() << "\n\n";

        for (auto& i: times_)
        {
            std::cout << i.first << ": " << std::setfill('-') << std::setw(20) << " " << i.second.count() << "s.\n"; 
        }
    }

    static void clear()
    {
        times_.clear();
    }

    ~Benchmark()
    {
        if (times_.find(name_) == times_.end())
        {
            times_[name_] = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - init_time_);
        }
        else
        {
            times_[name_] += std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - init_time_);
        }
    }
};


#define _CLASS_NAME_ typeid(*this).name()
#define _FUNCTION_NAME_ __builtin_FUNCTION()
#define _STR_(c_str) std::string(c_str)

#define _BENCHMARK_ Benchmark bench(_STR_(_));

// #endif


void a()
{
    _BENCHMARK_
    int calc = 100*100;
}

void b()
{
    _BENCHMARK_
    int calc = 100*100;
}


int main()
{
    for (int i = 0; i < 100; ++i)
    {
        a();
        b();
    }
    Benchmark::printRes();
}
