#pragma once

#include <vector>
#include <functional>
#include <queue>
#include <set>
#include <random>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <string>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <queue>

// #define BENCHMARK

#ifdef BENCHMARK

    struct Benchmark
    {
        #define _BENCHMARK std::chrono

        using TimePoint = _BENCHMARK ::steady_clock::time_point;

        struct BenchInfo
        {
            _BENCHMARK ::milliseconds time;
            long long num_calls = 0;
        };

        static inline std::unordered_map<std::string, BenchInfo> times_;
        std::string name_;
        TimePoint init_time_;

        Benchmark(const std::string& name): name_(name), init_time_(std::chrono::steady_clock::now()){}
        
        // static long long getTime(const std::string& name)
        // {
        //     if (times_.find(name) == times_.end())
        //     {
        //         return -1;
        //     }
        //     return times_[name].count();
        // }

        static void printRes()
        {
            std::cout << "\nNum functions: " << times_.size() << "\n\n";

            std::vector<std::pair<std::string, BenchInfo>> res;

            for (auto& i: times_)
            {
                res.emplace_back(i.first, i.second);
            }

            auto comp = [](std::pair<std::string, BenchInfo> a, std::pair<std::string, BenchInfo> b)
            {
                if (a.second.time == b.second.time)
                {
                    if (a.second.num_calls == b.second.num_calls)
                    {
                        return a.first < b.first;
                    }
                    return a.second.num_calls < b.second.num_calls;
                }
                return a.second.time < b.second.time;
            };

            std::sort(res.begin(), res.end(), comp);

            for (auto& i: res)
            {
                std::cout << std::setw(40) << std::setfill('-') << std::setiosflags(std::ios::left) << i.first << " " 
                    << "Time: " << std::setw(10) << std::setiosflags(std::ios::left) << i.second.time.count() << "ms. | "
                    << "Num calls: " << i.second.num_calls << '\n'; 
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
                times_[name_] = {_BENCHMARK ::duration_cast<_BENCHMARK ::milliseconds>(_BENCHMARK ::steady_clock::now() - init_time_), 1};
            }
            else
            {
                times_[name_].time += _BENCHMARK ::duration_cast<_BENCHMARK ::milliseconds>(_BENCHMARK ::steady_clock::now() - init_time_);
                ++times_[name_].num_calls;
            }
        }
    };

    // #define _CLASS_NAME_ typeid(*this).name()
    #define _FUNCTION_NAME_ __builtin_FUNCTION()
    // #define _STR_(c_str) std::string(c_str)

    #define _BENCHMARK_ Benchmark bench(_FUNCTION_NAME_);

#endif

namespace Graphs
{

    class Graph
    {
    public:
        using Vertex = long long;
        using Edge = std::pair<Vertex, Vertex>;
        using AdjacencyList = std::vector<std::vector<Vertex>>;
        using AdjacencyMatrix = std::vector<std::vector<bool>>;
        using ListEdges = std::vector<Edge>;
        using SetEdges = std::set<Edge>;

        template<typename KeyType, typename ValueType>
        struct BufInfo
        {
            BufInfo(long long buf_size = 100): max_buf_size_(buf_size){}

            void insert(const KeyType& key, const ValueType& value)
            {
                if (cur_buf_size_ == max_buf_size_)
                {
                    data_.erase(buf_queue_.front());
                    buf_queue_.pop();
                    --cur_buf_size_;
                }
                buf_queue_.push(key);
                data_[key] = value;
                ++cur_buf_size_;
            }

            ValueType* get(const KeyType& key)
            {
                auto iter = data_.find(key);
                if (iter == data_.end())
                {
                    return nullptr;
                }
                return &iter->second;
            }

            void clear()
            {
                cur_buf_size_ = 0;
                data_.clear();
                while(!buf_queue_.empty())
                {
                    buf_queue_.pop();
                }
            }
        
        protected:

            long long max_buf_size_;
            long long cur_buf_size_ = 0;
            std::map<KeyType, ValueType> data_;
            std::queue<KeyType> buf_queue_;
        };

        Graph(const AdjacencyList& adj): adj_(adj), graph_size_(adj.size()) {}

        Graph(const AdjacencyMatrix& adj): adj_(adj.size()), graph_size_(adj.size())
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            for (Vertex i = 0; i < adj.size(); ++i)
            {
                for (Vertex j = 0; j < adj[i].size(); ++j)
                {
                    if (adj[i][j]) 
                    {
                        adj_[i].emplace_back(j);
                        adj_[j].emplace_back(i);
                    }
                }
            }
        }

        Graph(const ListEdges& LE, size_t graph_size): graph_size_(graph_size), adj_(graph_size)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            for (const Edge& edge: LE)
            {
                adj_[edge.first].emplace_back(edge.second);
                adj_[edge.second].emplace_back(edge.first);
            }      
        }

        long long getNumberVertices()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            return graph_size_;
        }

        long long getNumberEdges()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (number_edges_ == -1)
            {
                number_edges_ = 0;
                for (auto& i: adj_)
                {
                    number_edges_ += i.size();
                }
                number_edges_ /= 2;
            }
            return number_edges_;
        }

        double getGraphDensity()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            return ((double) getNumberEdges() * 2) / ((double) (getNumberVertices()) * (getNumberVertices() - 1));
        }

        void clacEdges()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif
            
            if (getNumberEdges() != 0 && set_edges_.empty())
            {
                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    for (Vertex j: adj_[i])
                    {
                        set_edges_.insert({std::min(i, j), std::max(i, j)});
                    }
                }
            }
        }

        Edge getEdgeKey(Edge edge)
        {
            return {std::min(edge.first, edge.second), std::max(edge.first, edge.second)};
        }

        bool inGraph(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);

            clacEdges();
            if (set_edges_.find(edge) == set_edges_.end())
            {
                return false;
            }
            return true;
        }

        void calcWeakComponent()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (weak_components_.empty())
            {
                std::vector<bool> visited(getNumberVertices());

                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    if (!visited[i])
                    {
                        weak_components_.emplace_back();
                        std::vector<Vertex>& weak_component = weak_components_[weak_components_.size() - 1];
                        dfs(visited, i, [&weak_component](Vertex vertex){weak_component.emplace_back(vertex);});
                    }
                    std::sort(weak_components_[weak_components_.size() - 1].begin(), weak_components_[weak_components_.size() - 1].end());
                }
            }
            std::sort(weak_components_.begin(), weak_components_.end(), [](std::vector<Vertex>& v1, std::vector<Vertex>& v2){return v1.size() <= v2.size();});
        }

        long long getNumWeakComponents()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif
            
            calcWeakComponent();
            return weak_components_.size();
        }

        long long getMaxComponentSize()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif
            
            calcWeakComponent();
            return weak_components_[weak_components_.size() - 1].size();
        }

        std::vector<Vertex>& extractMaxWeakComponent()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcWeakComponent();
            return weak_components_[weak_components_.size() - 1];
        }

        void calcDistance(Vertex v)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (dists_.get(v) == nullptr)
            {
                std::vector<std::pair<Vertex, long long>> dists;
                std::vector<long long> tmp_dists(getNumberVertices(), -1);
                bfs(tmp_dists, v, [&dists, &tmp_dists](Vertex vertex){ dists.emplace_back(vertex, tmp_dists[vertex]);});

                dists_.insert(v, dists);
            }
        }

        std::vector<std::pair<Vertex, long long>>& getDistances(Vertex vertex)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcDistance(vertex);
            
            return *dists_.get(vertex);
        }

        long long getDistance(Vertex v1, Vertex v2)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (dists_.get(v2) != nullptr)
            {
                for (auto& i: *dists_.get(v2))
                {
                    if (i.first == v1)
                    {
                        return i.second;
                    }
                }
                return -1;
            }

            calcDistance(v1);
            for (auto& i: *dists_.get(v1))
            {
                if (i.first == v2)
                {
                    return i.second;
                }
            }
            return -1;
        }

        long long getEccentricity(Vertex vertex)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (eccentricities_.empty())
            {
                eccentricities_.assign(getNumberVertices(), -1);
            }
            if (eccentricities_[vertex] == -1)
            {
                long long res = 0;
                calcDistance(vertex);
                std::vector<std::pair<Vertex, long long>>& dists = *dists_.get(vertex);
                for (const auto& i: dists)
                {
                    res = std::max(res, i.second);
                }
                eccentricities_[vertex] = res;
            }
            return eccentricities_[vertex];
        } 

        long long getRadius()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (radius_ == -1)
            {
                if (getNumberVertices() > 500)
                {
                    calcRandomSubgraph();
                    radius_ = random_subgraph_->getRadius();
                }
                else
                {
                    long long res = LLONG_MAX;
                    for (Vertex i = 0; i < getNumberVertices(); ++i)
                    {
                        res = std::min(res, getEccentricity(i));
                    }
                    radius_ = res;
                }
            }
            return radius_;
        }

        long long getDiameter()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (diameter_ == -1)
            {
                if (getNumberVertices() > 1000)
                {
                    calcRandomSubgraph();
                    diameter_ = random_subgraph_->getDiameter();
                }
                else
                {
                    long long res = 0;
                    for (size_t i = 0; i < getNumberVertices(); ++i)
                    {
                        res = std::max(res, getEccentricity(i));
                    }
                    diameter_ = res;
                }
            }
            return diameter_;
        }

        AdjacencyList getSubgraph(const std::vector<Vertex>& vertices)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            std::set<Vertex> ver;
            std::vector<Vertex> rev(graph_size_);
            for (Vertex i = 0; i < vertices.size(); ++i)
            {
                ver.insert(vertices[i]);
                rev[vertices[i]] = i;
            }
            std::vector<std::vector<Vertex>> graph_data(vertices.size());
            for (Vertex i = 0; i < vertices.size(); ++i)
            {
                for (long long j = 0; j < adj_[vertices[i]].size(); ++j)
                {
                    if (ver.find(adj_[vertices[i]][j]) != ver.end())
                    {
                        graph_data[i].emplace_back(rev[adj_[vertices[i]][j]]);
                    }
                }
            }
            return graph_data;
        }

        long long getPercentile(double percent)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (getNumberVertices() > 1000)
            {
                calcRandomSubgraph();
                return random_subgraph_->getPercentile(percent);
            }
            if (percentile_count_.empty())
            {
                percentile_count_.assign(getNumberVertices(), 0);
                number_existing_paths_ = 0;
                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    calcDistance(i);
                    std::vector<std::pair<Vertex, long long>>& dists = *dists_.get(i);
                    for (auto& j: dists)
                    {
                        if (j.first > i)
                        {
                            ++number_existing_paths_;
                            ++percentile_count_[j.second];
                        }
                    }
                }
            }
            long long num = std::floor(((double)number_existing_paths_ * (percent / 100)) + 0.5);
            long long tmp = 0;
            if (num == 0)
            {
                return 0;
            }
            for (long long i = 1; i < getNumberVertices(); ++i)
            {
                tmp += percentile_count_[i];
                if (tmp >= num)
                {
                    return i;
                }
            }
            return getNumberVertices() - 1;
        }

        std::vector<Vertex> extractRandomSubgraph(long long graph_size)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (graph_size >= getNumberVertices()) 
            {
                std::vector<Vertex> vertices(getNumberVertices());
                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    vertices[i] = i;
                }
                return vertices;
            }

            std::vector<Vertex> vertices(getNumberVertices());
            for (size_t i = 0; i < getNumberVertices(); ++i)
            {
                vertices[i] = i;
            }
            std::shuffle(vertices.begin(), vertices.end(), std::default_random_engine());
            vertices.resize(graph_size);
            std::sort(vertices.begin(), vertices.end());
            return vertices;
        }

        std::vector<Vertex>& extractSnowballSampleSubgraph(Vertex start, long long graph_size = 100)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (local_extracted_subgraphs_.get(start) == nullptr)
            {
                graph_size = std::min(getNumberVertices(), graph_size);
                std::vector<std::pair<long long, Vertex>> tmp;
                calcDistance(start);
                std::vector<std::pair<Vertex, long long>>& dists = *dists_.get(start);
                for (auto& j: dists)
                {
                    tmp.emplace_back(j.second, j.first);
                }

                std::sort(tmp.begin(), tmp.end());
                bool flag = false;
                long long j = 0;
                local_extracted_subgraphs_.insert(start, {});

                for (long long i = 0; i < graph_size - 1; ++i)
                {
                    if (j == tmp.size())
                    {
                        break;
                    }
                    if (!flag && start < tmp[j].second)
                    {
                        local_extracted_subgraphs_.get(start)->emplace_back(start);
                        flag = true;
                        continue;
                    }
                    local_extracted_subgraphs_.get(start)->emplace_back(tmp[j].second);
                    ++j;
                }
                if(!flag)
                {
                    local_extracted_subgraphs_.get(start)->emplace_back(start);
                }
            }
            
            return *local_extracted_subgraphs_.get(start);
        }

        std::vector<Vertex> extractRundomSnowballsSubgraph(long long max_size)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (max_size >= getNumberVertices()) 
            {
                std::vector<Vertex> vertices(getNumberVertices());
                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    vertices[i] = i;
                }
                return vertices;
            }

            std::random_device rd;
            std::mt19937 gen(rd());
            long long max_num_iteration = max_size * 10;
            long long loc_size = std::sqrt(max_size);
            long long max_num_extraction = std::pow(max_size, 2. / 3.);

            std::set<Vertex> tmp;

            while(tmp.size() < max_size && max_num_iteration > 0 && max_num_extraction > 0)
            {
                Vertex start = 0;
                do
                {
                    if (max_num_iteration == 0)
                    {
                        break;
                    }
                    start = std::uniform_int_distribution<Vertex>(0, getNumberVertices() - 1)(gen);
                    --max_num_iteration;
                } 
                while (tmp.find(start) != tmp.end());
                
                if (max_num_iteration == 0)
                {
                    break;
                }
                long long tmp_loc_size = 0;
                for (auto& i: extractSnowballSampleSubgraph(start))
                {
                    if (tmp_loc_size > loc_size)
                    {
                        break;
                    }
                    tmp.insert(i);
                    ++tmp_loc_size;
                }
                --max_num_extraction;
            }
            std::vector<Vertex> res;
            for (auto& i: tmp)
            {
                res.emplace_back(i);
            }
            return res;
        }

        void calcR()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (R1_ == -1)
            {
                R1_ = R2_ = R3_ = Re_ = 0;
                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    R1_ += adj_[i].size();
                    R2_ += adj_[i].size() * adj_[i].size();
                    R3_ += adj_[i].size() * adj_[i].size() * adj_[i].size();
                    for (Vertex j: adj_[i])
                    {
                        Re_ += adj_[i].size() * adj_[j].size();
                    }
                }
            }
        }

        long long getR1()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcR();
            return R1_;
        }

        long long getR2()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcR();
            return R2_;
        }

        long long getR3()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcR();
            return R3_;
        }

        long long getRe()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcR();
            return Re_;
        }

        void calcSheredAdjacentVertices(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            if (shared_adjacent_vertices_.get(edge) == nullptr)
            {
                shared_adjacent_vertices_.insert(edge, getVertexIntersection(adj_[edge.first], adj_[edge.second]));
            }
        }

        std::vector<Vertex>& getSharedAdjacentVertices(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            calcSheredAdjacentVertices(edge);

            return *shared_adjacent_vertices_.get(edge);
        }

        long long getNumSharedAdjacentVertices(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            return getSharedAdjacentVertices(edge).size();
        }

        void calcAdjacentVertices(Edge edge)
        {
            edge = getEdgeKey(edge);
            if (adjacent_vertices_.get(edge) == nullptr)
            {
                adjacent_vertices_.insert(edge, getVertexCombining(adj_[edge.first], adj_[edge.second]));
            }
        }

        std::vector<Vertex>& getAdjacentVertices(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            calcAdjacentVertices(edge);

            return *adjacent_vertices_.get(edge);
        }

        long long getNumdAdjacentVertices(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            return getAdjacentVertices(edge).size();
        }

        long long getL(Vertex vertex)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (L_.empty())
            {
                L_.assign(getNumberVertices(), -1);
            }
            if (L_[vertex] == -1) 
            {
                L_[vertex] = 0;
                for (auto& i: adj_[vertex])
                {
                    L_[vertex] += getNumSharedAdjacentVertices({vertex, i});
                }
                L_[vertex] /= 2;
            }
            return L_[vertex];
        }

        double getClusteringCoefficient()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (clustering_coefficient_ == -1)
            {
                clustering_coefficient_ = 0;
                for (Vertex i = 0; i < getNumberVertices(); ++i)
                {
                    if (adj_[i].size() >= 2)
                    {
                        clustering_coefficient_ += (2.L * getL(i)) / ((long double) (adj_[i].size()) * (adj_[i].size() - 1));
                    }
                }
                clustering_coefficient_ /= (double) getNumberVertices();
            }
            return clustering_coefficient_;
        }

        double getPearsonCorrelationCoefficient()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            calcR();
            if (R3_ * R1_ == R2_ * R2_)
            {
                return 0.L;
            }
            return ((double)Re_ * R1_ - (double)R2_ * R2_) / ((double)R3_ * R1_ - (double)R2_ * R2_);
        }

        void dfs(std::vector<bool>& visited, Vertex vertex, std::function<void(Vertex)> fun = [](Vertex){})
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (visited[vertex]) return;
            visited[vertex] = true;
            fun(vertex);
            for (auto& i: adj_[vertex])
            {
                dfs(visited, i, fun);
            }
        }

        void bfs(std::vector<long long>& distance, Vertex vertex, std::function<void(Vertex)> fun = [](Vertex){})
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (distance[vertex] != -1) return;
            std::queue<Vertex> que;
            que.push(vertex);
            distance[vertex] = 0;
            fun(vertex);
            while (!que.empty())
            {
                Vertex tmp = que.front();
                que.pop();
                for (auto& i: adj_[tmp])
                {
                    if (distance[i] == -1)
                    {
                        distance[i] = distance[tmp] + 1;
                        que.push(i);
                        fun(i);
                    }
                }
            }
        }

        long long getVertexDegree(Vertex v)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            return adj_[v].size();
        }

        std::vector<Vertex> getVertexCombining(const std::vector<Vertex>& vs1, const std::vector<Vertex>& vs2)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            long long i = 0, j = 0;
            std::vector<Vertex> res;
            while(i < vs1.size() || j < vs2.size())
            {
                if (i == vs1.size())
                {
                    res.emplace_back(vs2[j]);
                    ++j;
                }
                else if (j == vs2.size())
                {
                    res.emplace_back(vs1[i]);
                    ++i;
                }
                else 
                {
                    if (vs1[i] == vs2[j])
                    {
                        res.emplace_back(vs1[i]);
                        ++i;
                        ++j;
                    }
                    else if (vs1[i] < vs2[j])
                    {
                        res.emplace_back(vs1[i]);
                        ++i;
                    }
                    else
                    {
                        res.emplace_back(vs2[j]);
                        ++j;
                    }
                }
            }
            return res;
        }

        std::vector<Vertex> getVertexIntersection(const std::vector<Vertex>& vs1, const std::vector<Vertex>& vs2)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            long long i = 0, j = 0;
            std::vector<Vertex> res;
            while(i < vs1.size() && j < vs2.size())
            {
                if (vs1[i] == vs2[j])
                {
                    res.emplace_back(vs1[i]);
                    ++i;
                    ++j;
                }
                else if (vs1[i] < vs2[j])
                {
                    ++i;
                }
                else
                {
                    ++j;
                }
            }
            return res;
        }

        std::vector<Vertex>& getNeighbors(Vertex v)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            return adj_[v];
        }

        void calcLocalSubgraph(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);

            if (local_subgraph_.get(edge) == nullptr)
            {
                std::vector<Vertex> loc_subgraph = getVertexCombining(extractSnowballSampleSubgraph(edge.first), extractSnowballSampleSubgraph(edge.second));
                Edge new_edge = {0, 0};
                for (long long i = 0; i < loc_subgraph.size(); ++i)
                {
                    if (loc_subgraph[i] == edge.first)
                    {
                        new_edge.first = i;
                    }
                    if (loc_subgraph[i] == edge.second)
                    {
                        new_edge.second = i;
                    }
                }
                local_subgraph_.insert(edge, {new Graph(getSubgraph(loc_subgraph)), new_edge});
            }
        }

        std::pair<Graph&, Edge> getLocalSubgraph(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            calcLocalSubgraph(edge);
            auto res = local_subgraph_.get(edge);

            return {*res->first, res->second};
        }

        void calcRandomSubgraph(long long size = 1000)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (!random_subgraph_)
            {
                random_subgraph_ = new Graph(getSubgraph(extractRundomSnowballsSubgraph(size)));
            }
        }

        long long getNumberNeighbors(Vertex vertex)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            return adj_[vertex].size();
        }

        double getJacquardCoefficient(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);
            
            if (getNumberNeighbors(edge.first) + getNumberNeighbors(edge.second) - getNumSharedAdjacentVertices(edge) == 0)
            {
                return 1.;
            }
            return (double)getNumSharedAdjacentVertices(edge) / (getNumberNeighbors(edge.first) + getNumberNeighbors(edge.second) - getNumSharedAdjacentVertices(edge));
        }

        double getAdarCoefficient(Edge edge)
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            edge = getEdgeKey(edge);

            double adar_coefficient = 0.;
            for (auto& i: getVertexIntersection(getNeighbors(edge.first), getNeighbors(edge.second)))
            {
                adar_coefficient += 1. / std::log(getNumberNeighbors(i) + 2);
            }
            return adar_coefficient;
        }

        ~Graph()
        {
            #ifdef BENCHMARK
                _BENCHMARK_
            #endif

            if (random_subgraph_)
            {
                delete random_subgraph_;
            }
        }

    protected:

        long long graph_size_;
        long long R1_ = -1, R2_ = -1, R3_ = -1, Re_ = -1;
        long long radius_ = -1;
        long long diameter_ = -1;
        long long number_existing_paths_ = -1;
        double clustering_coefficient_ = -1;
        long long number_edges_ = -1;

        AdjacencyList adj_;

        std::vector<long long> L_;
        std::vector<long long> eccentricities_;
        std::vector<long long> percentile_count_;
        SetEdges set_edges_;

        BufInfo<Edge, std::vector<Vertex>> shared_adjacent_vertices_;
        BufInfo<Edge, std::vector<Vertex>> adjacent_vertices_;
        BufInfo<Vertex, std::vector<std::pair<Vertex, long long>>> dists_;
        BufInfo<Vertex, std::vector<Vertex>> local_extracted_subgraphs_;
        BufInfo<Edge, std::pair<Graph*, Edge>> local_subgraph_;

        Graph* random_subgraph_ = nullptr;
        std::vector<std::vector<Vertex>> weak_components_;
        
    };


    class GraphTimestamp: virtual public Graph
    {
    public:
        using EdgeTS = std::pair<double, Edge>;
        using AdjacencyListTS = std::vector<std::vector<std::pair<Vertex, double>>>;
        using AdjacencyMatrixTS = std::vector<std::vector<double>>;
        using ListEdgesTS = std::vector<EdgeTS>; 

        GraphTimestamp(long long graph_size, const ListEdges& LE, const std::vector<double>& TS): Graph({}, graph_size), time_(0), LETS_(LE.size())
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            for (long long i = 0; i < LE.size(); ++i)
            {
                LETS_[i] = {TS[i], LE[i]};
                edges_time_[LE[i]] = TS[i];
            }
            std::sort(LETS_.begin(), LETS_.end());
        }

        void clear()
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            adj_.clear();
            number_edges_ = -1;
            weak_components_.clear();
            R1_ = R2_ = R3_ = Re_ = -1;
            L_.clear();
            radius_ = diameter_ = -1;
            eccentricities_.clear();
            percentile_count_.clear();
            number_existing_paths_ = clustering_coefficient_ = -1;
            set_edges_.clear();
            shared_adjacent_vertices_.clear();
            adjacent_vertices_.clear();
            dists_.clear();
            local_extracted_subgraphs_.clear();
            local_subgraph_.clear();
            
            if (random_subgraph_)
            {
                delete random_subgraph_;
            }
            random_subgraph_ = nullptr;
        }

        void initTime(double init_time)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            clear();
            
            long long i = 0;
            adj_.assign(getNumberVertices(), std::vector<Vertex>());
            while(i < LETS_.size() && init_time >= LETS_[i].first)
            {
                adj_[LETS_[i].second.first].emplace_back(LETS_[i].second.second);
                adj_[LETS_[i].second.second].emplace_back(LETS_[i].second.first);
                ++i;
            }
            time_ = init_time;
        }

        double getEdgeTime(Edge edge)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            if (edges_time_.find(edge) == edges_time_.end())
            {
                return -1.;
            }
            return edges_time_.at(edge);
        }

        double getFirstEdgeTime()
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            return LETS_[0].first;
        }

        double getLastEdgeTime()
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            return LETS_[LETS_.size() - 1].first;
        }

        double getNormTime(double time)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            return (time - getFirstEdgeTime()) / 1e10;
        }

        double getNormEdgeTime(Edge edge)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            double edge_time = getEdgeTime(edge);
            if (edge_time < 0)
            {
                return edge_time;
            }
            return getNormTime(edge_time);
        }

        double getNormTimePercent(double percent)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            double time = LETS_[size_t((LETS_.size() - 1) * percent / 100.)].first;
            return getNormTime(time);
        }

        void initRandomTime()
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            std::random_device rd;
            std::mt19937 gen(rd());

            long double new_time = LETS_[std::uniform_int_distribution<size_t>(0, LETS_.size() - 1)(gen)].first;

            initTime(new_time);
        }

        long long getFullNumEdges()
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            return LETS_.size();
        }

    protected:

    std::map<Edge, double> edges_time_;
    double time_;
    ListEdgesTS LETS_;
    };


    class GraphFeatureEdges: virtual public Graph
    {
    public:
        using EdgeWithFeatures = std::pair<Edge, std::vector<double>>;
        using FeatureVector = std::vector<double>;
        using FeatureMatrix = std::vector<FeatureVector>;

        GraphFeatureEdges (
            long long graph_size, 
            const ListEdges& LE, 
            std::vector<std::function<double(Graph&, Edge)>> funs): 
            Graph(LE, graph_size), funs_(funs) {}

        FeatureVector getFeatureVector(Edge edge)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            FeatureVector fv(funs_.size());

            for (long long i = 0; i < funs_.size(); ++i)
            {
                    fv[i] = funs_[i](*this, edge);
            }

            return fv;
        }

        std::pair<FeatureMatrix, ListEdges> getFeatureVectors(const ListEdges& edges)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            FeatureMatrix resMaxtrix(edges.size(), FeatureVector(funs_.size()));
            ListEdges resEdges(edges.size());
            for (long long i = 0; i < edges.size(); ++i)
            {
                resMaxtrix[i] = getFeatureVector(edges[i]);
                resEdges[i] = edges[i];
            }
            return {resMaxtrix, resEdges};
        }

        std::pair<FeatureMatrix, ListEdges> getRandomFeatureVectors(long long num)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            std::random_device rd;
            std::mt19937 gen(rd());

            Vertex v1, v2;
            std::set<Edge> tmp;
            ListEdges edges;
            long long tmp_number_edges = getNumberEdges();

            while (num > 0)
            {
                if (tmp_number_edges == (getNumberVertices()) * (getNumberVertices() - 1) / 2)
                {
                    break;
                }
                v1 = std::uniform_int_distribution<size_t>(0, graph_size_ - 2)(gen);
                v2 = std::uniform_int_distribution<size_t>(v1 + 1, graph_size_ - 1)(gen);
                Edge edge = {v1, v2};
                if (!inGraph(edge) && tmp.find(edge) == tmp.end())
                {
                    tmp.insert(edge);
                    ++tmp_number_edges;
                    edges.emplace_back(edge);
                    --num;
                }
            }

            return getFeatureVectors(edges);
        }

    protected:
        std::vector<std::function<double(Graph&, Edge)>> funs_;

    };


    class GraphFeatureEdgesTimestamp: public GraphTimestamp, public GraphFeatureEdges
    {
    public:
        GraphFeatureEdgesTimestamp (
            long long graph_size, 
            const ListEdges& LE,  
            const std::vector<double>& TS, 
            std::vector<std::function<double(Graph&, Edge)>> funs): 
                GraphTimestamp(graph_size, LE, TS),
                GraphFeatureEdges(graph_size, LE, funs),
                Graph(LE, graph_size) {}

        ListEdges extractListEdges(
            double init_time, 
            double predict_time, 
            long long true_num, 
            long long in_false_num, 
            long long out_false_num)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            ListEdges true_edges, in_false_edges, out_false_edges;

            for (long long i = 0; i < LETS_.size(); ++i)
            {
                if (LETS_[i].first > init_time && LETS_[i].first <= predict_time)
                {
                    true_edges.emplace_back(LETS_[i].second);
                }
                if (LETS_[i].first > predict_time)
                {
                    in_false_edges.emplace_back(LETS_[i].second);
                }
            }
            std::shuffle(true_edges.begin(), true_edges.end(), std::default_random_engine());
            std::shuffle(in_false_edges.begin(), in_false_edges.end(), std::default_random_engine());

            if (true_edges.size() > true_num)
            {
                true_edges.resize(true_num);
            }
            if (in_false_edges.size() > in_false_num)
            {
                in_false_edges.resize(in_false_num);
            }

            std::random_device rd;
            std::mt19937 gen(rd());
            std::set<Edge> false_edges;
            long long max_iter = out_false_num * 10;

            for (long long i = 0; i < out_false_num + in_false_num - in_false_edges.size(); ++i)
            {
                Edge edge;
                long double edge_time;
                do
                {
                    if (max_iter == 0)
                    {
                        break;
                    }
                    Vertex v1 = std::uniform_int_distribution<size_t>(0, getNumberVertices() - 2)(gen);
                    Vertex v2 = std::uniform_int_distribution<size_t>(v1 + 1,  getNumberVertices() - 1)(gen);
                    edge = {v1, v2};
                    edge_time = getEdgeTime(edge);
                    --max_iter;
                }
                while(edge_time != -1.L || false_edges.find(edge) != false_edges.end());
                
                if (max_iter == 0)
                {
                    break;
                }

                false_edges.insert(edge);
                out_false_edges.emplace_back(edge);
            }

            for (auto& i: in_false_edges)
            {
                true_edges.emplace_back(i);
            }

            for (auto& i: out_false_edges)
            {
                true_edges.emplace_back(i);
            }

            return true_edges;
        }

        std::vector<std::vector<double>> getData(double new_time, const ListEdges& edges)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif
            
            ListEdges res_edges;
            std::vector<std::set<Vertex>> tmp_adj(getNumberVertices());

            for (auto& i: edges)
            {
                tmp_adj[i.first].insert(i.second);
                tmp_adj[i.second].insert(i.first);
            }

            std::set<std::pair<long long, Vertex>> fr;
            for (long long i = 0; i < getNumberVertices(); ++i)
            {
                if (tmp_adj[i].size() > 0)
                {
                    fr.insert({-(long long)tmp_adj[i].size(), i});
                }
            }

            while(!fr.empty())
            {
                auto[fr_num, vertex] = *fr.begin();
                fr.erase(fr.begin());
                for (auto& j: tmp_adj[vertex])
                {
                    res_edges.emplace_back(vertex, j);
                    fr.erase({-(long long)tmp_adj[j].size(), j});
                    tmp_adj[j].erase(vertex);
                    if (tmp_adj[j].size() > 0)
                    {
                        fr.insert({-(long long)tmp_adj[j].size() + 1, j});
                    }
                }
            }

            if (time_ != new_time)
            {
                initTime(new_time);
            }
            double norm_time = getNormTime(new_time);
            auto[matrix, list_edges] = getFeatureVectors(res_edges);
            for (size_t i = 0; i < matrix.size(); ++i)
            {
                matrix[i].emplace_back(norm_time);
                matrix[i].emplace_back(getNormEdgeTime(list_edges[i]));
            }
            return matrix;
        }
    
        std::vector<std::vector<double>> getData(
            double init_time, 
            double predict_time, 
            long long true_num, 
            long long in_false_num, 
            long long out_false_num)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            return getData(init_time, extractListEdges(init_time, predict_time, true_num, in_false_num, out_false_num));
        }
    
        std::vector<std::vector<double>> getDataPart(
            double init_percent, 
            double predict_percent, 
            long long true_num, 
            long long in_false_num, 
            long long out_false_num)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            double init_time = LETS_[size_t((LETS_.size() - 1) * init_percent / 100.)].first;
            double predict_time = LETS_[size_t((LETS_.size() - 1) * predict_percent / 100.)].first;
            return getData(init_time, predict_time, true_num, in_false_num, out_false_num);
        }
    
        std::vector<std::vector<double>> getDataWithCalc(
            double init_percent, 
            double predict_percent, 
            long long true_num, 
            long long in_false_num, 
            long long out_false_num)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            double init_time = LETS_[size_t((LETS_.size() - 1) * init_percent / 100.)].first;
            double predict_time = LETS_[size_t((LETS_.size() - 1) * predict_percent / 100.)].first;

            ListEdges edges = extractListEdges(init_time, predict_time, true_num, in_false_num, out_false_num);

            long long iterations = std::sqrt(edges.size()) + 1;
            
            std::vector<ListEdges> parts_edges;
            std::vector<std::vector<double>> res;

            for (long long i = 0; i < iterations; ++i)
            {
                parts_edges.emplace_back();
                for (long long j = 0; j < iterations; ++j)
                {
                    if (i * iterations + j >= edges.size())
                    {
                        break;
                    }
                    parts_edges[parts_edges.size() - 1].emplace_back(edges[i * iterations + j]);
                }
            }

            auto allbegin = std::chrono::steady_clock::now();
            double norm_predict_time = getNormTime(predict_time);

            for (long long i = 0; i < parts_edges.size(); ++i)
            {
                auto begin = std::chrono::steady_clock::now();
                for (auto& i: getData(init_time, parts_edges[i]))
                {
                    res.emplace_back(i);
                    res[res.size() - 1].emplace_back(norm_predict_time);
                    if (res[res.size() - 1][res[res.size() - 1].size() - 2] != -1. && 
                        res[res.size() - 1][res[res.size() - 1].size() - 2] <= norm_predict_time)
                    {
                        res[res.size() - 1].emplace_back(1.);
                    }
                    else
                    {
                        res[res.size() - 1].emplace_back(0.);
                    }
                }
                auto end = std::chrono::steady_clock::now();
                std::chrono::seconds elapsed_ms = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
                std::chrono::seconds all_elapsed_ms = std::chrono::duration_cast<std::chrono::seconds>(end - allbegin);
                long double per = (long double)(i + 1) * 100.L / ((long double)(iterations));
                std::cout << per << "%" << " ---- Complete | " << elapsed_ms.count() << "s. | Total: " << all_elapsed_ms.count() << "s." << std::endl;
            }

            return res;
        }

        std::vector<std::vector<double>> getDataPartWithCalc(
            double init_percent, 
            double predict_percent, 
            long long true_num, 
            long long in_false_num, 
            long long out_false_num)
        {
            // #ifdef BENCHMARK
            //     _BENCHMARK_
            // #endif

            long long init_num = (LETS_.size() - 1) * init_percent / 100.;
            long long predict_num = (LETS_.size() - 1) * predict_percent / 100.;

            long long iterations = std::sqrt(true_num + in_false_num + out_false_num) + 1;

            double step_num = (double)(predict_num - init_num) / (double)iterations;

            std::vector<std::vector<double>> res;

            auto allbegin = std::chrono::steady_clock::now();

            for (long long i = 0; i < iterations; ++i)
            {
                auto begin = std::chrono::steady_clock::now();
                double tmp_init_time = LETS_[size_t(init_num + step_num * i)].first;
                double tmp_predict_time = LETS_[size_t(init_num + step_num * (i + 1))].first;
                for (auto& j: getData(tmp_init_time, extractListEdges(
                    tmp_init_time,
                    tmp_predict_time,
                    true_num / iterations + 1,
                    in_false_num / iterations + 1,
                    out_false_num / iterations + 1)))
                {
                    double norm_time = getNormTime(tmp_predict_time);
                    res.emplace_back(j);
                    res[res.size() - 1].emplace_back(norm_time);
                    if (res[res.size() - 1][res[res.size() - 1].size() - 2] != -1. && 
                        res[res.size() - 1][res[res.size() - 1].size() - 2] <= norm_time)
                    {
                        res[res.size() - 1].emplace_back(1.);
                    }
                    else
                    {
                        res[res.size() - 1].emplace_back(0.);
                    }

                }
                auto end = std::chrono::steady_clock::now();
                std::chrono::seconds elapsed_ms = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
                std::chrono::seconds all_elapsed_ms = std::chrono::duration_cast<std::chrono::seconds>(end - allbegin);
                long double per = (long double)(i + 1) * 100.L / ((long double)(iterations));
                std::cout << per << "%" << " ---- Complete | " << elapsed_ms.count() << "s. | Total: " << all_elapsed_ms.count() << "s." << std::endl;
            }

            #ifdef BENCHMARK
                Benchmark::printRes();
            #endif

            return res;
        }
    };

}
