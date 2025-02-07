#ifndef DE_BRUIJN_H
#define DE_BRUIJN_H

#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

struct read {
    //std::string label;
    std::string seq;
    //std::string quality;
};

struct contig {
    std::string label;
    std::string seq;
};

struct edge {
    unsigned int sourceId;
    unsigned int destinationId;

    bool operator<(const edge& a) const {
        if (sourceId != a.sourceId)
            return sourceId < a.sourceId;
        return destinationId < a.destinationId;
    }

    bool operator==(const edge& a) const {
        return sourceId == a.sourceId && destinationId == a.destinationId;
    }
};

template <>
struct std::hash<edge> {
    std::size_t operator()(const edge& e) const noexcept {
        return (std::hash<unsigned int>()(e.sourceId) ^ (std::hash<unsigned int>()(e.destinationId) << 1));
    }
};

// Get time since the program started as a string in a nice format ( [hr:min:sec.ms] )
std::string elapsedTime(const std::chrono::steady_clock::time_point&);

class DeBruijnGraph {
private:
    std::unordered_map<std::string, unsigned int> vertices;
    std::unordered_map<edge, unsigned int> edges;
    std::unordered_map<unsigned int, std::vector<unsigned int>> al;
    
    void depthFirstSearch(unsigned int, std::unordered_set<unsigned int>&, std::vector<unsigned int>&);
    std::vector<unsigned int> traceEulerianPath(unsigned int);
public:
    unsigned int k;
    DeBruijnGraph(const std::vector<read>&, unsigned int, std::chrono::steady_clock::time_point);
};

#endif