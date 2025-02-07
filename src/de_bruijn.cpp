#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <filesystem>
#include <fstream>
#include <de_bruijn.hpp>

DeBruijnGraph::DeBruijnGraph(const std::vector<read>& reads, unsigned int k, std::chrono::steady_clock::time_point st) {
    (*this).k = k;

    std::ostringstream dirPathSS;
    dirPathSS << k << "mer";
    std::string dirPath = dirPathSS.str();
    
    if (!std::filesystem::exists(dirPath)) {
        std::cout << "  --> Creating directory " << dirPath << std::endl;
        if (!std::filesystem::create_directory(dirPath))
            throw std::runtime_error("Could not create directory");
    } else {
        std::cout << "  --> Directory " << dirPath << " already exists, wiping previous data... " << std::endl;
        for (const auto& entry : std::filesystem::directory_iterator(dirPath)) {
            if (std::filesystem::is_regular_file(entry)) { // Only delete regular files
                std::filesystem::remove(entry);
                std::cout << "      * Deleted file: " << entry.path() << std::endl;
            }
        }
    }

    std::unordered_map<std::string, unsigned int>::iterator searchV;
    std::unordered_map<edge, unsigned int>::iterator searchN;
    std::string km1L;
    std::string km1R;
    unsigned int km1LId;
    unsigned int km1RId;
    unsigned int km1RIdPrev;
    edge newEdge;
    
    float splitPointPercent = 0.0f;
    std::size_t splitPointReads = 0.0f;
    std::chrono::steady_clock::time_point startSplitPoint = std::chrono::steady_clock::now();

    for (unsigned int i = 0; i < reads.size(); i++) {
        unsigned int readLen = reads[i].seq.length();

        if (readLen < k) {
            throw std::runtime_error("K-mer size is larger than the read length");
        }

        // split read into k-1 mers
        for (unsigned int j = 0; j < readLen - k + 1; j += 2) {
            km1L = reads[i].seq.substr(j, k - 1); // (k - 1)-mer left
            km1R = reads[i].seq.substr(j + 1, k - 1); // (k - 1)-mer right

            searchV = (*this).vertices.find(km1L);
            if (searchV == (*this).vertices.end()) {
                km1LId = (*this).vertices.size();
                (*this).vertices[km1L] = (*this).vertices.size();
            } else {
                km1LId = searchV->second;
            }

            searchV = (*this).vertices.find(km1R);
            if (searchV == (*this).vertices.end()) {
                km1RId = (*this).vertices.size();
                (*this).vertices[km1R] = (*this).vertices.size();
            } else {
                km1RId = searchV->second;
            }

            // Make edge between previous km1RId and current km1LId
            if (j > 0) {
                newEdge = {.sourceId = km1RIdPrev, .destinationId = km1LId};
                searchN = (*this).edges.find(newEdge);

                if (searchN == (*this).edges.end()) {
                    (*this).edges[newEdge] = 1;
                } else {
                    (*this).edges[newEdge] += 1;
                }
            }

            
            // Turn km1LId and km1RId into edge
            newEdge = {.sourceId = km1LId, .destinationId = km1RId};
            searchN = (*this).edges.find(newEdge);
            
            if (searchN == (*this).edges.end()) {
                (*this).edges[newEdge] = 1;
            } else {
                (*this).edges[newEdge] += 1;
            }
            
            km1RIdPrev = km1RId;
        }

        // Every 5% of reads read a status update
        if (i == splitPointReads) {
            std::cout << "  --> " << splitPointPercent * 100 << "% of reads processed";

            if (splitPointPercent == 0.0f) {
                std::cout << " (est. ?? hr ?? min ?? sec remaining)" << std::endl;
            }
            else {
                std::chrono::seconds::rep secElapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - startSplitPoint).count();
                std::chrono::seconds::rep secEstTimeRemaining = (secElapsed / splitPointPercent) * (1.0f - splitPointPercent);    // Estimated seconds remaining

                std::chrono::hours::rep hr = secEstTimeRemaining / 3600;
                secEstTimeRemaining -= 3600 * hr;
                std::chrono::minutes::rep min = secEstTimeRemaining / 60;
                secEstTimeRemaining -= 60 * min;

                std::cout << " (est. " << hr << " hr " << min << " min " << secEstTimeRemaining << " sec remaining)" << std::endl;
            }
            
            splitPointPercent += 0.05f;
            splitPointReads = static_cast<std::size_t>(std::round(static_cast<float>(reads.size()) * splitPointPercent));
        }
    }

    std::cout << "  --> 100% of reads have finished processing!\t" << (*this).vertices.size() << " unique (k-1)-mers" << std::endl;
    std::cout << "  --> Getting (k-1)-mers data ready to save to file..." << std::endl;

    // Write map of vertices to file where the value at line number = the key
    std::string mapFilePath = dirPath + "/" + "map.txt";
    unsigned int maxLine = (*this).vertices.size();

    std::ofstream mapFileOut(mapFilePath);
    if (!mapFileOut.is_open())
        throw std::runtime_error("Could not open map file to write to");

    std::vector<std::string> lines(maxLine - 1, ""); 

    for (auto it = (*this).vertices.begin(); it != (*this).vertices.end();) {
        if (it->second < lines.size())
            lines[it->second] = it->first;
        it = (*this).vertices.erase(it);
    }

    std::cout << "  --> Writing (k-1)-mers data to file..." << std::endl;

    unsigned int lineCntr = 0;

    for (auto& line : lines) {
        mapFileOut << line << "\n";

        if (lineCntr % 5000 == 0) {
            mapFileOut.flush();
        }

        line = "";
        ++lineCntr;
    }

    mapFileOut.close();
    lines.clear();

    std::cout << "  --> Assembling adjacency graph from edges" << std::endl;

    // Build adjacency graph from edges
    std::unordered_map<unsigned int, std::vector<unsigned int>>::iterator searchAL;

    unsigned int s; // source id
    unsigned int d; // destination id

    for (auto it = (*this).edges.begin(); it != (*this).edges.end(); ) {
        s = it->first.sourceId;
        d = it->first.destinationId;

        //std::cout << s << " --> " << d << std::endl;

        searchAL = (*this).al.find(s);
        if (searchAL == (*this).al.end()) {
            (*this).al[s] = std::vector<unsigned int>();
        }

        searchAL = (*this).al.find(d);
        if (searchAL == (*this).al.end()) {
            (*this).al[d] = std::vector<unsigned int>();
        }

        (*this).al[s].push_back(d);
        it = (*this).edges.erase(it);
    }

    std::cout << "  --> Identifying disjoint subgraphs" << std::endl;

    // Identify disjoint subgraphs; each of these are contigs
    std::unordered_set<unsigned int> visited; // Use unordered_set for faster lookups
    std::vector<std::vector<unsigned int>> components;

    if ((*this).al.empty())
        throw std::runtime_error("Adjacency list is empty. Exiting...");

    // Perform DFS on every unvisited node
    for (auto it = (*this).al.begin(); it != (*this).al.end(); ++it) {
        if (visited.find(it->first) == visited.end()) {
            std::vector<unsigned int> component;
            (*this).depthFirstSearch(it->first, visited, component);
            components.push_back(component);
        }
    }

    unsigned int contigCntr = 0;

    std::ifstream mapFileIn(mapFilePath);

    std::string contigFilePath = dirPath + "/" + "contigs.fasta";

    // Process each disjoint subgraph (component) to reconstruct contigs
    std::cout << "  --> Assembling contigs from disjoint subgraphs (this may take a while)" << std::endl;
    for (const auto& component : components) {
        if (component.empty()) continue;

        std::ofstream contigFile(contigFilePath, std::ios::app);

        if (!contigFile.is_open())
            throw std::runtime_error("Could not open output file");

        // Find a starting node for the Eulerian path in the current component
        unsigned int startNode = component[0];
        for (unsigned int node : component) {
            if (al[node].size() > 0) {
                startNode = node;
                break;
            }
        }

        
        std::vector<unsigned int> eulerianPath = traceEulerianPath(startNode);
        std::ostringstream contigStream;
        std::string contigStr;
        std::ostringstream contigLabel;

        if (eulerianPath.size() > 1) {
            for (auto it = eulerianPath.rbegin(); it != eulerianPath.rend(); ++it) {
                std::string km1mer;

                // get km1mer from file mapFile
                mapFileIn.seekg((*it) * ((*this).k), std::ios::beg);
                std::getline(mapFileIn, km1mer);

                if (km1mer.length() > 0) {
                    if (it == eulerianPath.rbegin()) contigStream << km1mer;
                    else contigStream << km1mer.back();
                }
            }

            contigStr = contigStream.str();
            contigStr.erase(std::remove(contigStr.begin(), contigStr.end(), '\0'), contigStr.end());

            if (contigStr.length() > 0) {
                ++contigCntr;
                contigLabel << ">contig." << contigCntr;
                
                // write contig to file
                contigFile << contigLabel.str() << "\n" << contigStr << std::endl;
            }

        }
    
        contigFile.close();
    }

    mapFileIn.close();
    std::cout << "  --> Assembled " << contigCntr << " contigs." << std::endl;
}

// Get time since the program started as a string in a nice format ( [hr:min:sec.ms] )
std::string elapsedTime(const std::chrono::steady_clock::time_point& st) {
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    std::chrono::milliseconds::rep ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - st).count();
    
    std::chrono::hours::rep hr = ms / 3600000;
    ms -= 3600000 * hr;

    std::chrono::minutes::rep min = ms / 60000;
    ms -= 60000 * min;

    std::chrono::seconds::rep sec = ms / 1000;
    ms -= 1000 * sec;

    std::ostringstream output;
    output << "["
        << hr << ":"
        << std::setw(2) << std::setfill('0') << min << ":"
        << std::setw(2) << std::setfill('0') << sec << "."
        << std::setw(3) << std::setfill('0') << ms
        << "]";

    return output.str();
}

// Recursive DFS kept throwing a seg fault, I am thinking it was because of stack size limit being exceeded on windows????
void DeBruijnGraph::depthFirstSearch(unsigned int start_node, std::unordered_set<unsigned int>& visited, std::vector<unsigned int>& component) {
    std::stack<unsigned int> stack; // Create a stack to simulate recursion
    stack.push(start_node);

    while (!stack.empty()) {
        unsigned int node = stack.top();
        stack.pop();

        if (visited.find(node) == visited.end()) {
            visited.insert(node);
            component.push_back(node);

            auto it = (*this).al.find(node);
            if (it->second.empty()) continue;
            if (it != (*this).al.end()) {
                for (auto neighbor_it = it->second.rbegin(); neighbor_it != it->second.rend(); ++neighbor_it) {
                    if (visited.find(*neighbor_it) == visited.end()) {
                        stack.push(*neighbor_it);
                    }
                }
            }
        }
    }
}

std::vector<unsigned int> DeBruijnGraph::traceEulerianPath(unsigned int startNode) {
    std::vector<unsigned int> path;
    std::stack<unsigned int> stack;
    stack.push(startNode);

    while (!stack.empty() && (*this).al[startNode].size() > 0) {
        unsigned int currNode = stack.top();

        // If the current node has unused edges
        if ((*this).al[currNode].size() > 0) {
            // Traverse the next edge
            unsigned int nextNode = (*this).al[currNode].back();
            (*this).al[currNode].pop_back(); // Remove the edge to prevent reuse
            stack.push(nextNode);
        } else {
            // If no unused edges are left, backtrack
            path.push_back(currNode);
            stack.pop();
        }
    }

    return path;
}
