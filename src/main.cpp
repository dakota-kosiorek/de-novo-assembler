#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <limits.h>
#include <de_bruijn.hpp>

// Multithread points
//      splitting reads into kmers?

std::string elapsedTime(const std::chrono::steady_clock::time_point&);
bool onlyValidNucleotides(const std::string&);

int main(int argc, char* argv[]) {
    std::ostringstream errorMsg;
    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

    // TODO: Get arguments

    // FOR USER: put your filename here
    char fileName[] = "YOUR FILE HERE";

    try {
        std::cout << elapsedTime(startTime) << " Opening file '" << fileName << "'..." << std::endl; 
        std::ifstream fastqFile(fileName);

        // See if the file was opened successfully
        if (fastqFile.is_open() == 0) {
            errorMsg  << "Could not open file '" << fileName << "'";
            throw std::runtime_error(errorMsg.str());
        }

        // Make sure the fastq file is structured correctly
        std::cout << elapsedTime(startTime) << " Checking file integrity and loading reads..." << std::endl;

        int fastqStatus = 0; // 0: Good, 1: bad label, 2: bad seq, 3: bad seperator, 4: bad qscores 
        unsigned int lineCount = 0;

        std::string line = "";
        std::string label = "";
        std::string seq = "";
        std::string qScores = "";
        unsigned int shortestReadLen = UINT_MAX;
        unsigned int longestReadLen = 0;
        unsigned int totalReadLen = 0;

        std::vector<read> reads;

        while (std::getline(fastqFile, line)) {
            if (lineCount % 4 == 0) {               // Label
                label = line;
                if (label[0] != '@') {
                    fastqStatus = 1;
                    break;
                }

            } else if ((lineCount - 1) % 4 == 0) {  // Sequence
                seq = line;
                if (!onlyValidNucleotides(seq)) {
                    fastqStatus = 2;
                    //break;
                } else {
                    unsigned int seqLen = seq.length();

                    if (seqLen < shortestReadLen) 
                        shortestReadLen = seq.length();

                    if (seqLen > longestReadLen)
                        longestReadLen = seq.length();

                    totalReadLen += seqLen;
                }
            } else if ((lineCount - 2) % 4 == 0) {  // Seperator
                if (line[0] != '+' || line.length() != 1) {
                    fastqStatus = 3;
                    break;
                }

            } else {                        // Q scores
                qScores = line;
                if (qScores.length() != seq.length()) {
                    fastqStatus = 4;
                    break;
                }

                //reads.push_back({label, seq, qScores});
                if (fastqStatus == 0) reads.push_back({seq});
                if (fastqStatus == 2) fastqStatus = 0; 
            }

            lineCount++;
        }

        fastqFile.close();

        if (fastqStatus != 0 || lineCount % 4 != 0) {
            std::string lineError;

            if (fastqStatus != 0) {
                switch (fastqStatus) {
                    case 1: lineError = "Label does not start with @"; break;
                    case 2: lineError = "Incorrect nucleotides found in sequence (not in {A, T, C, G})"; break;
                    case 3: lineError = "Seperator not '+'"; break; 
                    case 4: lineError = "Length of Q scores != length of sequence"; break;
                    default: lineError = ""; break;
                }

                errorMsg << lineError << " (line " << lineCount << ")"; 
            } else {
                errorMsg << "all reads do not have all fastq 4 elements (label, seq, seperator, Q scores)";
            }
            
            throw std::runtime_error(errorMsg.str());
        }

        reads.shrink_to_fit();
        std::cout << elapsedTime(startTime) << " FASTQ file is correctly formatted\n";
        std::cout << "  --> " << reads.size() << " reads\n";
        std::cout << "  --> Total read length: " << totalReadLen << " nt\n";
        std::cout << "  --> Shortest read length: " << shortestReadLen << " nt\n";
        std::cout << "  --> Longest read length: " << longestReadLen << " nt" << std::endl;

        // Make De Bruijn graph of reads (k-mer <= shortestReadLen)
        unsigned int k = 21;
        std::cout << elapsedTime(startTime) << " Assembling De Bruijn graph (k = " << k << ")\n";
        std::cout << "  --> Total (k-1)-mers to generate: " << totalReadLen - reads.size() * (k - 1) << std::endl;

        DeBruijnGraph dbg(reads, k, startTime);    // Generates the De Bruijn graph

        std::cout << elapsedTime(startTime) << " Finished assembling De Bruijn graph (k = " << k << ")" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << elapsedTime(startTime) << " ERROR: " << e.what() << std::endl;
    }

    return 0;
}

// Check to see if string is only made of valid nucleotides
bool onlyValidNucleotides(const std::string& s) {
    bool isValid = true;
    
    for (unsigned int i = 0; (i < s.length() && isValid); i++) {
        char c = toupper(s[i]);
        if (!(c == 'A' || c == 'T' || c == 'G' || c == 'C')) {
            isValid = false;
        }
    }

    return isValid;
}