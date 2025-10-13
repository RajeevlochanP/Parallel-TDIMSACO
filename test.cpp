#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept> // For std::invalid_argument in stoi

// Using specific declarations is often safer than 'using namespace std;'
using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;

// The function you provided is already well-written. No changes are needed here.
vector<vector<char>> readGridFromFile(const string& filePath, int gridIndex) {
    vector<vector<char>> grid;
    
    ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open file at path: " << filePath << endl;
        return grid;
    }

    string line;
    bool gridFound = false;

    while (getline(inputFile, line)) {
        if (line.rfind("GRID ", 0) == 0) {
            stringstream headerStream(line);
            string gridKeyword, nKeyword;
            int currentIndex, n;
            char equalsChar;

            headerStream >> gridKeyword >> currentIndex >> nKeyword >> equalsChar >> n;

            if (currentIndex == gridIndex) {
                gridFound = true;
                
                grid.resize(n);

                for (int i = 0; i < n; ++i) {
                    if (getline(inputFile, line)) {
                        stringstream rowStream(line);
                        char cellValue;
                        while (rowStream >> cellValue) {
                            grid[i].push_back(cellValue);
                        }
                    } else {
                        cerr << "Error: File format is incorrect. Reached end of file while reading grid data." << endl;
                        return {};
                    }
                }
                break;
            }
        }
    }

    inputFile.close();
    
    if (!gridFound) {
        cerr << "Error: GRID " << gridIndex << " was not found in the file." << endl;
    }

    return grid;
}


int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <path_to_file> <grid_index>\n";
        return 1;
    }

    // --- IMPROVEMENT: Use std::string and std::stoi for safer parsing ---
    string filePath = argv[1];
    int gridIndex;

    try {
        gridIndex = std::stoi(argv[2]); // Converts string to int, throws exception on failure
    } catch (const std::invalid_argument& e) {
        cerr << "Error: Invalid grid index. Please provide an integer." << endl;
        return 1;
    }

    // Getting grid from file
    vector<vector<char>> test = readGridFromFile(filePath, gridIndex);
    
    // --- IMPROVEMENT: Use range-based for loops for safer and cleaner iteration ---
    // This also avoids the potential crash from test[0] if 'test' is empty.
    if (!test.empty()) {
        for (const auto& row : test) {
            for (const auto& cell : row) {
                cout << cell << " ";
            }
            cout << endl;
        }
    }
    
    return 0;
}