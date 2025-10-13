#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

int main() {
    std::ofstream file("input.csv"); // open file for writing
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return 1;
    }

    srand(static_cast<unsigned>(time(nullptr))); // seed random generator

    int n = 500; // number of lines
    for (int i = 0; i < n; ++i) {
        int num = 1+ (rand() % 100);  // random integer between 0 and 100
        double prob = static_cast<double>(rand()) / RAND_MAX * 0.5; // random double in [0, 0.5)
        
        file << num << "," << std::fixed << std::setprecision(2) << prob << "\n";
    }

    file.close();
    std::cout << "✅ Data written successfully to input.csv\n";
    return 0;
}