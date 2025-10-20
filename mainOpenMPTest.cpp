#include <bits/stdc++.h>
using namespace std;

struct Position {
    int row, col;
    Position(int r = -1, int c = -1) : row(r), col(c) {}
    bool operator==(const Position& o) const {
        return row == o.row && col == o.col;
    }
    string toString() const {
        return "(" + to_string(row) + "," + to_string(col) + ")";
    }
    friend ostream& operator<<(ostream& os, const Position& p) {
        os << p.toString();
        return os;
    }
};

class GraphT {
public:
    vector<vector<char>> graph;
    vector<vector<double>> tdiValues;
    vector<vector<Position>> index_data;
    int size;

    GraphT(int size,vector<vector<char>> graph) {
        generateGraphT(size,graph);
        restartTdiValues();
        this->size = size;
    }

    void generateGraphT(int size,vector<vector<char>> graph) {
        this->graph = graph;
        this->tdiValues = vector<vector<double>>(size, vector<double>(size, 0.0));
        this->index_data = vector<vector<Position>>(size, vector<Position>(size, Position(-1, -1)));
        this->size = size;
    }

    void restartTdiValues() {
        double initialTdiValue = 10 * GraphT::distance(Position(0,0), Position((int)tdiValues.size()-1, (int)tdiValues[0].size()-1));
        //highly parallelizable
        for (size_t i = 0; i < tdiValues.size(); ++i) {
            for (size_t j = 0; j < tdiValues[0].size(); ++j) {
                tdiValues[i][j] = initialTdiValue;
                index_data[i][j] = Position(-1, -1);
            }
        }
        if (!tdiValues.empty()){
            tdiValues[size-1][size-1] = 0;
        }
    }

    static double distance(const Position& p1, const Position& p2) {
        double dr = (double)(p1.row - p2.row);
        double dc = (double)(p1.col - p2.col);
        return sqrt(dr*dr + dc*dc);
    }

    //this is not parllelizable as it contains too much interaction with shared memory but can be runned with multiple ants at a time
    void updateTdiValues(const vector<Position>& path) {
        if (path.size() < 2) return;
        for (int idx = (int)path.size() - 2; idx > -1; --idx) {
            const Position &pcur = path[idx];
            const Position &pnext = path[idx+1];
            double newVal = GraphT::distance(pcur, pnext) + this->tdiValues[pnext.row][pnext.col];
            if (newVal < this->tdiValues[pcur.row][pcur.col]) {
                this->tdiValues[pcur.row][pcur.col] = newVal;
                this->index_data[pcur.row][pcur.col].row = pnext.row;
                this->index_data[pcur.row][pcur.col].col = pnext.col;
            }
        }
        return;
    }
};

class AntT {
public:
    vector<Position> path;
    double pathCost;
    GraphT* grid;
    int stepSize;
    double alpha, beta;

    mt19937 rng;

    AntT(GraphT* grid, int stepSize, double alpha, double beta, unsigned int seed) {
        this->grid = grid;
        this->path.clear();
        this->pathCost = 0.0;
        this->stepSize = stepSize;
        this->alpha = alpha;
        this->beta = beta;
        this->rng = mt19937(seed);
    }

    bool isAllowed(const Position& p1, const Position& p2) {
        // cout << p1 << "  " << p2 << endl;
        double i1 = p1.row + 0.5;
        double j1 = p1.col + 0.5;
        double i2 = p2.row + 0.5;
        double j2 = p2.col + 0.5;

        if (i1 < 0 || i2 < 0 || j1 < 0 || j2 < 0 ||
            i1 >= (double)grid->graph.size() || i2 >= (double)grid->graph.size() ||
            j1 >= (double)grid->graph[0].size() || j2 >= (double)grid->graph[0].size()) {
            return false;
        }
        // if(p2.row==1 && p2.col==0){
        //     cout<< grid->graph[(int)floor(i2)][(int)floor(j2)] << " , " << grid->graph[(int)floor(i2)][(int)floor(j2)];
        // }
        if ((grid->graph[(int)floor(i1)][(int)floor(j1)]==1) || (grid->graph[(int)floor(i2)][(int)floor(j2)]==1)) {
            // if(p2.row==1 && p2.col==0){
            //     cout<< " hii ";
            // }
            return false;
        }

        double dx = i2 - i1;
        double dy = j2 - j1;
        int steps = (int) (max(fabs(dx), fabs(dy)) * 10000.0);
        if (steps <= 0) steps = 1;
        double xinc = dx / (double)steps;
        double yinc = dy / (double)steps;
        //parllelizable just checking at all i1,j1 doing it from k=0 to steps
        for (int k = 0; k < steps; ++k) {
            i1 += xinc;
            j1 += yinc;
            int i, j;
            if (fabs(ceil(i1) - i1) < 1e-9) {
                i = (int)ceil(i1);
            } else {
                i = (int)floor(i1);
            }
            if (fabs(ceil(j1) - j1) < 1e-9) {
                j = (int)ceil(j1);
            } else {
                j = (int)floor(j1);
            }
            if (fabs(i1 - i) < 1e-9 || fabs(j1 - j) < 1e-9) {
                continue;
            }

            if (i >= 0 && i < (int)grid->graph.size() && j >= 0 && j < (int)grid->graph[0].size()) {
                if (grid->graph[i][j]==1) {
                    return false;
                }
            } else {
                return false;
            }
        }

        return true;
    }

    // Return Position(-1,-1) on failure
    Position nextPosition(const Position& currPosition) {
        int row = currPosition.row;
        int col = currPosition.col;
        Position target(grid->size - 1, grid->size - 1);
        vector<Position> allowed;
        #pragma omp parallel for num_threads(2) collapse(2)
        for (int i = row - stepSize; i <= row + stepSize; ++i) {
            for (int j = col - stepSize; j <= col + stepSize; ++j) {
                if (i == row && j == col) continue;
                Position p(i, j);
                if (isAllowed(currPosition, p)) {
                    #pragma omp critical
                    {
                        allowed.push_back(p);
                    }
                }
            }
        }
        if (allowed.empty()) {
            cout << "No solution for this graph" << endl;
            return Position(-1, -1);
        }

        vector<double> probability(allowed.size());
        double sum = 0.0;
        for (size_t idx = 0; idx < allowed.size(); ++idx) {
            int i = allowed[idx].row;
            int j = allowed[idx].col;
            double term1 = pow((1.0 / GraphT::distance(allowed[idx], target)), alpha);
            double term2 = pow((1.0 / grid->tdiValues[i][j]), beta);
            probability[idx] = term1 + term2;
            sum += probability[idx];
        }
        for (size_t i = 0; i < probability.size(); ++i) {
            probability[i] = probability[i] / sum;
        }
        for (size_t i = 1; i < probability.size(); ++i) {
            probability[i] = probability[i-1] + probability[i];
        }

        uniform_real_distribution<double> dist(0.0, 1.0);
        double randv = dist(rng);

        if (randv < probability[0]) {
            return allowed[0];
        }
        for (size_t i = 0; i < probability.size() - 1; ++i) {
            if (randv >= probability[i] && randv < probability[i+1]) {
                return allowed[i+1];
            }
        }
        return allowed.back();
    }

    void findSolution() {
        Position goalPosition(grid->size - 1, grid->size - 1);
        this->path.push_back(Position(0,0));
        for (size_t i = 0; !(this->path[i] == goalPosition); ++i) {
            Position nxt = nextPosition(this->path[i]);
            if (nxt.row == -1 && nxt.col == -1) {
                cout << "nextPosition failed (no allowed). Aborting findSolution for this ant." << endl;
                return;
            }
            this->path.push_back(nxt);
            this->pathCost += GraphT::distance(this->path[this->path.size()-2], this->path[this->path.size()-1]);
        }
    }

    void restartPath() {
        this->path.clear();
        this->pathCost = 0.0;
        return;
    }
};

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
            string gridKeyword;
            int currentIndex, n = 0;
            double p = 0.0;
            long long seed = 0;

            stringstream ss(line);
            ss >> gridKeyword >> currentIndex;

            string token;
            while (ss >> token) {
                if (token.rfind("n=", 0) == 0) n = stoi(token.substr(2));
                else if (token.rfind("p=", 0) == 0) p = stod(token.substr(2));
                else if (token.rfind("seed=", 0) == 0) seed = stoll(token.substr(5));
            }

            if (currentIndex == gridIndex) {
                gridFound = true;
                grid.assign(n, vector<char>());

                for (int i = 0; i < n; ++i) {
                    if (!getline(inputFile, line)) {
                        cerr << "Error: Unexpected EOF while reading GRID " << gridIndex << endl;
                        return {};
                    }
                    stringstream row(line);
                    char cell;
                    while (row >> cell)
                        grid[i].push_back(cell);
                }

                // done reading desired grid
                cout << "Read GRID " << gridIndex
                     << " successfully: n=" << n
                     << " p=" << p
                     << " seed=" << seed << endl;
                break;
            }
        }
    }

    if (!gridFound)
        cerr << "Error: GRID " << gridIndex << " not found in " << filePath << endl;

    return grid;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <path to grids.txt> <grid_index>\n" << std::flush;
        return 1;
    }
    // taking required inputs
    char* filePath = argv[1];
    int gridIndex = atoi(argv[2]);
    // getting grid from file
    vector<vector<char>> test=readGridFromFile(filePath,gridIndex);


    GraphT grid((int)test.size(),test);
    vector<Position> solution;
    double solutionCost = 0.0;
    vector<vector<Position>> solutions;
    int noOfAnts = 10, noOfIterations = 20, stepSize = 3;
    double alpha = 1.5, beta = 0.8;
    vector<double> solutionsCost(noOfIterations);
    vector<unique_ptr<AntT>> ants;
    ants.reserve(noOfAnts);
    Position goalPosition(grid.size - 1, grid.size - 1);

    unsigned int baseSeed = 69;

    //ants creation and pushing into vectors i think with reduction i can do this in O(log(ants count))
    for (int i = 0; i < noOfAnts; ++i) {
        ants.emplace_back(make_unique<AntT>(&grid, stepSize, alpha, beta, baseSeed+i));
    }

    for (int iter = 0; iter < noOfIterations; ++iter) {
        //completely independent so parllelize not even shared memory zero communication needed but after completion barrier is needed
        #pragma omp parallel for num_threads(10)
        for (int j = 0; j < noOfAnts; ++j) {
            ants[j]->findSolution();
        }
        //must not parllelize as it contains too much interaction with shared memory
        for (int j = 0; j < noOfAnts; ++j) {
            grid.updateTdiValues(ants[j]->path);
        }
        for (int j = 0; j < noOfAnts; ++j) {
            ants[j]->restartPath();
        }
        cout << iter << ", "<< std::flush;

        solutions.push_back(vector<Position>());
        solutions.back().push_back(Position(0, 0));
        while (!(solutions.back().back() == goalPosition)) {
            Position last = solutions.back().back();
            Position next = grid.index_data[last.row][last.col];
            solutions.back().push_back(next);
            if (next.row == -1 && next.col == -1) break;
        }

        double sc = 0.0;
        for (size_t k = 1; k < solutions.back().size(); ++k) {
            sc += GraphT::distance(solutions.back()[k-1], solutions.back()[k]);
        }
        solutionsCost[iter] = sc;
    }

    solution.push_back(Position(0, 0));
    while (!(solution.back() == goalPosition)) {
        Position last = solution.back();
        Position next = grid.index_data[last.row][last.col];
        solution.push_back(next);
        if (next.row == -1 && next.col == -1) break;
    }
    //printing the solution do not remove it
    cout << "Solution : " << solution[0]<< std::flush;
    for (size_t j = 1; j < solution.size(); ++j) {
        cout << "," << solution[j]<< std::flush;
        solutionCost += GraphT::distance(solution[j-1], solution[j]);
    }

    cout << "\nSolutionCost : " << solutionCost << endl;

    return 0;
}
