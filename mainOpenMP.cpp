#include <bits/stdc++.h>
#include <omp.h>
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
};

class GraphT {
public:
    vector<vector<char>> graph;
    vector<vector<double>> tdiValues;
    vector<vector<Position>> index_data;
    int size;
    // void merge(const std::vector<GraphT>& localGrids);

    GraphT(int size) {
        generateGraphT(size);
        restartTdiValues();
        this->size = size;
    }

    void generateGraphT(int size) {
        this->graph = vector<vector<char>>(size, vector<char>(size, 0));
        this->tdiValues = vector<vector<double>>(size, vector<double>(size, 0.0));
        this->index_data = vector<vector<Position>>(size, vector<Position>(size, Position(-1, -1)));
        this->size = size;

        vector<pair<int,int>> trueIndices = {
            { 1, 4 }, { 1, 5 }, { 2, 4 }, { 2, 4 }, { 2, 5 }, { 3, 4 }, { 3, 5 }, { 1, 11 }, { 2, 11 }, { 3, 11 },
            { 1, 13 }, { 1, 14 }, { 2, 13 }, { 2, 14 }, { 3, 13 }, { 3, 14 }, { 2, 17 }, { 2, 18 }, { 6, 3 },
            { 6, 4 }, { 7, 3 }, { 7, 4 }, { 8, 3 }, { 8, 4 }, { 6, 8 }, { 6, 9 }, { 7, 8 }, { 7, 9 }, { 8, 8 },
            { 8, 9 }, { 6, 15 }, { 6, 16 }, { 6, 17 }, { 7, 15 }, { 7, 16 }, { 7, 17 }, { 8, 15 }, { 8, 16 },
            { 8, 17 }, { 10, 3 }, { 10, 4 }, { 10, 5 }, { 10, 9 }, { 11, 7 }, { 11, 8 }, { 11, 9 }, { 11, 10 },
            { 11, 11 }, { 11, 12 }, { 11, 13 }, { 12, 7 }, { 12, 8 }, { 12, 9 }, { 12, 10 }, { 12, 11 }, { 12, 12 },
            { 12, 13 }, { 13, 7 }, { 13, 8 }, { 13, 9 }, { 13, 10 }, { 13, 11 }, { 13, 12 }, { 13, 13 }, { 12, 5 },
            { 13, 5 }, { 14, 5 }, { 15, 5 }, { 16, 5 }, { 14, 2 }, { 14, 3 }, { 15, 2 }, { 15, 3 }, { 17, 2 },
            { 17, 3 }, { 18, 2 }, { 18, 3 }, { 15, 7 }, { 16, 7 }, { 17, 7 }, { 18, 7 }, { 15, 12 }, { 15, 13 },
            { 15, 14 }, { 16, 12 }, { 16, 13 }, { 16, 14 }, { 17, 12 }, { 17, 13 }, { 17, 14 }, { 18, 12 },
            { 18, 13 }, { 18, 14 }, { 17, 16 }, { 17, 17 }, { 18, 16 }, { 18, 17 }
        };

        for (auto &pr : trueIndices) {
            int r = pr.first;
            int c = pr.second;
            if (r >= 0 && r < size && c >= 0 && c < size) {
                this->graph[r][c] = 1;
            }
        }
    }

    void restartTdiValues() {
        double initialTdiValue = 10 * GraphT::distance(Position(0,0), Position((int)tdiValues.size()-1, (int)tdiValues[0].size()-1));
        //highly parallelizable
        #pragma omp parallel for collapse(2) num_threads(20)
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
    void merge(const std::vector<GraphT>& localGrids) {
        if (localGrids.empty()) return;
        int rows = this->tdiValues.size();
        int cols = this->tdiValues[0].size();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double bestVal = this->tdiValues[i][j];
                Position bestPos = this->index_data[i][j];
                for (const auto& g : localGrids) {
                    double val = g.tdiValues[i][j];
                    if (val < bestVal) {
                        bestVal = val;
                        bestPos = g.index_data[i][j];
                    }
                }
                this->tdiValues[i][j] = bestVal;
                this->index_data[i][j] = bestPos;
            }
        }
    }
};


class AntT {
public:
    vector<Position> path;
    double pathCost;
    GraphT* grid;
    int stepSize;
    double alpha, beta;

    static mt19937 rng;

    AntT(GraphT* grid, int stepSize, double alpha, double beta) {
        this->grid = grid;
        this->path.clear();
        this->pathCost = 0.0;
        this->stepSize = stepSize;
        this->alpha = alpha;
        this->beta = beta;
    }

    bool isAllowed(const Position& p1, const Position& p2) {
        double i1 = p1.row + 0.5;
        double j1 = p1.col + 0.5;
        double i2 = p2.row + 0.5;
        double j2 = p2.col + 0.5;

        if (i1 < 0 || i2 < 0 || j1 < 0 || j2 < 0 ||
            i1 >= (double)grid->graph.size() || i2 >= (double)grid->graph.size() ||
            j1 >= (double)grid->graph[0].size() || j2 >= (double)grid->graph[0].size()) {
            return false;
        }

        if (grid->graph[(int)floor(i1)][(int)floor(j1)] ||
            grid->graph[(int)floor(i2)][(int)floor(j2)]) {
            return false;
        }

        double dx = i2 - i1;
        double dy = j2 - j1;
        int steps = (int) (max(fabs(dx), fabs(dy)) * 1000.0);
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
        #pragma omp parallel num_threads(20)
        {
            vector<Position> local_allowed;
        
            #pragma omp for nowait collapse(2)
            for (int i = row - stepSize; i <= row + stepSize; ++i) {
                for (int j = col - stepSize; j <= col + stepSize; ++j) {
                    if (i == row && j == col) continue;
                    Position p(i, j);
                    if (isAllowed(currPosition, p)) {
                        local_allowed.push_back(p);
                    }
                }
            }
        
            #pragma omp critical
            allowed.insert(allowed.end(), local_allowed.begin(), local_allowed.end());
        }

        if (allowed.empty()) {
            cerr << "No solution for this graph" << endl;
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
                cerr << "nextPosition failed (no allowed). Aborting findSolution for this ant." << endl;
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

mt19937 AntT::rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count()); //for generating high quality random numbers

int main() {
    GraphT grid(20);
    vector<Position> solution;
    double solutionCost = 0.0;
    vector<vector<Position>> solutions;
    int noOfAnts = 50, noOfIterations = 100, stepSize = 3;
    double alpha = 1.5, beta = 0.8;
    vector<double> solutionsCost(noOfIterations);
    vector<unique_ptr<AntT>> ants;
    ants.reserve(noOfAnts);
    Position goalPosition(grid.size - 1, grid.size - 1);

    //ants creation and pushing into vectors i think with reduction i can do this in O(log(ants count))
    for (int i = 0; i < noOfAnts; ++i) {
        ants.emplace_back(make_unique<AntT>(&grid, stepSize, alpha, beta));
    }

    for (int iter = 0; iter < noOfIterations; ++iter) {
        //completely independent so parllelize not even shared memory zero communication needed but after completion barrier is needed
        #pragma omp parallel for schedule(dynamic) num_threads(20)
        for (int j = 0; j < noOfAnts; ++j) {
        ants[j]->findSolution();
        }
        // //must not parllelize as it contains too much interaction with shared memory
        // for (int j = 0; j < noOfAnts; ++j) {
        //     grid.updateTdiValues(ants[j]->path);
        // }
        std::vector<GraphT> localGrids(omp_get_max_threads(), grid);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            #pragma omp for schedule(dynamic)
            for (int j = 0; j < noOfAnts; ++j){
                localGrids[tid].updateTdiValues(ants[j]->path);
            }
        }

        //merge results
        grid.merge(localGrids);
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
    for (size_t j = 1; j < solution.size(); ++j) {
        solutionCost += GraphT::distance(solution[j-1], solution[j]);
    }

    cout << "\nTdi value of (0,0) :" << grid.tdiValues[0][0] << endl;
    return 0;
}
