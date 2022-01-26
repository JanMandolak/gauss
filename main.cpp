#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

void print_help() {
    std::cout << "Usage: Gauss [OPTION]...  [FILE]..." << std::endl
              << "Options: " << std::endl <<
              "-h, --help               to print help manual" << std::endl <<
              "-o, --one-thread         use one thread only" << std::endl <<
              "-m, --multiple-thread    use multiple threads, if available" << std::endl;
}

void print(vector< vector<double> > Ab) {
    int n = Ab.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            cout << Ab[i][j] << "\t";
            if (j == n-1) {
                cout << "| ";
            }
        }
        cout << "\n";
    }
//    cout << endl;
}

void swap_rows(vector< vector<double> > &Ab, const int row1, const int row2) {
    int n = Ab.size();

    for (int k=0; k<n+1;k++) {
        double tmp = Ab[row2][k];
        Ab[row2][k] = Ab[row1][k];
        Ab[row1][k] = tmp;
    }
}

int find_max_under(const vector< vector<double> > Ab, const int row, const int col) {
    int n = Ab.size();

    double max_value = abs(Ab[row][col]);
        int max_row = row;
        for (int k=row+1; k<n; k++) {
            if (abs(Ab[k][col]) > max_value) {
                max_value = abs(Ab[k][col]);
                max_row = k;
            }
        }
        return max_row;
}

void null_under(vector< vector<double> > &Ab, const int upper_row, const int left_col) {
    int n = Ab.size();

    for (int k=upper_row+1; k<n; k++) {
        double multiplicator = -Ab[k][left_col] / Ab[upper_row][left_col];
        for (int j=upper_row; j<n+1; j++) {
            if (upper_row==j) {
                Ab[k][j] = 0;
            } else {
                Ab[k][j] += multiplicator * Ab[upper_row][j];
            }
        }
    }
}

vector<double> provide_one_solution(const vector< vector<double> > & upper, int dim, vector<int> pivots) {
    int n = upper.size();
    int m = n-dim;
    vector<double> line(m+1,0);
    vector<vector<double>> temp(m, line);

    for (int col=0; col<m; col++) {
        for (int row=0; row<m; row++) {
            temp[row][col] = upper[row][pivots[col]];
        }
    }
    for (int row=0; row<m; row++) {
        temp[row][m] = upper[row][n];
    }


    vector<double> partial_solution(m);
    for (int i=m-1; i>=0; i--) {
        partial_solution[i] = temp[i][m] / temp[i][i];
        for (int k=i-1;k>=0; k--) {
            temp[k][m] -= temp[k][i] * partial_solution[i];
        }
    }
    vector<double> solution(n,0);
    for (int i=0; i< m; i++) { solution[pivots[i]] = partial_solution[i];}

    return solution;
}


bool frobenius(const vector<vector<double>> Ab, int dim) {
    int n = Ab.size();

    for (int row=n-1; row > n-1-dim; row--) {
        if (Ab[row][n] != 0) {
            return false;
        }
    }

    return true;
}

tuple<vector<vector<double>>, vector<int>> get_upper(vector<vector<double>> Ab) {
    int n = Ab.size();
    int cur_row = 0;
    vector<int> pivots;

    for (int cur_col=0; cur_col < n; cur_col++) {
        int max_row = find_max_under(Ab, cur_row, cur_col);
        if (Ab[max_row][cur_col] == 0 ) {
            continue;}
        if (max_row != cur_row){
            swap_rows(Ab, cur_row, max_row);
        }
        null_under(Ab, cur_row, cur_col);
        pivots.push_back(cur_col);
        cur_row++;
    }


    return make_tuple(Ab, pivots);
}


vector<vector<double>> find_base(vector<vector<double>> upper, int dim, vector<int> pivots) {
    int n = upper.size();
    vector<double> line(n);
    vector<vector<double>> base;
    vector<vector<double>> temp;
    int col = 0;
    for (int row=0; row<dim; row++) {
        while (binary_search(pivots.begin(), pivots.end(), col)) {col++;}
        temp = upper;
        for (int row_=0; row_<n- dim; row_++) {temp[row_][n] = -upper[row_][col];}
        base.push_back(provide_one_solution(temp, dim,pivots));
        base[row][col] = 1;
        col++;
    }
    return base;
}

vector<tuple<int, vector<vector<double>>>> read_matrices(const string &path) {
        vector<tuple<int, vector<vector<double>>>> problems;
        fstream infile;
        infile.open(path);
        if (infile.is_open()) {
                int n;
                while (true) {
                    if (!(infile >> n)) break;
                    vector<double> line(n + 1, 0);
                    vector<vector<double>> Ab(n, line);

                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n + 1; j++) {
                            infile >> Ab[i][j];
                        }
                    }
                    problems.emplace_back(n, Ab);
                }
            return problems;
        } else {
            std::cout << "File " << path << " NOT found!" << std::endl;
            return problems;
        }
        infile.close();
}

bool in_options(std::vector<std::string> &args, const std::string &str) {
    return std::find(args.begin(), args.end(), str) != args.end();
}


int main(int argc, char *argv[]) {
//    if (argc < 2) {
    if (false) {
        print_help();
        return 0;}
    else {
        std::vector<std::string> arguments(argv + 1, argv + argc);
        if (in_options(arguments, "--help") || in_options(arguments, "-h")) {
            print_help();
//            add readme
            return 0;
        }

//        if (in_options(arguments, "--one-thread") || in_options(arguments, "-o")) {
        if (true) {
            std::string path;
            path = "multi.txt";
//            std::cin >> path;
            vector<tuple<int, vector<vector<double>>>> problems = read_matrices(path);

            while (!problems.empty())
            {
                int n;
                vector<vector<double>> Ab;
                tie(n, Ab) = problems.back();
                problems.pop_back();

                print(Ab);

                std::tuple<int, double> a;
                vector<vector<double>> upper;
                vector<int> pivots;
                tie(upper, pivots) = get_upper(Ab);

                int dimension;
                dimension = n - pivots.size();


                if (!frobenius(upper, dimension)) {
                    cout << "There is no solution. \n\n";
                }
                else {
                    vector<double> x(n);
                    vector<vector<double>> base;
                    x = provide_one_solution(upper, dimension, pivots);

                    if (dimension > 0) {
                        base = find_base(upper, dimension, pivots);
                    }

                    cout << "Result:\t";
                    cout << "K={[";
                    for (int i = 0; i < n; i++) {
                        if (i < n - 1) { cout << x[i] << " "; }
                        else { cout << x[i] << "]"; }
                    }
                    if (dimension == 0) { cout << "}"; }
                    else {
                        cout << " + <";
                        for (int s = 0; s < dimension; s++) {
                            cout << "[";
                            for (int i = 0; i < n; i++) {
                                if (i < n - 1) { cout << base[s][i] << " "; }
                                else { cout << base[s][i] << "]"; }
                            }
                            if (s != dimension - 1) { cout << ", "; }
                        }
                        cout << ">";
                    }
                    cout << "}\n" << endl;
                }


            }
            return 0;
        }
    }
}