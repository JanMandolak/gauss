#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

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
    cout << endl;
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

vector<double> provide_the_only_solution(vector< vector<double> > &Ab) {
    int n = Ab.size();

    vector<double> solution(n);
    for (int i=n-1; i>=0; i--) {
        solution[i] = Ab[i][n] / Ab[i][i];
        for (int k=i-1;k>=0; k--) {
            Ab[k][n] -= Ab[k][i] * solution[i];
        }
    }
    return solution;
}


int frobenius(const vector<vector<double>> Ab) {
    int n = Ab.size();
    int row = n-1;
    int dim = 0;

    bool null_row;
    while (row > -1) {
        null_row = true;
        for (int col=row; col < n; col++) {
            if (Ab[row][col] != 0) {
                null_row = false;
                break;
            }
        }
        if (null_row) {
            dim++;
        } else {
            break;}
        row--;
    }
    for (row=n-1; row > n-1-dim; row--) {
        if (Ab[row][n] != 0) {
            dim = -1;
        }
    }

    return dim;
}

vector<vector<double>> get_upper(vector<vector<double>> Ab) {
    int n = Ab.size();
    int cur_row = 0;

    for (int cur_col=0; cur_col < n; cur_col++) {
        int max_row = find_max_under(Ab, cur_col, cur_col);
        if (Ab[max_row][cur_col] == 0 ) {
            continue;}
        if (max_row != cur_row){
            swap_rows(Ab, cur_row, max_row);
        }
        null_under(Ab, cur_row, cur_col);
        cur_row++;
    }

    return Ab;
}
//
//vector<int> find_loose(vector<vector<double>> upper, int dim) {
//    vector<int> loose(dim);
//    int n = upper.size();
//
//
//    }
//    return vector<int>();
//}



int main() {
    int n;
    cin >> n;

    vector<double> line(n+1,0);
    vector<vector<double>> Ab(n, line);

    // Read input data
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            cin >> Ab[i][j];
        }
    }

    for (int i=0; i<n; i++) {
        cin >> Ab[i][n];
    }

    print(Ab);

    std::tuple<int, double> a;

    vector<vector<double>> upper = get_upper(Ab);

    int dimension = frobenius(upper);

    if (dimension == -1) {
        cout << "There is no solution. " ;
        return 0;
    }
    if (dimension == 0) {

        vector<double> x(n);
        x = provide_the_only_solution(upper);

        cout << "Result:\t";
        for (int i=0; i<n; i++) {
            cout << x[i] << " ";
        }
        cout << endl;
        return 0;
    }

    if (dimension > 0) {


        vector<double> x(n);
        vector<int> loose_indices(dimension);

//        loose_indices = find_loose(upper, dimension);
        x = provide_the_only_solution(upper);

        cout << "Result:\t";
        for (int i=0; i<n; i++) {
            cout << x[i] << " ";
        }
        cout << endl;
        return 0;
    }


}