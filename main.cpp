#include <vector>
#include <iostream> 

using namespace std;

class Matrix {
public:
    Matrix(int n = 0) : n_(n), field_(n, vector<int>(n)) {

    }

    Matrix(const Matrix &other) {
        n_ = other.n_;
        field_ = vector<vector<int> >(n_, vector<int>(n_));
        for (int i = 0; i < n_; i++) {
            for (int j = 0; j < n_; j++) {
                field_[i][j] = other.field_[i][j];
            }
        }
    }

    Matrix& operator = (const Matrix &other) {
        n_ = other.n_;
        field_ = vector<vector<int> >(n_, vector<int>(n_));
        for (int i = 0; i < n_; i++) {
            for (int j = 0; j < n_; j++) {
                field_[i][j] = other.field_[i][j];
            }
        }
        return *this;
    }

    int get_determinant() {
        return find_determinant(field_, n_);
    }

    int get_size() {
        return n_;
    }

    Matrix apply_permutations(vector<int> &white_permutation, 
                              vector<int> &black_permutation) {
        vector<int> white_numbers;
        vector<int> black_numbers;
        int n = n_;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + j) % 2 == 0) {
                    white_numbers.push_back(field_[i][j]);
                } else {
                    black_numbers.push_back(field_[i][j]);
                }
            }
        }
        Matrix new_matrix(n);
        int white_used = 0;
        int black_used = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + j) % 2 == 0) {
                    new_matrix.field_[i][j] = white_numbers[white_permutation[white_used++]];
                } else {
                    new_matrix.field_[i][j] = black_numbers[black_permutation[black_used++]];
                }
            }
        }
        return new_matrix;
    }

    friend istream& operator >> (istream &in, Matrix &m);
    friend ostream& operator << (ostream &out, Matrix m);

private:
    int find_determinant(vector<vector<int> > &field, int n) {
        if (n == 1)
            return field[0][0];
        else if (n == 2)
            return field[0][0] * field[1][1] - field[0][1] * field[1][0];
        else {
            int det = 0;
            for (int k = 0; k < n; k++) {
                vector<vector<int> > adjunct(n - 1, vector<int>(n - 1));
                for (int i = 1; i < n; i++) {
                    int t = 0;
                    for (int j = 0; j < n; j++) {
                        if (j == k)
                            continue;
                        adjunct[i - 1][t] = field[i][j];
                        t++;
                    }
                }
                det += (k % 2 == 0 ? 1 : -1) * field[0][k] * find_determinant(adjunct, n - 1);
            }
            return det;
        }
    }

    int n_;
    vector<vector<int> > field_;
};

istream& operator >> (istream &in, Matrix &m) {
    in >> m.n_;
    int n = m.n_;
    m.field_ = vector<vector<int> > (n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> m.field_[i][j];
        }
    }
}

ostream& operator << (ostream &out, Matrix m) {
    int n = m.n_;
    out << n << "\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            out << m.field_[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void my_swap(vector<int> &a, int i, int j) {
    int s = a[i];
    a[i] = a[j];
    a[j] = s;
}

bool next_perm(vector<int> &a, int n) {
    int j = n - 2;
    while (j != -1 && a[j] >= a[j + 1]) {
        j--;
    }
    if (j == -1)
        return false;
    int k = n - 1;
    while (a[j] >= a[k]) {
        k--;
    }
    my_swap(a, j, k);
    int l = j + 1, r = n - 1;
    while (l < r) {
        my_swap(a, l++, r--);
    }
    return true;
}

void initialize_perm(vector<int> &permutation, int n) {
    for (int i = 0; i < n; i++) {
        permutation[i] = i;
    }
}

int main() {
    Matrix matrix;
    cin >> matrix;
    int n = matrix.get_size();
    int n_white = (n * n + 1) / 2;
    int n_black = (n * n) / 2;

    vector<int> white_permutation(n_white);
    vector<int> black_permutation(n_black);
    initialize_perm(white_permutation, n_white);
    initialize_perm(black_permutation, n_black);
    Matrix best_matrix = matrix;
    int best_det = matrix.get_determinant();

    while (true) {
        if (!next_perm(white_permutation, n_white)) {
            initialize_perm(white_permutation, n_white);
            if (!next_perm(black_permutation, n_black)) {
                break;
            }
        }
        Matrix another = matrix.apply_permutations(white_permutation, black_permutation);
        if (another.get_determinant() > best_det) {
            best_matrix = another;
            best_det = another.get_determinant();
        }
    }
    cout << "Best determinant : " << best_det << "\n" << best_matrix << "\n";
}
