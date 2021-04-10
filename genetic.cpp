#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;

const int POPULATION = 100;
const int GENERATIONS = 10000;
const double RADIATION = 0.05;

double det(const vector<vector<int> >& matrix) {
    int n = matrix.size();
    if (n == 1) {
        return matrix[0][0];
    } else if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    } else {
        int det_ = 0;
        for (int k = 0; k < n; k++) {
            vector<vector<int> > adjunct(n - 1, vector<int>(n - 1));
            for (int i = 1; i < n; i++) {
                int t = 0;
                for (int j = 0; j < n; j++) {
                    if (j == k) {
                        continue;
                    }
                    adjunct[i - 1][t] = matrix[i][j];
                    t++;
                }
            }
            det_ += (k % 2 == 0 ? 1 : -1) * matrix[0][k] * det(adjunct);
        }
        return det_;
    }
}

vector<vector<int> > matrix_from_person(const vector<vector<int> >& person, vector<vector<int> > matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int tmp = matrix[i][j];
            matrix[i][j] = matrix[person[i][j] / n][person[i][j] % n];
            matrix[person[i][j] / n][person[i][j] % n] = tmp;
        }
    }
    return matrix;
}

double fitness(const vector<vector<int> >& person, vector<vector<int> > matrix) {
    return det(matrix_from_person(person, matrix));
}

double rnd(double r) {
    return r * rand() / (1ll << 32);
}

int choose(const vector<double> pref, double r) {
    double point = rnd(r);
    int id = lower_bound(pref.begin(), pref.end(), point) - pref.begin();
    if (id == pref.size()) {
        id--;
    }
    return id;
}

vector<vector<int> > love(const vector<vector<int> >& mom, const vector<vector<int> >& dad) {
    int n = mom.size();
    vector<vector<int> > baby(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            baby[i][j] = mom[i][j];
            if (rand() & 1) {
                baby[i][j] = dad[i][j];
            }
            if (rnd(1) < RADIATION) {
                do {
                    baby[i][j] = rand() % (n * n - i * n - j) + i * n + j;
                } while ((baby[i][j] / n + baby[i][j] % n + i + j) & 1);
            }
        }
    }
    return baby;
}

vector<vector<int> > solve(vector<vector<int> > in) {
    int n = in.size();
    vector<vector<vector<int> > > population(POPULATION, vector<vector<int> >(n, vector<int>(n, n * n - 1)));
    for (int i = 0; i < POPULATION; i++) {
        for (int j = 0; j < n; j++) {
            for (int l = 0; l < n; l++) {
                do {
                    population[i][j][l] = rand() % (n * n - j * n - l) + j * n + l;
                } while ((population[i][j][l] / n + population[i][j][l] % n + j + l) & 1);
            }
        }
    }
    vector<vector<int> > the_fittest = population[0];
    double max_fit = fitness(the_fittest, in);
    
    for (int g = 0; g < GENERATIONS; g++) {
        vector<vector<vector<int> > > new_population(POPULATION);
        vector<double> fits(POPULATION);
        double r = 0;
        for (int i = 0; i < POPULATION; i++) {
            r += (fits[i] = fitness(population[i], in));
            if (fits[i] > max_fit) {
                max_fit = fits[i];
                the_fittest = population[i];
            }
        }
        vector<double> fits_pref;
        double sum = 0;
        for (int i = 0; i < POPULATION; i++) {
            sum += fits[i];
            fits_pref.push_back(sum);
        }
        
        for (int i = 0; i < POPULATION; i++) {
            int mom = choose(fits_pref, r);
            int dad = choose(fits_pref, r);
            new_population[i] = love(population[mom], population[dad]);
        }

        swap(new_population, population);
    }

    return matrix_from_person(the_fittest, in);
}

int main() {
    srand(time(NULL));
    int n;
    cin >> n;
    vector<vector<int> > in(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> in[i][j];
        }
    }
    vector<vector<int> > out = solve(in);
    cout << "Best determinant : " << det(out) << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << out[i][j] << " ";
        }
        cout << endl;
    }
}