#include <bits/stdc++.h>
#include "json.hpp" // include nlohmann/json header file
using namespace std;
using json = nlohmann::json;

// Convert string in given base -> decimal
long long convertBase(const string &value, int base) {
    long long result = 0;
    for (char c : value) {
        int digit;
        if (isdigit(c)) digit = c - '0';
        else digit = tolower(c) - 'a' + 10;
        result = result * base + digit;
    }
    return result;
}

// Gaussian Elimination
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        // Pivoting
        int maxRow = i;
        for (int k = i+1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i]))
                maxRow = k;
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        // Normalize pivot
        double pivot = A[i][i];
        if (fabs(pivot) < 1e-12) continue;
        for (int j = i; j < n; j++) A[i][j] /= pivot;
        b[i] /= pivot;

        // Eliminate below
        for (int k = i+1; k < n; k++) {
            double factor = A[k][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
    return x;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Read full JSON from stdin
    string input, line;
    while (getline(cin, line)) input += line;
    json j = json::parse(input);

    int n = j["keys"]["n"];
    int k = j["keys"]["k"];
    int m = k - 1; // polynomial degree

    // Decode roots
    vector<long long> roots;
    for (int i = 1; i <= n; i++) {
        if (j.contains(to_string(i))) {
            int base = stoi(j[to_string(i)]["base"].get<string>());
            string val = j[to_string(i)]["value"];
            roots.push_back(convertBase(val, base));
        }
    }

    // Build Vandermonde matrix for first k roots
    vector<vector<double>> A(m, vector<double>(m, 0));
    vector<double> b(m, 0);
    for (int i = 0; i < m; i++) {
        long long x = roots[i];
        double term = 1.0;
        for (int j = 0; j < m; j++) {
            A[i][j] = term;
            term *= x;
        }
        b[i] = -pow((double)x, m);
    }

    // Solve for coefficients c0..c(m-1), fix cm = 1
    vector<double> coeffs = gaussianElimination(A, b);
    coeffs.push_back(1.0);

    // Output
    cout << "Polynomial coefficients (from c0 + c1*x + ... + cm*x^m):\n";
    for (int i = 0; i < coeffs.size(); i++) {
        cout << "c" << i << " = " << coeffs[i] << "\n";
    }

    return 0;
}
