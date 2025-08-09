#include <bits/stdc++.h>
#include <nlohmann/json.hpp>
#include <fstream>
using namespace std;
using json = nlohmann::json;

// Function to convert a number from given base to decimal
long long baseToDecimal(const string &numStr, int base) {
    long long result = 0;
    for (char c : numStr) {
        int digit;
        if (isdigit(c)) digit = c - '0';
        else digit = toupper(c) - 'A' + 10;
        result = result * base + digit;
    }
    return result;
}

// Gaussian elimination for solving system of linear equations
vector<double> gaussianSolve(vector<vector<double>> a, vector<double> b) {
    int n = a.size();
    for (int i = 0; i < n; i++) {
        // Pivot
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
            if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                maxRow = k;
        swap(a[i], a[maxRow]);
        swap(b[i], b[maxRow]);

        // Eliminate
        for (int k = i + 1; k < n; k++) {
            double factor = a[k][i] / a[i][i];
            for (int j = i; j < n; j++)
                a[k][j] -= factor * a[i][j];
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }
    return x;
}

int main() {
    // Read JSON file
    ifstream inFile("input.json");
    json j;
    inFile >> j;

    int n = j["keys"]["n"];
    int k = j["keys"]["k"];

    // Collect given points
    vector<pair<long long, long long>> points;
    for (auto &el : j.items()) {
        if (el.key() == "keys") continue;
        long long x = stoll(el.key());
        long long y = baseToDecimal(el.value()["value"], stoi(el.value()["base"]));
        points.push_back({x, y});
    }

    // Sort points by x
    sort(points.begin(), points.end());

    // We need first (n) points to solve
    if ((int)points.size() < n) {
        cout << "Not enough points to solve polynomial\n";
        return 0;
    }

    // Build Vandermonde matrix
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);
    for (int i = 0; i < n; i++) {
        double xi = points[i].first;
        double power = 1;
        for (int j = 0; j < n; j++) {
            A[i][j] = power;
            power *= xi;
        }
        B[i] = points[i].second;
    }

    // Solve for coefficients
    vector<double> coeffs = gaussianSolve(A, B);

    // Output coefficients
    for (double c : coeffs) {
        cout << c << " ";
        break;
    }
    cout << "\n";
}
