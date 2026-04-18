#include <bits/stdc++.h>
using namespace std;

// Safe string → double
bool safeParse(const string& s, double& value) {
    try {
        size_t idx;
        value = stod(s, &idx);
        return idx == s.size();
    } catch (...) {
        return false;
    }
}

// Read CSV (first column only, skip header/bad rows)
vector<double> readCSV(string filename) {
    ifstream file(filename);
    vector<double> data;

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string val;
        getline(ss, val, ',');

        double num;
        if (safeParse(val, num)) {
            data.push_back(num);
        }
    }

    return data;
}

// Compute MSE
double computeMSE(const vector<double>& A, const vector<double>& B) {
    int n = min(A.size(), B.size());
    if (n == 0) return 0.0;

    double mse = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = A[i] - B[i];
        mse += diff * diff;
    }
    return mse / n;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: ./error k\n";
        return 1;
    }

    int k = stoi(argv[1]);

    vector<vector<double>> allData(k + 1);

    for (int i = 1; i <= k; i++) {
        string filename = to_string(i) + ".csv";
        allData[i] = readCSV(filename);

        if (allData[i].empty()) {
            cerr << "Warning: " << filename << " has no valid data\n";
        }
    }

    // Pairwise MSE
    for (int i = 1; i <= k; i++) {
        for (int j = i + 1; j <= k; j++) {
            double mse = computeMSE(allData[i], allData[j]);
            cout << i << "," << j << " -> " << mse << endl;
        }
    }

    return 0;
}