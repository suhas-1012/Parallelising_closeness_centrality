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
    if (argc < 3) {
        cout << "Usage: ./error k\n";
        return 1;
    }

    int a = stoi(argv[1]);
    int b = stoi(argv[2]);

    vector<vector<double>> allData(b-a + 2);

    for (int i = a; i <= b; i++) {
        string filename = to_string(i) + ".csv";
        allData[i-a+1] = readCSV(filename);

        if (allData[i-a+1].empty()) {
            cerr << "Warning: " << filename << " has no valid data\n";
        }
    }

    // Pairwise MSE
    for (int i = a; i <= b; i++) {
        for (int j = i + 1; j <= b; j++) {
            double mse = computeMSE(allData[i-a+1], allData[j-a+1]);
            cout << i << "," << j << " -> " << mse << endl;
        }
    }

    return 0;
}