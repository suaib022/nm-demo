#include <bits/stdc++.h>
using namespace std;

void print_matrix(vector<vector<double>> &a, ostream &out)
{
    int n = a.size();
    int m = a[0].size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            out << fixed << setprecision(6) << a[i][j] << " ";
        }
        out << "\n";
    }
    out << "\n";
}

int gauss_jordan(vector<vector<double>> a, vector<double> &root, ostream &out)
{
    int n = a.size();

    out << "Initial Augmented Matrix:\n";
    print_matrix(a, out);

    for (int i = 0; i < n; i++)
    {
        int pivot = i;

        for (int r = i + 1; r < n; r++)
            if (fabs(a[r][i]) > fabs(a[pivot][i]))
                pivot = r;

        if (fabs(a[pivot][i]) < 1e-12)
            continue;

        swap(a[i], a[pivot]);

        double div = a[i][i];
        for (int c = 0; c <= n; c++)
            a[i][c] /= div;

        for (int r = 0; r < n; r++)
        {
            if (r == i)
                continue;
            double factor = a[r][i];
            for (int c = 0; c <= n; c++)
                a[r][c] -= factor * a[i][c];
        }

        out << "After pivot normalization in row " << i + 1 << ":\n";
        print_matrix(a, out);
    }

    for (int i = 0; i < n; i++)
    {
        bool zeroRow = true;
        for (int j = 0; j < n; j++)
            if (fabs(a[i][j]) > 1e-12)
                zeroRow = false;

        if (zeroRow && fabs(a[i][n]) > 1e-12)
            return -1;
        if (zeroRow && fabs(a[i][n]) < 1e-12)
            return 0;
    }

    root.assign(n, 0.0);
    for (int i = 0; i < n; i++)
        root[i] = a[i][n];

    return 1;
}

int main()
{
    ifstream fin("inputgauss.txt");
    ofstream fout("outputgauss.txt");

    int T;
    fin >> T;

    for (int t = 1; t <= T; t++)
    {
        int n;
        fin >> n;

        vector<vector<double>> a(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> a[i][j];

        fout << "System " << t << " (Gauss-Jordan)\n";

        vector<double> root;
        int status = gauss_jordan(a, root, fout);

        if (status == 1)
        {
            fout << "Unique Solution:\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << root[i] << "\n";
        }
        else if (status == 0)
            fout << "Infinite Solutions\n";
        else
            fout << "No Solution\n";

        fout << "\n\n";
    }

    return 0;
}
