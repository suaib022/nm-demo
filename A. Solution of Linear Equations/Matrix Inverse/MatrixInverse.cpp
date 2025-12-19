#include <bits/stdc++.h>

using namespace std;

void getCofactor(const vector<vector<double>> &A, vector<vector<double>> &temp, int p, int q, int n)
{
    int i = 0, j = 0;
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double determinant(const vector<vector<double>> &A, int n)
{
    if (n == 1)
        return A[0][0];

    double det = 0;
    vector<vector<double>> temp(n, vector<double>(n));
    int sign = 1;

    for (int f = 0; f < n; f++)
    {
        getCofactor(A, temp, 0, f, n);
        det += sign * A[0][f] * determinant(temp, n - 1);
        sign = -sign;
    }
    return det;
}

void adjoint(const vector<vector<double>> &A, vector<vector<double>> &adj, int n, ofstream &outfile)
{
    if (n == 1)
    {
        adj[0][0] = 1;
        outfile << "\nCofactor Matrix:" << endl
                << fixed << setprecision(0) << 1.0 << endl;
        outfile << "\nAdjoint Matrix:" << endl
                << fixed << setprecision(0) << 1.0 << endl;
        return;
    }

    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n));
    vector<vector<double>> cof(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            getCofactor(A, temp, i, j, n);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            cof[i][j] = (sign) * (determinant(temp, n - 1));
            adj[j][i] = cof[i][j];
        }
    }

    outfile << "\nCofactor Matrix:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            outfile << fixed << setprecision(0) << cof[i][j] << "\t";
        outfile << endl;
    }

    outfile << "\nAdjoint Matrix:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            outfile << fixed << setprecision(0) << adj[i][j] << "\t";
        outfile << endl;
    }
}

bool inverse(const vector<vector<double>> &A, vector<vector<double>> &inv, int n, ofstream &outfile)
{
    double det = determinant(A, n);
    outfile << "\nDeterminant: " << fixed << setprecision(0) << det << endl;

    if (abs(det) < 1e-9)
    {
        return false;
    }

    vector<vector<double>> adj(n, vector<double>(n));
    adjoint(A, adj, n, outfile);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] / det;

    return true;
}

int main()
{
    ifstream infile("inputinverse.txt");
    ofstream outfile("outputinverse.txt");

    if (!infile || !outfile)
        return 1;

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<vector<double>> A(n, vector<double>(n));
        vector<vector<double>> inv(n, vector<double>(n));

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                infile >> A[i][j];
            }
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Given Matrix:" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                outfile << fixed << setprecision(0) << A[i][j] << "\t";
            }
            outfile << endl;
        }

        if (inverse(A, inv, n, outfile))
        {
            outfile << "\nInverse Matrix:" << endl;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    outfile << fixed << setprecision(6) << inv[i][j] << "\t";
                }
                outfile << endl;
            }
        }
        else
        {
            outfile << "\nSingular Matrix, Inverse does not exist." << endl;
        }
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}
