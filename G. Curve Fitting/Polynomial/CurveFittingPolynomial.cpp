#include <bits/stdc++.h>

using namespace std;

void gaussianElimination(vector<vector<double>> &a, vector<double> &b, vector<double> &x, int n)
{
    for (int i = 0; i < n; i++)
    {
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (abs(a[k][i]) > abs(a[maxRow][i]))
            {
                maxRow = k;
            }
        }

        swap(a[i], a[maxRow]);
        swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; k++)
        {
            double factor = a[k][i] / a[i][i];
            for (int j = i; j < n; j++)
            {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}

int main()
{
    ifstream infile("inputpoly.txt");
    ofstream outfile("outputpoly.txt");

    if (!infile || !outfile)
        return 1;

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i];
        }

        int degree;
        infile >> degree;
        int m = degree + 1;

        vector<vector<double>> A(m, vector<double>(m, 0.0));
        vector<double> B(m, 0.0);
        vector<double> coeffs(m);

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                {
                    sum += pow(x[k], i + j);
                }
                A[i][j] = sum;
            }
            double sum = 0;
            for (int k = 0; k < n; k++)
            {
                sum += y[k] * pow(x[k], i);
            }
            B[i] = sum;
        }

        gaussianElimination(A, B, coeffs, m);

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Data Points: " << n << ", Degree: " << degree << endl;
        outfile << "Polynomial Fit Equation: y = ";
        for (int i = 0; i < m; i++)
        {
            if (i == 0)
                outfile << fixed << setprecision(4) << coeffs[i];
            else if (coeffs[i] >= 0)
                outfile << " + " << coeffs[i] << "x^" << i;
            else
                outfile << " - " << abs(coeffs[i]) << "x^" << i;
        }
        outfile << endl;
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}
