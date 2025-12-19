#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream infile("inputlu.txt");
    ofstream outfile("outputlu.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                infile >> a[i][j];
            infile >> b[i];
        }

        vector<vector<double>> l(n, vector<double>(n, 0));
        vector<vector<double>> u(n, vector<double>(n, 0));

        for (int i = 0; i < n; i++)
        {
            for (int k = i; k < n; k++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += l[i][j] * u[j][k];
                u[i][k] = a[i][k] - sum;
            }

            for (int k = i; k < n; k++)
            {
                if (i == k)
                    l[i][i] = 1;
                else
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += l[k][j] * u[j][i];

                    if (abs(u[i][i]) > 1e-9)
                        l[k][i] = (a[k][i] - sum) / u[i][i];
                    else
                        l[k][i] = 0;
                }
            }
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Lower Triangular Matrix :" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                outfile << fixed << setprecision(3) << setw(8) << l[i][j] << " ";
            outfile << endl;
        }
        outfile << endl;

        outfile << "Upper Triangular Matrix :" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                outfile << fixed << setprecision(3) << setw(8) << u[i][j] << " ";
            outfile << endl;
        }

        vector<double> y(n);
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += l[i][j] * y[j];
            y[i] = b[i] - sum;
        }

        bool unique = true;
        bool inconsistent = false;

        for (int i = 0; i < n; i++)
        {
            if (abs(u[i][i]) < 1e-9)
            {
                unique = false;
                bool allZero = true;
                for (int j = i + 1; j < n; j++)
                {
                    if (abs(u[i][j]) > 1e-9)
                    {
                        allZero = false;
                        break;
                    }
                }

                if (allZero && abs(y[i]) > 1e-9)
                {
                    inconsistent = true;
                    break;
                }
            }
        }

        if (unique)
        {
            vector<double> x(n);
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i + 1; j < n; j++)
                    sum += u[i][j] * x[j];
                x[i] = (y[i] - sum) / u[i][i];
            }

            outfile << "\nSolution:" << endl;
            for (int i = 0; i < n; i++)
                outfile << "x" << i + 1 << " = " << fixed << setprecision(4) << x[i] << endl;

            outfile << "\nThe system has a unique solution." << endl;
        }
        else
        {
            if (inconsistent)
                outfile << "\nThe system has no solution." << endl;
            else
                outfile << "\nThe system has infinite solutions." << endl;
        }
        outfile << endl
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}
