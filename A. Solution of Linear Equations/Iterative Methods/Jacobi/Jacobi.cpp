#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream infile("inputjacobi.txt");
    ofstream outfile("outputjacobi.txt");

    if (!infile || !outfile)
        return 1;
    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);
        vector<double> x(n, 0.0);
        vector<double> new_x(n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                infile >> a[i][j];
            }
            infile >> b[i];
        }

        double tol;
        int max_iter;
        infile >> tol >> max_iter;

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Iteration\t";
        for (int i = 0; i < n; i++)
            outfile << "x" << i + 1 << "\t\t";
        outfile << "Error" << endl;
        outfile << endl;

        int iter = 0;
        double error = 0.0;

        do
        {
            error = 0.0;
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += a[i][j] * x[j];
                    }
                }
                new_x[i] = (b[i] - sum) / a[i][i];
                double diff = abs(new_x[i] - x[i]);
                if (diff > error)
                {
                    error = diff;
                }
            }

            x = new_x;
            iter++;

            outfile << iter << "\t\t";
            for (int i = 0; i < n; i++)
                outfile << fixed << setprecision(3) << x[i] << "\t";
            outfile << scientific << setprecision(3) << error << endl;

        } while (error > tol && iter < max_iter);

        outfile << endl;
        if (iter < max_iter)
        {
            outfile << "Converged in " << iter << " iterations." << endl;
        }
        else
        {
            outfile << "Maximum iterations reached." << endl;
        }

        outfile << "Solution:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << "x" << i + 1 << " = " << fixed << setprecision(3) << x[i] << endl;
        }
        outfile << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}
