#include <bits/stdc++.h>

using namespace std;

double calculateU(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u + i);
    return temp;
}

double fact(int n)
{
    double f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream infile("inputbackward.txt");
    ofstream outfile("outputbackward.txt");

    if (!infile || !outfile)
    {
        return 1;
    }

    int n;
    int caseNum = 1;

    while (infile >> n)
    {
        vector<double> x(n);
        vector<vector<double>> y(n, vector<double>(n));

        for (int i = 0; i < n; i++)
        {
            infile >> x[i] >> y[i][0];
        }

        for (int i = 1; i < n; i++)
        {
            for (int j = n - 1; j >= i; j--)
                y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Backward Difference Table:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << setw(10) << x[i];
            for (int j = 0; j <= i; j++)
                outfile << setw(10) << fixed << setprecision(4) << y[i][j];
            outfile << endl;
        }

        double value;
        infile >> value;

        double sum = y[n - 1][0];
        double u = (value - x[n - 1]) / (x[1] - x[0]);
        double last_term = 0.0;

        for (int i = 1; i < n; i++)
        {
            last_term = (calculateU(u, i) * y[n - 1][i]) / fact(i);
            sum = sum + last_term;
        }

        outfile << endl
                << "Value at " << value << " is " << fixed << setprecision(6) << sum << endl;

        if (sum != 0)
        {
            double relative_error = abs(last_term / sum) * 100.0;
            outfile << "Approximate Relative Error: " << relative_error << "%" << endl;
        }
        else
        {
            outfile << "Approximate Error: " << abs(last_term) << " (Cannot calculate percentage relative error as sum is 0)" << endl;
        }

        outfile << endl
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}
