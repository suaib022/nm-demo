#include <bits/stdc++.h>
using namespace std;

double calculateU(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u - i);
    return temp;
}

int fact(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream infile("inputnewtforward.txt");
    ofstream outfile("outputnewtforward.txt");

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
            for (int j = 0; j < n - i; j++)
                y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }

        outfile << "Case " << caseNum++ << ":" << endl;
        outfile << "Difference Table:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << setw(10) << x[i];
            for (int j = 0; j < n - i; j++)
                outfile << setw(10) << fixed << setprecision(4) << y[i][j];
            outfile << endl;
        }

        double value;
        infile >> value;

        double sum = y[0][0];
        double u = (value - x[0]) / (x[1] - x[0]);

        for (int i = 1; i < n; i++)
        {
            sum = sum + (calculateU(u, i) * y[0][i]) / fact(i);
        }

        outfile << endl
                << "Value at " << value << " is " << fixed << setprecision(3) << sum << endl;
        outfile << endl;
    }
    infile.close();
    outfile.close();
    return 0;
}