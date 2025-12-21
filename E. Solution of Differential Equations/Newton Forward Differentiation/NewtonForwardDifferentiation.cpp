#include <bits/stdc++.h>
using namespace std;
typedef vector<double> Poly;

Poly multiply(const Poly &p, double constant)
{
    Poly res(p.size() + 1, 0.0);
    for (size_t i = 0; i < p.size(); i++)
    {
        res[i + 1] += p[i];
        res[i] += constant * p[i];
    }
    return res;
}

Poly differentiate(const Poly &p)
{
    if (p.size() <= 1)
        return {0.0};
    Poly res;
    for (size_t i = 1; i < p.size(); i++)
    {
        res.push_back(p[i] * i);
    }
    return res;
}

double evaluate(const Poly &p, double x)
{
    double res = 0;
    double x_pow = 1;
    for (double c : p)
    {
        res += c * x_pow;
        x_pow *= x;
    }
    return res;
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
    ifstream infile("inputdiff.txt");
    ofstream outfile("outputdiff.txt");

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
        outfile << "Forward Difference Table:" << endl;
        for (int i = 0; i < n; i++)
        {
            outfile << setw(10) << x[i];
            for (int j = 0; j < n - i; j++)
                outfile << setw(10) << fixed << setprecision(4) << y[i][j];
            outfile << endl;
        }

        double value;
        int order;
        infile >> value >> order;

        double h = x[1] - x[0];
        double u = (value - x[0]) / h;

        double derivative_val = 0;
        Poly P = {1.0};

        for (int k = 0; k < n; k++)
        {
            if (k > 0)
            {
                P = multiply(P, -(double)(k - 1));
            }

            if (k >= order)
            {
                Poly P_deriv = P;
                for (int m = 0; m < order; m++)
                {
                    P_deriv = differentiate(P_deriv);
                }

                double term_val = evaluate(P_deriv, u);
                derivative_val += (y[0][k] / fact(k)) * term_val;
            }
        }

        derivative_val /= pow(h, order);

        outfile << endl
                << "For x = " << value << ":" << endl;
        outfile << "Derivative Order: " << order << endl;
        outfile << "Result: " << fixed << setprecision(3) << derivative_val << endl;
        outfile << endl
                << endl;
    }
    infile.close();
    outfile.close();

    return 0;
}