#include <bits/stdc++.h>
using namespace std;

double f(const vector<double> &coeffs, double x)
{
    double result = 0.0;
    int n = coeffs.size();
    for (int i = 0; i < n; i++)
    {
        result += coeffs[i] * pow(x, n - i - 1);
    }
    return result;
}

bool stop(double prev, double curr, double tol)
{
    return fabs(prev - curr) < tol;
}

vector<double> Bisection(const vector<double> &coeffs, double l, double u, double tol, double step)
{
    vector<pair<double, double>> intervals;
    vector<double> roots;

    for (double x = l; x <= u - step; x += step)
    {
        double fx = f(coeffs, x);
        double fx_next = f(coeffs, x + step);

        if (fabs(fx) < 1e-12)
        {
            roots.push_back(x);
            continue;
        }

        if (fx * fx_next < 0)
        {
            intervals.push_back({x, x + step});
        }
    }

    for (auto interval : intervals)
    {
        double a = interval.first;
        double b = interval.second;
        double fa = f(coeffs, a);
        double fb = f(coeffs, b);

        if (fa * fb > 0)
            continue;

        double prev = a;
        double mid = (a + b) / 2.0;
        int iter = 0;

        while (!stop(prev, mid, tol))
        {
            double fmid = f(coeffs, mid);
            if (fabs(fmid) < 1e-12)
                break;

            if (fa * fmid < 0)
            {
                b = mid;
                fb = fmid;
            }
            else
            {
                a = mid;
                fa = fmid;
            }

            prev = mid;
            mid = (a + b) / 2.0;
            iter++;
        }

        roots.push_back(mid);
    }

    sort(roots.begin(), roots.end());
    roots.erase(unique(roots.begin(), roots.end(), [tol](double a, double b)
                       { return fabs(a - b) < tol; }),
                roots.end());

    return roots;
}

string print_function(const vector<double> &coeffs)
{
    stringstream out;
    int n = coeffs.size();
    out << "f(x)=";
    for (int i = 0; i < n; i++)
    {
        if (i != 0 && coeffs[i] >= 0)
            out << "+";
        out << coeffs[i] << "x^" << n - i - 1 << " ";
    }
    return out.str();
}

int main()
{
    ifstream infile("inputbisection.txt");
    ofstream outfile("outputbisection.txt");

    if (!infile)
    {
        return 1;
    }

    int numProblems;
    infile >> numProblems;

    for (int p = 0; p < numProblems; p++)
    {
        int degree;
        infile >> degree;

        vector<double> coeffs(degree + 1);
        for (int i = 0; i <= degree; i++)
        {
            infile >> coeffs[i];
        }

        double l, u, tol, step;
        infile >> l >> u >> tol >> step;

        outfile << "Problem " << p + 1 << endl;
        outfile << "Function: " << print_function(coeffs) << endl;
        outfile << "Interval: [" << l << ", " << u << "], Tolerance: " << tol << ", Step: " << step << endl;

        vector<double> roots = Bisection(coeffs, l, u, tol, step);

        outfile << "Roots found : " << endl;
        if (!roots.empty())
        {
            for (double root : roots)
            {
                outfile << fixed << setprecision(6) << root << endl;
            }
        }
        outfile << "\n\n"
                << endl;
    }

    infile.close();
    outfile.close();
    return 0;
}
