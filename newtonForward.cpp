#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

double calculateU(double u, int n) {
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u - i);
    return temp;
}

int fact(int n) {
    int f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    int n;
    cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n);
    vector<vector<double>> y(n, vector<double>(n));

    cout << "Enter data points (x and y):" << endl;
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i][0];
    }

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++)
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
    }

    cout << endl << "Difference Table:" << endl;
    for (int i = 0; i < n; i++) {
        cout << setw(10) << x[i];
        for (int j = 0; j < n - i; j++)
            cout << setw(10) << fixed << setprecision(2) << y[i][j];
        cout << endl;
    }

    double value;
    cout << endl << "Enter the value to interpolate: ";
    cin >> value;

    double sum = y[0][0];
    double u = (value - x[0]) / (x[1] - x[0]);

    for (int i = 1; i < n; i++) {
        sum = sum + (calculateU(u, i) * y[0][i]) / fact(i);
    }

    cout << endl << "Value at " << value << " is " << fixed << setprecision(4) << sum << endl;
    double initial_ans = sum;

    char choice;
    cout << endl << "Use additional data point? (y/n): ";
    cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        n++;
        x.resize(n);
        y.resize(n, vector<double>(n));

        cout << "Enter the new data point (x and y): ";
        cin >> x[n - 1] >> y[n - 1][0];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < n - i; j++)
                y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }

        cout << endl << "Updated Difference Table:" << endl;
        for (int i = 0; i < n; i++) {
            cout << setw(10) << x[i];
            for (int j = 0; j < n - i; j++)
                cout << setw(10) << fixed << setprecision(2) << y[i][j];
            cout << endl;
        }

        sum = y[0][0];
        u = (value - x[0]) / (x[1] - x[0]);

        for (int i = 1; i < n; i++) {
            sum = sum + (calculateU(u, i) * y[0][i]) / fact(i);
        }

        cout << endl << "New Value at " << value << " is " << fixed << setprecision(4) << sum << endl;
        cout << "Error: " << abs(sum - initial_ans) << endl;
    }

    return 0;
}

/*
Sample Input/Output:

Enter number of data points: 5
Enter data points (x and y):
35 31
45 73
55 124
65 159
75 190

Difference Table:
     35.00     31.00     42.00      9.00    -25.00     37.00
     45.00     73.00     51.00    -16.00     12.00
     55.00    124.00     35.00     -4.00
     65.00    159.00     31.00
     75.00    190.00

Enter the value to interpolate: 40

Value at 40.0000 is 47.8672

Use additional data point? (y/n): y
Enter the new data point (x and y): 85 215

Updated Difference Table:
     35.00     31.00     42.00      9.00    -25.00     37.00    -51.00
     45.00     73.00     51.00    -16.00     12.00    -14.00
     55.00    124.00     35.00     -4.00     -2.00
     65.00    159.00     31.00     -6.00
     75.00    190.00     25.00
     85.00    215.00

New Value at 40.0000 is 46.4727
Error: 1.3945
*/