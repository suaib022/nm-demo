#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

double calculateU(double u, int n) {
    double temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u + i);
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
        for (int j = n - 1; j >= i; j--)
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
    }

    cout << endl << "Backward Difference Table:" << endl;
    for (int i = 0; i < n; i++) {
        cout << setw(10) << x[i];
        for (int j = 0; j <= i; j++)
            cout << setw(10) << fixed << setprecision(2) << y[i][j];
        cout << endl;
    }

    double value;
    cout << endl << "Enter the value to interpolate: ";
    cin >> value;

    double sum = y[n - 1][0];
    double u = (value - x[n - 1]) / (x[1] - x[0]);

    for (int i = 1; i < n; i++) {
        sum = sum + (calculateU(u, i) * y[n - 1][i]) / fact(i);
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
            for (int j = n - 1; j >= i; j--)
                y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }

        cout << endl << "Updated Backward Difference Table:" << endl;
        for (int i = 0; i < n; i++) {
            cout << setw(10) << x[i];
            for (int j = 0; j <= i; j++)
                cout << setw(10) << fixed << setprecision(2) << y[i][j];
            cout << endl;
        }

        sum = y[n - 1][0];
        u = (value - x[n - 1]) / (x[1] - x[0]);

        for (int i = 1; i < n; i++) {
            sum = sum + (calculateU(u, i) * y[n - 1][i]) / fact(i);
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

Backward Difference Table:
     35.00     31.00
     45.00     73.00     42.00
     55.00    124.00     51.00      9.00
     65.00    159.00     35.00    -16.00    -25.00
     75.00    190.00     31.00     -4.00     12.00     37.00

Enter the value to interpolate: 70

Value at 70.0000 is 172.8047

Use additional data point? (y/n): y
Enter the new data point (x and y): 85 215

Updated Backward Difference Table:
   35.0000     31.00
     45.00     73.00     42.00
     55.00    124.00     51.00      9.00
     65.00    159.00     35.00    -16.00    -25.00
     75.00    190.00     31.00     -4.00     12.00     37.00
     85.00    215.00     25.00     -6.00     -2.00    -14.00    -51.00

New Value at 70.0000 is 174.1992
Error: 1.3945
*/