// highload-mathlibs.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <conio.h>
#include <math.h>
const int n = 6;
const int N = 10000;

using namespace std;

double f(double x)
{
    return (2 * x - 1) * (1. / pow((x * (1 - x)), 5.)) + 0.1 * pow(x, 2.);
}

double Gauss(double a, double b)
{
    const double Xi[n] = { -0.9324700, -0.6612094, -0.2386142, 0.2386142, 0.6612094, 0.9324700 };
    const double Ci[n] = { 0.1713245, 0.3607616, 0.4679140, 0.4679140, 0.3607616, 0.1713245 };

    double ra = (b - a) / 2;
    double su = (a + b) / 2;
    double Q, S = 0.0;
    for (int i = 0; i < n; i++)
    {
        Q = su + ra * Xi[i];
        S += Ci[i] * f(Q);
    }
    return ra * S;
}

int main()
{
    double a = 2.0;
    double b = 10.0;
    double s = 0.0;
    for (int i = 0; i < N; ++i)
    {
        s += Gauss(a + i * (b - a) / N, a + (i + 1) * (b - a) / N);
    }
    cout << "I = " << s << endl;

    return 0;
}
