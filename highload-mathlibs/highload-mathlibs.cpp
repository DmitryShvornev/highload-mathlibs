// highload-mathlibs.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <iomanip>
#include <conio.h>
#include <math.h>
const int n = 6;

const int N = 10000;

using namespace std;

// функция f  в методе Гаусса 

double f(double x)
{
    return (2 * x - 1) * (1. / pow((x * (1 - x)), 5.)) + 0.1 * pow(x, 2.);
}

// функция f  в процессе ортогонализации Грамма-Шмидта 

double f_ort(double x, int n)
{
	return  pow(x, n);
}

// Реализация метода Гаусса

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

//вычисление интеграла методом Гаусса

void integral_calculation() {
	double a = 0.0;
	double b = 10000.0;
	double s = 0.0;
	for (int i = 0; i < N; ++i)
	{
		s += Gauss(a + i * (b - a) / N, a + (i + 1) * (b - a) / N);
	}
	cout << "I = " << s << endl;
}

// процесс ортогонализации Грамма-Шмидта:
 
void orthogonalization() {
	//сетка 
	int n_ort = 100;
	double f_norm, a = -1;
	double b = 1;
	double h = (b - a) / n_ort;

	int N_ort = 50;
	double ** ort_func = new double*[N_ort];
	double * alpha_i_j = new double[N_ort];
	for (int i = 0; i < N_ort; ++i)
		ort_func[i] = new double[n_ort];
	
	auto scalar_prod = [&](int i_idx, int j_idx) {
		double summ = 0;
		for (int i = 0; i < n_ort ; ++i) 
			summ += ort_func[i_idx][i] * ort_func[j_idx][i] * h;
		return summ;
	};

	cout << "Gram-Schmidt orthogonalization process:" << endl;

	for (int i = 0; i < N_ort; ++i) {
		for (int j = 0; j < n_ort ; ++j)
			ort_func[i][j] = f_ort(a + j * h + h / 2, i);
		for (int j = 0; j < i; ++j)
			alpha_i_j [j]= -scalar_prod(i, j);
		for (int j = 0; j < i; ++j) {
			for (int k = 0; k < n_ort ; ++k)
				ort_func[i][k] += alpha_i_j[j] * ort_func[j][k];
		}
		f_norm = sqrt(scalar_prod(i, i));
		for (int j = 0; j < n_ort ; ++j)
			ort_func[i][j] = ort_func[i][j] / f_norm;
	}

	// вывод на печать значений скалярных произведений - функции высокой степени уже не ортогональны...
	// надо бы сделать адекватный вывод
	for (int i = 0; i < N_ort; ++i) {
		cout << "scalar_prod for f_" << i << " : ";
		cout << fixed;
		cout.precision(5);
		for(int j = i; j< N_ort; ++j)
			cout << scalar_prod(i, j)<< " ";
		cout << endl;
	}
	

	delete[] alpha_i_j;
	for (int i = 0; i < N_ort; ++i) {
		delete[] ort_func[i];
	}
	delete[] ort_func;
}

void print_matrix(double** mat, int n) {
	for (int i = 0; i < n; i++)
	{
		cout << "row" << " " << i << ":";
		for (int j = 100; j < 2 * n; j++)
		{
			cout << setw(2) << mat[i][j] << " ";
		}
		cout << endl;
		cout << endl;
	}
}

// вычисление обратной матрицы методом гаусса

void matrix_inverse_by_gauss() {
	int N_inv = 100;
	double coeff;

	double **A_matr = new double *[N_inv];
	double** A_matr_copy = new double* [N_inv];
	double** E = new double* [N_inv];
	for (int i = 0; i < N_inv; ++i) {
		A_matr[i] = new double[2 * N_inv];
		A_matr_copy[i] = new double[N_inv];
		E[i] = new double[N_inv];
	}
		
	for (int i = 0; i < N_inv; ++i)
		for (int j = 0; j < N_inv; ++j) {
			A_matr[i][j] = pow(double(i) / 100, j);
			A_matr_copy[i][j] = pow(double(i) / 100, j);
			if (i == j) {
				A_matr[i][N_inv + j] = 1;
			}
			else {
				A_matr[i][N_inv + j] = 0;
			}
		}
	// вычисление обратной матрицы
	for (int i = 0; i < N_inv; ++i) {
		coeff = A_matr[i][i];
		for (int j = 0; j < 2 * N_inv; ++j) {
			A_matr[i][j] = A_matr[i][j] / coeff;
		}
		for (int j = 0; j < N_inv; ++j) {
			if (i != j) {
				coeff = A_matr[j][i];
				for (int k = 0; k < 2 * N_inv; ++k) {
					A_matr[j][k] -= coeff * A_matr[i][k];
				}	
			}
		}
	}

	// проверка-перемножение с исходной - оценка влияния погрешностей.
	for (int i = 0; i < N_inv; i++)
	{
		for (int j = 0; j < N_inv; j++)
		{
			for (int k = 0; k < N_inv; k++)
			{
				E[i][j] += A_matr[i][k + N_inv] * A_matr_copy[k][j];
			}
		}

	}

	print_matrix(E, N_inv);

	for (int i = 0; i < N_inv; ++i) {
		delete[] A_matr[i];
		delete[] A_matr_copy[i];
		delete[] E[i];
	}
	delete[] A_matr;
	delete[] A_matr_copy;
	delete[] E;
}


int main()
{

	//integral_calculation();
	//orthogonalization();
	matrix_inverse_by_gauss();

    return 0;
}
