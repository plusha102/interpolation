#include <iostream>
#include <cstdio>
#include <conio.h>
#include <cstdlib>
#include <vector>
#include <array>

using namespace std;

int m, n;
double w, l;

//квадратичная интерполяция
//формулы для коэффициантов a1, a2, a3, a4 выведены из уравнения
//z = a0 + a1*x + a2*y + a3*x*y
//по четырем точкам


double fa0(int i, int j, vector<double>& input_z)
{
	return i * j*(input_z[i*m + j] - input_z[i*m + j + 1] - input_z[(i + 1)*m + j] + input_z[(i + 1)*m + j + 1]) + i * (input_z[i*m + j] - input_z[i*m + j + 1]) + j * (input_z[i*m + j] - input_z[(i + 1)*m + j]) + input_z[i*m + j];
}

double fa1(int i, int j, vector<double>& input_z)
{
	return -n * (j*(input_z[i*m + j] - input_z[i*m + j + 1] - input_z[(i + 1)*m + j] + input_z[(i + 1)*m + j + 1]) + input_z[i*m + j] - input_z[i*m + j + 1]) / l;
}

double fa2(int i, int j, vector<double>& input_z)
{
	return -m * (i*(input_z[i*m + j] - input_z[i*m + j + 1] - input_z[(i + 1)*m + j] + input_z[(i + 1)*m + j + 1]) + input_z[i*m + j] - input_z[(i + 1)*m + j]) / w;
}

double fa3(int i, int j, vector<double>& input_z)
{
	return m * n*(input_z[i*m + j] - input_z[i*m + j + 1] - input_z[(i + 1)*m + j] + input_z[(i + 1)*m + j + 1]) / (l*w);
}

class Q_Interpol
{
	int nodes_x, nodes_y;
	double widht, lenght;

	vector<double> a0;
	vector<double> a1;
	vector<double> a2;
	vector<double> a3;
	vector<double> input_z;

public: 
	
	Q_Interpol(int nx, int ny, double w, double l, vector<double>& input_z)
{
	nodes_x = nx;
	nodes_y = ny;
	widht = w;
	lenght = l;

	for (int i = 0; i < nodes_x; i += 2) // двойной шаг по строчкам
	{
		for (int j = 0; j < nodes_y; j += 2) // двойной шаг по столбцам
		{
			a0.push_back(fa0(i, j, input_z));
			a1.push_back(fa1(i, j, input_z));
			a2.push_back(fa2(i, j, input_z));
			a3.push_back(fa3(i, j, input_z));

			if (i%nodes_x & (i == nodes_x - 2)) // если m - нечетное число, то последний квадрат в столбце строится по m-1 и m
			{
				a0.push_back(fa0(i + 1, j, input_z));
				a1.push_back(fa1(i + 1, j, input_z));
				a2.push_back(fa2(i + 1, j, input_z));
				a3.push_back(fa3(i + 1, j, input_z));
			}

			if (j%nodes_y & (j == nodes_y - 2)) // если n - нечетное число, то последний квадрат в строчке строится по n-1 и n
			{
				a0.push_back(fa0(i, j + 1, input_z));
				a1.push_back(fa1(i, j + 1, input_z));
				a2.push_back(fa2(i, j + 1, input_z));
				a3.push_back(fa3(i, j + 1, input_z));

				if (i%nodes_x & (i == nodes_x - 2)) // n и m - нечетные числа, последний квадрат матрицы строится по m-1, m, n-1, n
				{
					a0.push_back(fa0(i + 1, j + 1, input_z));
					a1.push_back(fa1(i + 1, j + 1, input_z));
					a2.push_back(fa2(i + 1, j + 1, input_z));
					a3.push_back(fa3(i + 1, j + 1, input_z));
				}
			}
		}
	}
}
		//поиск z по x и y
		double z(double x, double y)
		{
			int dx = int((x*nodes_x) / widht);
			int dy = int((y*nodes_y) / lenght);

			return a0[dx*nodes_x + dy] + a1[dx*nodes_x + dy] * x + a2[dx*nodes_x + dy] * y + a3[dx*nodes_x + dy] * x*y;

			throw "error in q_interpolation";
		}
};

//линейная интерполяция
//коэффициенты a и b выведены из уравнения
// y = a*x + b
//по двум точкам

class L_Interpol
{
	int nodes;
	double lenght;

	vector<double> a;
	vector<double> b;
	vector<double> input_x;
	vector<double> input_y;

public:

	L_Interpol(int n, double l, vector<double>& input_x, vector<double>& input_y)
	{
		a.resize(n);
		b.resize(n);

		this->input_x.resize(n);
		this->input_y.resize(n);

		nodes = n;
		lenght = l;

		this->input_x = input_x;
		this->input_y = input_y;

		for (int i = 0; i < nodes - 1; i++)
		{
			a[i] = (input_y[i + 1] - input_y[i]) / (input_x[i + 1] - input_x[i]);
			b[i] = input_y[i] - a[i] * input_x[i];
		}
	}

	double y(double x)
	{
		int i;
		for (i = 0; i < nodes - 1; i++)
		{
			if (x == input_x[i])
				return input_y[i];
		}

		for (i = 0; i < nodes - 1; i++)
		{
			if (x >= input_x[i] && x <= input_x[i + 1])
				return a[i] * x + b[i];
		}

		throw "error in l_interpolation";
	}
};


//тест
void main()
{
	m = 2, n = 5, w = 2., l = 5.;

	vector<double> x;
	vector<double> y;
	vector<double> z;

	x.resize(n);
	y.resize(n);
	z.resize(n*m);

	for (int i = 0; i < n; i++)
	{
		x.push_back(rand() % 10);
		y.push_back(rand() % 10);
	}

	for (int i = 0; i < n*m; i++)
	{
		z.push_back(rand() % 10);
	}

	Q_Interpol test1(m, n, w, l, z);
	L_Interpol test2(n, l, x, y);

	_getch();
	
}
