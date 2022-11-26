#pragma once
#include<iostream>
#include<cmath>
#include"Complex.hpp"
#include"Matrix.hpp"
using namespace std;

//声明
Matrix<double> HessenBerg_Transformation(Matrix<double>);
double* HessenMatQR(Matrix<double>, int, double);
Complex* EigenValue(const Matrix<double>&);

//定义
Matrix<double> HessenBerg_Transformation(Matrix<double> A)
{
	if (A.m != A.n)
	{
		cout << "HessenBerg_Transformation Failed.Not a square matrix!\n";
		return Matrix<double>();
	}
	int n = A.n;
	int i, j, k;
	double d, t;
	for (k = 2; k <= n - 1; k++)
	{
		d = 0.0;
		for (j = k; j <= n; j++)
		{
			t = A(j, k - 1);
			if (fabs(t) > fabs(d))
			{
				d = t;
				i = j;
			}
		}
		if (fabs(d) + 1.0 != 1.0)
		{
			if (i != k)
			{
				for (j = k - 1; j <= n; j++)
				{
					swap(A(i, j), A(k, j));
				}
				for (j = 1; j <= n; j++)
				{
					swap(A(j, i), A(j, k));
				}
			}
			for (i = k + 1; i <= n; i++)
			{
				t = A(i, k - 1) / d;
				A(i, k - 1) = 0.0;
				for (j = k; j <= n; j++)
				{
					A(i, j) = A(i, j) - t * A(k, j);
				}
				for (j = 1; j <= n; j++)
				{
					A(j, k) = A(j, k) + t * A(j, i);
				}
			}
		}
	}
	return A;
}

double* HessenMatQR(Matrix<double> A, int jt = 60, double eps = 1e-6)
{
	//A为Hessenberg矩阵
	//jt 为最大迭代次数,eps为精度
	int n = A.n;
	int m, it, i, j, k, l;//m行数,it已经迭代次数
	double b, c, w, g, xy, p, q, r, x, s, e, f, z, y;
	double* u = new double[n];//保存实部
	double* v = new double[n];//保存虚部
	it = 0; m = n;
	while (m > 0)
	{
		l = m;
		while ((l > 1) && (fabs(A(l, l - 1)) > eps * fabs(A(l - 1, l - 1) + fabs(A(l, l))))) l--;
		if (l == m)
		{
			u[m - 1] = A(m, m);
			v[m - 1] = 0;
			m--; it = 0;
		}
		else if (l == m - 1)
		{
			b = -(A(m, m) + A(m - 1, m - 1));
			c = A(m, m) * A(m - 1, m - 1) - A(m, m - 1) * A(m - 1, m);
			w = b * b - 4.0 * c;
			y = sqrt(fabs(w));
			if (w > 0.0)
			{
				xy = 1.0;
				if (b < 0.0) xy = -1.0;
				u[m - 1] = (-b - xy * y) / 2;
				u[m - 2] = c / u[m - 1];
				v[m - 1] = 0;
				v[m - 2] = 0;
			}
			else
			{
				u[m - 2] = u[m - 1] = -b / 2.0;
				v[m - 1] = y / 2.0;
				v[m - 2] = -v[m - 1];
			}
			m = m - 2;
			it = 0;
		}
		else
		{
			if (it >= jt)
			{
				cout << "Failed,return NULL\n";
				return NULL;
			}
			it++;
			for (j = l + 2; j <= m; j++) A(j, j - 2) = 0;
			for (j = l + 3; j <= m; j++) A(j, j - 3) = 0;
			for (k = l; k <= m - 1; k++)
			{
				if (k != l)
				{
					p = A(k, k - 1);
					q = A(k + 1, k - 1);
					r = 0.0;
					if (k != m - 1) r = A(k + 2, k - 1);
				}
				else
				{
					x = A(m, m) + A(m - 1, m - 1);
					y = A(m, m) * A(m - 1, m - 1) - A(m, m - 1) * A(m - 1, m);
					p = A(l, l) * (A(l, l) - x) + A(l, l + 1) * A(l + 1, l) + y;
					q = A(l + 1, l) * (A(l, l) + A(l + 1, l + 1) - x);
					r = A(l + 1, l) * A(l + 2, l + 1);
				}
				if (fabs(p) + fabs(q) + fabs(r) != 0.0)
				{
					xy = 1.0;
					if (p < 0.0) xy = -1;
					s = xy * sqrt(p * p + q * q + r * r);
					if (k != l) A(k, k - 1) = -s;
					e = -q / s;
					f = -r / s;
					x = -p / s;
					y = -x - f * r / (p + s);
					g = e * r / (p + s);
					z = -x - e * q / (p + s);
					for (j = k; j <= m; j++)
					{
						p = x * A(k, j) + e * A(k + 1, j);
						q = e * A(k, j) + y * A(k + 1, j);
						r = f * A(k, j) + g * A(k + 1, j);
						if (k != m - 1)
						{
							p = p + f * A(k + 2, j);
							q = q + g * A(k + 2, j);
							r = r + z * A(k + 2, j);
							A(k + 2, j) = r;
						}
						A(k + 1, j) = q;
						A(k, j) = p;
					}
					j = k + 3;
					if (j >= m) j = m;
					for (i = l; i <= j; i++)
					{
						p = x * A(i, k) + e * A(i, k + 1);
						q = e * A(i, k) + y * A(i, k + 1);
						r = f * A(i, k) + g * A(i, k + 1);
						if (k != m - 1)
						{
							p = p + f * A(i, k + 2);
							q = q + g * A(i, k + 2);
							r = r + z * A(i, k + 2);
							A(i, k + 2) = r;
						}
						A(i, k) = p;
						A(i, k + 1) = q;
					}
				}
			}
		}
	}
	double* uv = new double[2 * n];
	for (i = 0; i < n; i++)
	{
		uv[i] = u[i];
		uv[i + n] = v[i];
	}//0~n-1个元素为特征值的实部,n~2n-1个元素为特征值的虚部, 特征值为uv[k]+uv[k+n]i
	return uv;
}


Complex* EigenValue(const Matrix<double>& A)
{
	double* p = HessenMatQR(HessenBerg_Transformation(A));
	Complex* c = new Complex[A.n];
	for (int i = 0; i < A.n; i++)
	{
		c[i] = Complex(p[i], p[i + A.n]);
	}
	if (p) delete p;//防止内存泄漏
	return c;
}