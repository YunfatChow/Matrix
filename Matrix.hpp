#pragma once
#ifndef Matrix_H_
#define Matrix_H_
//Matrix.h
#include<iostream>
#include<fstream>
#include<cmath>
#include"Complex.hpp"

template <typename ElemType = double>
class Matrix
{
public:
	int n, m;
	ElemType** element;
	Matrix(int = 0, int = 0, int = 0);//第三个参数判断初始化类型，0为零矩阵，1为单位阵
	Matrix(const Matrix<ElemType>&);
	~Matrix();
	void Scan();
	void Print();
	void row_exchange(int, int);//行交换
	void row_add(int, int, ElemType);//给第i行加上第j行的k倍
	void row_multiply(int, ElemType);//行乘法
	void down_em(int = 0);//向下消元 
	void up_em(int = -2);//向上消元
	Matrix<ElemType> Augement(Matrix<ElemType>);//增广矩阵
	ElemType det();
	ElemType& operator()(int, int);
	Matrix<ElemType> operator+(const Matrix<ElemType>&) const;
	Matrix<ElemType> operator-(const Matrix<ElemType>&) const;
	Matrix<ElemType> operator*(const Matrix<ElemType>&) const;
	Matrix<ElemType> operator=(const Matrix<ElemType>&);
	Matrix<ElemType> transpose();
	bool NullRowCheck(int i);//检查第i行是否为零行
	ElemType trace();//求迹
};
template <typename ElemType>
std::istream& operator>>(std::istream&, Matrix<ElemType>&);

template <typename ElemType>
std::ostream& operator<<(std::ostream&, const Matrix<ElemType>&);

typedef Matrix<double> RealMatrix;
typedef Matrix<Complex> ComplexMatrix;

//QR_Decomposition.h
Matrix<double>* QR(const Matrix<double>&);//S为实矩阵
//--------------------------------------------------

//Matrix.cpp
template <typename ElemType>
Matrix<ElemType>::Matrix(int m, int n, int f) :m(m), n(n)
{
	if (m == 0 || n == 0)
	{
		element = NULL;
		return;
	}
	element = new ElemType * [m];
	for (int i = 0; i < m; i++) element[i] = new ElemType[n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			element[i][j] = 0;
		}
	}
	if (f == 1) for (int i = 0; i < m && i < n; i++) element[i][i] = 1;
}
template <typename ElemType>
Matrix<ElemType>::Matrix(const Matrix<ElemType>& A) :m(A.m), n(A.n)
{
	element = new ElemType * [m];
	for (int i = 0; i < m; i++) element[i] = new ElemType[n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			element[i][j] = A.element[i][j];
		}
	}
}
template <typename ElemType>
Matrix<ElemType>::~Matrix()
{
	for (int i = 0; i < m; i++) if (element[i] != NULL) delete[] element[i];
	if (element != NULL) delete[] element;
	element = NULL;
}
template <typename ElemType>
ElemType Matrix<ElemType>::det()
{
	if (n != m)
	{
		cout << "只有方阵才能求行列式\n";
		return 0;
	}
	Matrix B = *this;
	B.down_em();
	ElemType d = 1;
	for (int i = 0; i < n; i++) d = d * B.element[i][i];
	return d;
}
template <typename ElemType>
ElemType Matrix<ElemType>::trace()
{
	if (m != n)
	{
		cout << "只有方阵才有迹\n";
		return 0;
	}
	ElemType t = 0;
	for (int i = 0; i < n; i++) t = t + element[i][i];
	return t;
}
template <typename ElemType>
void Matrix<ElemType>::Scan()
{
	cin >> m >> n;
	element = new ElemType * [m];
	for (int i = 0; i < m; i++) element[i] = new ElemType[n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cin >> element[i][j];
		}
	}
}
template <typename ElemType>
void Matrix<ElemType>::Print()
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << element[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}
template <typename ElemType>
void Matrix<ElemType>::row_exchange(int i, int j)//交换第i，j行
{
	if (i<1 || i>m || j<1 || j>n)
	{
		cout << "Row exchange--Index out of range.\n";
		return;
	}
	i--;
	j--;
	for (int k = 0; k < n; k++)
	{
		std::swap(element[i][k], element[j][k]);
	}
}
template <typename ElemType>
void Matrix<ElemType>::row_add(int i, int j, ElemType k)//给第i行加上k倍的第j行
{
	if (i<1 || i>m || j<1 || j>n)
	{
		cout << "Row add--Index out of range.\n";
		return;
	}
	i--;
	j--;
	for (int p = 0; p < n; p++)
	{
		element[i][p] = element[i][p] + k * element[j][p];
	}
}
template <typename ElemType>
void Matrix<ElemType>::row_multiply(int i, ElemType k)//给第i行数乘k
{
	if (i<1 || i>m)
	{
		cout << "Row multiply--Index out of range.\n";
		return;
	}
	i--;
	for (int p = 0; p < n; p++)
	{
		element[i][p] = element[i][p] * k;
	}
}
template <typename ElemType>
bool Matrix<ElemType>::NullRowCheck(int i)
{
	if (i<1 || i>m)
	{
		cout << "NullRowCheck--Index out of range.\n";
		return 0;
	}
	i--;
	for (int j = 0; j < n; j++)
	{
		if (element[i][j]) return false;
	}
	return true;
}
template <typename ElemType>
Matrix<ElemType> Matrix<ElemType>::Augement(Matrix<ElemType> B)//返回[A|B]的增广矩阵
{
	if (m != B.m) return *this;//行数不相等不可增广
	Matrix<ElemType> C(m, B.n + n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			C.element[i][j] = element[i][j];
		}
		for (int j = n; j < n + B.n; j++)
		{
			C.element[i][j] = B.element[i][j - B.n];
		}
	}
	return C;
}
template <typename ElemType>
void Matrix<ElemType>::down_em(int r)
{
	if (r >= m - 1) return;
	if (element[r][r] == 0)
	{
		int p = r;
		while (p < m && element[p][r] == 0) p++;
		if (p >= m) return;
		row_exchange(r + 1, p + 1);
	}
	for (int i = r + 1; i < m; i++)//用第r行消去r+1至m行的第r个元素
	{
		ElemType k = element[i][r] / element[r][r];
		row_add(i + 1, r + 1, -k);
	}
	down_em(r + 1);
}
template <typename ElemType>
void Matrix<ElemType>::up_em(int r)//针对上三角矩阵 
{
	if (r == -2) r = m - 1;//默认参数，从第m行开始消元
	if (r <= 0) return;//r==0 消元至第一行
	while (element[r][r] == 0) r--;
	for (int i = 0; i <= r - 1; i++)//用第r行消去第1至r-1行的第r个元素
	{
		ElemType k = element[i][r] / element[r][r];
		row_add(i + 1, r + 1, -k);
	}
	up_em(r - 1);
}
template <typename ElemType>
Matrix<ElemType> inverse(Matrix<ElemType> A)
{
	if (A.m != A.n)
	{
		cout << "只有方阵才存在逆矩阵\n";
		return A;
	}
	Matrix<ElemType> B(A.m, A.n, 1);
	A = A.Augement(B);
	A.down_em();
	for (int i = 0; i < A.m; i++)
	{
		if (A.element[i][i] == 0)
		{
			cout << "矩阵A不可逆\n";
			return B;//返回值为单位阵
		}
	}
	A.up_em();
	for (int i = 0; i < B.n; i++)
	{
		for (int j = 0; j < B.n; j++)
		{
			B.element[i][j] = A.element[i][j + B.n] / A.element[i][i];//取出消元后的矩阵B并单位化A 
		}
	}
	return B;
}
template <typename ElemType>
Matrix<ElemType> Matrix<ElemType>::operator+(const Matrix<ElemType>& B) const
{
	Matrix<ElemType> C(m, n, 0);
	if (n != B.n || m != B.m)
	{
		cout << "矩阵不同型不可相加\n";
		return C;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			C.element[i][j] = element[i][j] + B.element[i][j];
		}
	}
	return C;
}

template <typename ElemType>
Matrix<ElemType> Matrix<ElemType>::operator-(const Matrix<ElemType>& B)const
{
	Matrix C(m, n, 0);
	if (n != B.n || m != B.m)
	{
		cout << "矩阵不同型不可相减\n";
		return C;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			C.element[i][j] = element[i][j] - B.element[i][j];
		}
	}
	return C;
}

template <typename ElemType>
Matrix<ElemType> Matrix<ElemType>::operator*(const Matrix<ElemType>& B) const
{
	Matrix<ElemType> C(m, B.n, 0);//C初始化为零矩阵 
	if (n != B.m)
	{
		std::cout << "该乘法不合法\n";
		return C;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < B.n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				C.element[i][j] = C.element[i][j] + element[i][k] * B.element[k][j];
			}
		}
	}
	return C;
}
template <typename ElemType>
Matrix<ElemType> Matrix<ElemType>::operator=(const Matrix<ElemType>& A)
{
	if (element != NULL)
	{
		for (int i = 0; i < m; i++) if (element[i] != NULL) delete[] element[i];
		if (element) delete element;
	}
	m = A.m;
	n = A.n;
	element = new ElemType * [m];
	for (int i = 0; i < m; i++) element[i] = new ElemType[n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			element[i][j] = A.element[i][j];
		}
	}
	return *this;
}
template <typename ElemType>
ElemType& Matrix<ElemType>::operator()(int i, int j)
{
	try
	{
		if (i<1 || i>m || j<1 || j>n)
		{
			throw "Error:Index(i,j) out of range.\n";
		}
	}
	catch (const char* s)
	{
		std::cout << s;
	}
	return element[i - 1][j - 1];
}
template <typename ElemType>
Matrix<ElemType> Matrix<ElemType>::transpose()
{
	Matrix<ElemType> X(n, m, 0);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			X.element[i][j] = element[j][i];
		}
	}
	return X;
}

template <typename ElemType>
std::istream& operator>>(std::istream& fin, Matrix<ElemType>& A)
{
	A.Scan();
	return fin;
}
template <typename ElemType>
std::ostream& operator<<(std::ostream& fout, const Matrix<ElemType>& A)
{
	for (int i = 0; i < A.m; i++)
	{
		for (int j = 0; j < A.n; j++)
		{
			cout << A.element[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";
	return fout;
}
template <typename ElemType>
void SolveEquation(Matrix<ElemType> A)
{
	A.down_em();
	A.up_em();
	for (int i = 0; i < A.m; i++)
	{
		if (A.element[i][i] != 0) A.row_multiply(i, 1 / A.element[i][i]);
	}
	int j = A.m;
	for (int i = 0; i < A.m && i < j; i++)//将零行移到最下方 
	{
		if (A.NullRowCheck(i))//第i行为零行 
		{
			A.row_exchange(i, j);
			j--;
		}
	}
	A.Print();
}
//显式实例化
template std::istream& operator>>(std::istream&, Matrix<double>&);
template std::ostream& operator<<(std::ostream&, const Matrix<double>&);

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
template class Matrix<Complex>;

//QR_Decomposition.cpp
Matrix<double>* QR(const Matrix<double>& S)//S为实矩阵
{
	int m = S.m;
	int n = S.n;
	Matrix<double> Q(m, m, 1);//初始化Q为单位阵
	Matrix<double> A(S);
	int i, j, k, nn, jj;
	double u, alpha, w, t;
	nn = n;
	if (m == n) nn--;
	for (k = 1; k <= nn; k++)
	{
		u = 0.0;
		for (i = k; i <= m; i++)
		{
			w = fabs(double(A(i, k)));
			if (w > u) u = w;
		}
		alpha = 0.0;
		for (i = k; i <= m; i++)
		{
			t = A(i, k) / u;
			alpha += t * t;
		}
		if (A(k, k) > 0.0) u = -u;
		alpha = u * sqrt(alpha);
		if (fabs(alpha) + 1.0 == 1.0)
		{
			printf("Fail\n");
			return NULL;
		}
		u = sqrt(2.0 * alpha * (alpha - A(k, k)));
		if ((u + 1.0) != 1.0)
		{
			A(k, k) = (A(k, k) - alpha) / u;
			for (i = k + 1; i <= m; i++)
			{
				A(i, k) /= u;
			}
		}
		for (j = 1; j <= m; j++)
		{
			t = 0;
			for (jj = k; jj <= m; jj++)
			{
				t += A(jj, k) * Q(jj, j);
			}
			for (i = k; i <= m; i++)
			{
				Q(i, j) -= 2 * t * A(i, k);
			}
		}
		for (j = k + 1; j <= n; j++)
		{
			t = 0;
			for (jj = k; jj <= m; jj++)
			{
				t += A(jj, k) * A(jj, j);
			}
			for (i = k; i <= m; i++)
			{
				A(i, j) -= 2 * t * A(i, k);
			}
		}
		A(k, k) = alpha;
		for (i = k + 1; i <= m; i++)
		{
			A(i, k) = 0;
		}
	}
	Q = Q.transpose();
	Matrix<double>* qr = new Matrix<double>[2];
	qr[0] = Q;
	qr[1] = A;
	return qr;
}
#endif