#pragma once

// Complex.h
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
class Complex
{
public:
    double real, imaginary;
    Complex();
    Complex(int); //构造函数们
    Complex(float);
    Complex(double);
    Complex(int, int);
    Complex(int, float);
    Complex(int, double);
    Complex(float, int);
    Complex(float, float);
    Complex(float, double);
    Complex(double, int);
    Complex(double, float);
    Complex(double, double);

    Complex(const Complex &);
    ~Complex();
    double Norm() const;
    double arg() const;

    Complex operator-() const;
    Complex operator+(const Complex &) const;
    Complex operator-(const Complex &) const;
    Complex operator*(const Complex &) const;
    Complex operator/(const Complex &) const;
    Complex operator^(const int &) const;
    Complex operator^(const float &) const;
    Complex operator^(const double &) const;
    Complex &operator+=(const Complex &);
    Complex &operator-=(const Complex &);
    Complex &operator*=(const Complex &);
    Complex &operator/=(const Complex &);
    Complex &operator=(const Complex &);
    bool operator==(const int &) const;
    bool operator==(const float &) const;
    bool operator==(const double &) const;
    bool operator==(const Complex &) const;
    bool operator!=(const int &)const;
    bool operator!=(const float &)const;
    bool operator!=(const double &)const;
    bool operator!=(const Complex &) const;

    template <typename T>
    Complex operator+(const T &) const; // T模板重载int float double
    template <typename T>
    Complex operator-(const T &) const;
    template <typename T>
    Complex operator*(const T &) const;
    template <typename T>
    Complex operator/(const T &) const;
    template <typename T>
    Complex &operator+=(const T &);
    template <typename T>
    Complex &operator-=(const T &);
    template <typename T>
    Complex &operator*=(const T &);
    template <typename T>
    Complex &operator/=(const T &);
    template <typename T>
    Complex &operator=(const T &);
    template <typename T>
    bool operator==(const T &) const;
    template <typename T>
    bool operator!=(const T &) const;
    Complex Conjugate() const;
    operator double();
    friend ostream &operator<<(ostream &, const Complex &);
    friend istream &operator>>(istream &, Complex &);
};
template <typename T>
Complex operator+(const T &, Complex);
template <typename T>
Complex operator-(const T &, Complex);
template <typename T>
Complex operator*(const T &, Complex);
template <typename T>
Complex operator/(const T &, Complex);

Complex pow(const Complex &, int n);   //整数次幂
Complex pow(const Complex &, float p); //实数次幂
Complex pow(const Complex &, double p);
Complex *ntrt(const Complex &, int n); // n次方根
Complex exp(const Complex &);//复数指数
Complex log(const Complex &, float e = M_E);//复数的自然对数,默认底数为e
Complex sin(const Complex&);//正弦
Complex cos(const Complex&);//余弦
// Complex.cpp
Complex::Complex() : real(0), imaginary(0){};
Complex::Complex(int x) : real(x), imaginary(0){};   //构造函数们
Complex::Complex(float x) : real(x), imaginary(0){}; //后面不能带默认参数，否则会出现多义性
Complex::Complex(double x) : real(x), imaginary(0){};
Complex::Complex(int x, int y) : real(x), imaginary(y){};
Complex::Complex(int x, float y) : real(x), imaginary(y){};
Complex::Complex(int x, double y) : real(x), imaginary(y){};
Complex::Complex(float x, int y) : real(x), imaginary(y){};
Complex::Complex(float x, float y) : real(x), imaginary(y){};
Complex::Complex(float x, double y) : real(x), imaginary(y){};
Complex::Complex(double x, int y) : real(x), imaginary(y){};
Complex::Complex(double x, float y) : real(x), imaginary(y){};
Complex::Complex(double x, double y) : real(x), imaginary(y){};

Complex::Complex(const Complex &z) : real(z.real), imaginary(z.imaginary) {}

Complex::~Complex() {}

double Complex::Norm() const
{
    return sqrt(real * real + imaginary * imaginary);
}
double Complex::arg() const
{
    return atan2(imaginary, real);
}

Complex Complex::operator-() const
{
    return Complex(-real, -imaginary);
}

Complex Complex::operator+(const Complex &z) const
{
    return Complex(real + z.real, imaginary + z.imaginary);
}

Complex Complex::operator-(const Complex &z) const
{
    return Complex(real - z.real, imaginary - z.imaginary);
}

Complex Complex::operator*(const Complex &z) const
{
    return Complex(real * z.real - imaginary * z.imaginary, real * z.imaginary + imaginary * z.real);
}

Complex Complex::operator/(const Complex &z) const
{
    try
    {
        if (z.real == 0 && z.imaginary == 0)
        {
            throw "Divided by zero!\n";
        }
    }
    catch (const char *s)
    {
        cout << s;
    }
    double norm2 = z.real * z.real + z.imaginary * z.imaginary;
    return Complex((real * z.real + imaginary * z.imaginary) / norm2, (imaginary * z.real - real * z.imaginary) / norm2);
}
Complex Complex::operator^(const int &n) const
{
    return pow(*this, n);
}
Complex Complex::operator^(const float &p) const
{
    return pow(*this, p);
}
Complex Complex::operator^(const double &p) const
{
    return pow(*this, p);
}
Complex &Complex::operator+=(const Complex &z)
{
    real += z.real;
    imaginary += z.imaginary;
    return *this;
}

Complex &Complex::operator-=(const Complex &z)
{
    real -= z.real;
    imaginary -= z.imaginary;
    return *this;
}

Complex &Complex::operator*=(const Complex &z)
{
    double x = real * z.real - imaginary * z.imaginary;
    double y = real * z.imaginary + imaginary * z.real;
    real = x;
    imaginary = y;
    return *this;
}

Complex &Complex::operator/=(const Complex &z)
{
    try
    {
        if (z.real == 0 && z.imaginary == 0)
        {
            throw "Divided by zero!\n";
        }
    }
    catch (const char *s)
    {
        std::cout << s;
    }
    double norm2 = z.real * z.real + z.imaginary * z.imaginary;
    double x = (real * z.real + imaginary * z.imaginary) / norm2;
    double y = (imaginary * z.real - real * z.imaginary) / norm2;
    real = x;
    imaginary = y;
    return *this;
}

Complex &Complex::operator=(const Complex &z)
{
    real = z.real;
    imaginary = z.imaginary;
    return *this;
}
bool Complex::operator==(const int &n) const
{
    if (imaginary != 0)
        return false;
    return real == n;
}
bool Complex::operator==(const float &n) const
{
    if (imaginary != 0)
        return false;
    return real == n;
}
bool Complex::operator==(const double &n) const
{
    if (imaginary != 0)
        return false;
    return real == n;
}
bool Complex::operator==(const Complex &z) const
{
    return (real == z.real) && (imaginary == z.imaginary);
}
bool Complex::operator!=(const int &n)const
{
    if(imaginary != 0) return true;
    return real != n;
}
bool Complex::operator!=(const float &n)const
{
    if(imaginary != 0) return true;
    return real != n;
}
bool Complex::operator!=(const double &n)const
{
    if(imaginary != 0) return true;
    return real != n;
}
bool Complex::operator!=(const Complex &z) const
{
    return (real != z.real) || (imaginary != z.imaginary);
}

template <typename T>
Complex Complex::operator+(const T &x) const
{
    return Complex(real + x, imaginary);
}

template <typename T>
Complex Complex::operator-(const T &x) const
{
    return Complex(real - x, imaginary);
}

template <typename T>
Complex Complex::operator*(const T &x) const
{
    return Complex(real * x, imaginary * x);
}

template <typename T>
Complex Complex::operator/(const T &x) const
{
    try
    {
        if (x == 0)
        {
            throw "Divided by zero!\n";
        }
    }
    catch (const char *s)
    {
        std::cout << s;
    }
    return Complex(real / x, imaginary / x);
}

template <typename T>
Complex &Complex::operator+=(const T &x)
{
    real += x;
    return *this;
}

template <typename T>
Complex &Complex::operator-=(const T &x)
{
    real -= x;
    return *this;
}

template <typename T>
Complex &Complex::operator*=(const T &x)
{
    real *= x;
    imaginary *= x;
    return *this;
}

template <typename T>
Complex &Complex::operator/=(const T &x)
{
    real /= x;
    imaginary /= x;
    return *this;
}

template <typename T>
Complex &Complex::operator=(const T &x)
{
    real = x;
    return *this;
}

template <typename T>
bool Complex::operator==(const T &x) const //判断实数
{
    return (imaginary == 0) && (real == x);
}

template <typename T>
bool Complex::operator!=(const T &x) const //判断不等于实数
{
    return (imaginary != 0) || (real != x);
}

template <typename T>
Complex operator+(const T &x, Complex z)
{
    return Complex(x + z.real, z.imaginary);
}

template <typename T>
Complex operator-(const T &x, Complex z)
{
    return Complex(x - z.real, z.imaginary);
}

template <typename T>
Complex operator*(const T &x, Complex z)
{
    return Complex(x * z.real, x * z.imaginary);
}

template <typename T>
Complex operator/(const T &x, Complex z)
{
    return Complex(x) / z;
}

Complex Complex::Conjugate() const
{
    return Complex(real, -imaginary);
}
Complex::operator double()
{
    return Norm();
}
Complex pow(const Complex &z, int n)
{
    double r = z.Norm();
    double arg = z.arg();
    if (r + 1.0 == 1.0)
        return Complex(0.0);
    r = pow(r, n);
    return Complex(r * cos(n * arg), r * sin(n * arg));
}
Complex pow(const Complex &z, float p)
{
    double r = z.Norm();
    double arg = z.arg();
    if (r + 1.0 == 1.0)
        return Complex(0.0);
    r = pow(r, p);
    return Complex(r * cos(p * arg), r * sin(p * arg));
}
Complex pow(const Complex &z, double p)
{
    double r = z.Norm();
    double arg = z.arg();
    if (r + 1.0 == 1.0)
        return Complex(0.0);
    r = pow(r, p);
    return Complex(r * cos(p * arg), r * sin(p * arg));
}
Complex *ntrt(const Complex &z, int n)
{
    double r, arg, t;
    r = z.Norm();
    arg = z.arg();
    if (r + 1.0 == 1.0)
    {
        Complex *p = new Complex[1];
        p[0] = 0;
        return p;
    }
    Complex *p = new Complex[n];
    for (int i = 0; i < n; i++)
    {
        t = (2.0 * i * M_PI + arg) / n;
        p[i] = Complex(r * cos(t), r * sin(t));
    }
    return p;
}
Complex exp(const Complex &z)
{
    double r = exp(z.real);
    return Complex(r * cos(z.imaginary), r * sin(z.imaginary));
}
Complex log(const Complex &z, float e)
{
    return Complex(log(z.Norm()) / log(e), z.arg() / log(e));
}
Complex sin(const Complex& z)
{
    return Complex(sin(z.real) * cosh(z.imaginary), cos(z.real) * sinh(z.imaginary));
}
Complex cos(const Complex& z)
{
    return Complex(cos(z.real) * cosh(z.imaginary), -sin(z.real) * sinh(z.imaginary));
}
ostream &operator<<(ostream &os, const Complex &z) //重载<<输出
{
    //double eps = 1e-14;
    if (z.real == 0 && z.imaginary == 0)
    {
        os << "0";
        return os;
    }
    if (z.real != 0)
    {
        os << z.real;
    }
    if (z.imaginary != 0)
    {
        if (z.imaginary > 0 && z.real != 0)
            os << "+";
        if (z.imaginary == -1)
            os << "-";
        else if (z.imaginary != 1)
            os << z.imaginary;
        os << "i";
    }
    return os;
}
istream &operator>>(istream &is, Complex &z) //重载>>输入，要求格式为 a+bi (a可以为0，b可以为1)
{
    string s;
    char c;
    c = is.get();
    while (c == '\n')
        c = is.get();
    while (c != ' ' && c != '\n')
    {
        s += c;
        c = is.get();
    }
    int pos = int(s.find('i')); //判断i是否存在
    if (pos == s.npos)          //不存在i,即为纯实数
    {
        z.real = stod(s);
        z.imaginary = 0;
        return is;
    }
    //存在i，输入中存在虚部 (i必须写在虚部的结尾)
    int sgnpos = 0; //找到虚部前符号位置
    int sgn = 0;
    for (int i = 1; i < s.size(); i++)
    {
        if (s[i] == '+')
        {
            sgnpos = i;
            sgn = 1;
        }
        if (s[i] == '-')
        {
            sgnpos = i;
            sgn = -1;
        }
    }
    if (sgnpos == 0) //第二位起没有发现符号，则为纯虚数拷
    {
        if (s[0] == '-') //首位为负号
        {
            if (s[1] == 'i') //输入了-i
            {
                z.imaginary = -1;
            }
            else
                z.imaginary = -stod(s.substr(1)); //从第1位开始获取数字
        }
        else if (s[0] == 'i') //单个字符i 即只输入了i
        {
            z.imaginary = 1;
        }
        else
        {
            z.imaginary = stod(s);
        }
        z.real = 0; //实部为0
    }
    else
    {
        z.real = stod(s.substr(0, sgnpos));
        if (s[sgnpos + 1] == 'i') //符号位后是i，则虚部系数为1，乘上符号
        {
            z.imaginary = sgn;
        }
        else
            z.imaginary = sgn * stod(s.substr(sgnpos + 1));
    }
    return is;
}

template Complex Complex::operator+(const int &) const; //显式实例化 T=int
template Complex Complex::operator-(const int &) const;
template Complex Complex::operator*(const int &) const;
template Complex Complex::operator/(const int &) const;
template Complex &Complex::operator+=(const int &);
template Complex &Complex::operator-=(const int &);
template Complex &Complex::operator*=(const int &);
template Complex &Complex::operator/=(const int &);
template Complex &Complex::operator=(const int &);
template bool Complex::operator==(const int &) const;
template bool Complex::operator!=(const int &) const;

template Complex operator+(const int &, Complex);
template Complex operator-(const int &, Complex);
template Complex operator*(const int &, Complex);
template Complex operator/(const int &, Complex);

template Complex Complex::operator+(const float &) const; //显式实例化 T=float
template Complex Complex::operator-(const float &) const;
template Complex Complex::operator*(const float &) const;
template Complex Complex::operator/(const float &) const;
template Complex &Complex::operator+=(const float &);
template Complex &Complex::operator-=(const float &);
template Complex &Complex::operator*=(const float &);
template Complex &Complex::operator/=(const float &);
template Complex &Complex::operator=(const float &);
template bool Complex::operator==(const float &) const;
template bool Complex::operator!=(const float &) const;
template Complex operator+(const float &, Complex);
template Complex operator-(const float &, Complex);
template Complex operator*(const float &, Complex);
template Complex operator/(const float &, Complex);

template Complex Complex::operator+(const double &) const; //显式实例化 T=double
template Complex Complex::operator-(const double &) const;
template Complex Complex::operator*(const double &) const;
template Complex Complex::operator/(const double &) const;
template Complex &Complex::operator+=(const double &);
template Complex &Complex::operator-=(const double &);
template Complex &Complex::operator*=(const double &);
template Complex &Complex::operator/=(const double &);
template Complex &Complex::operator=(const double &);
template bool Complex::operator==(const double &) const;
template bool Complex::operator!=(const double &) const;

template Complex operator+(const double &, Complex); //显式实例化
template Complex operator-(const double &, Complex);
template Complex operator*(const double &, Complex);
template Complex operator/(const double &, Complex);

