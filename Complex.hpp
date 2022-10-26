#pragma once
#ifndef Complex_H_
#define Complex_H_ 
//Complex.h
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
using namespace std;
class Complex{
public:
    double real,imaginary;
    Complex();
    Complex(int);//���캯����
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

    
    Complex(const Complex&);
    ~Complex();
    double Norm();
    Complex operator-()const;
    Complex operator+(const Complex&)const;
    Complex operator-(const Complex&)const;
    Complex operator*(const Complex&)const;
    Complex operator/(const Complex&)const;
    Complex& operator+=(const Complex&);
    Complex& operator-=(const Complex&);
    Complex& operator*=(const Complex&);
    Complex& operator/=(const Complex&);
    Complex& operator=(const Complex&);
    bool operator==(const Complex&)const;
    bool operator!=(const Complex&)const;

    template <typename T> Complex operator+(const T&)const;//Tģ������int float double
    template <typename T> Complex operator-(const T&)const;
    template <typename T> Complex operator*(const T&)const;
    template <typename T> Complex operator/(const T&)const;
    template <typename T> Complex& operator+=(const T&);
    template <typename T> Complex& operator-=(const T&);
    template <typename T> Complex& operator*=(const T&);
    template <typename T> Complex& operator/=(const T&);
    template <typename T> Complex& operator=(const T&);
    template <typename T> bool operator==(const T&)const;
    template <typename T> bool operator!=(const T&)const;
    Complex Conjugate()const;
    operator double();
    friend ostream& operator<<(ostream&, const Complex&);
    friend istream& operator>>(istream&, Complex&);
};
template <typename T> Complex operator+(const T&, Complex);
template <typename T> Complex operator-(const T&, Complex);
template <typename T> Complex operator*(const T&, Complex);
template <typename T> Complex operator/(const T&, Complex);

//Complex.cpp
Complex::Complex() :real(0), imaginary(0) {};
Complex::Complex(int x) :real(x), imaginary(0) {};//���캯����
Complex::Complex(float x) :real(x), imaginary(0) {};//���治�ܴ�Ĭ�ϲ������������ֶ�����
Complex::Complex(double x) :real(x), imaginary(0) {};
Complex::Complex(int x , int y ) :real(x), imaginary(y) {};
Complex::Complex(int x , float y ) :real(x), imaginary(y) {};
Complex::Complex(int x , double y ) :real(x), imaginary(y) {};
Complex::Complex(float x , int y ) :real(x), imaginary(y) {};
Complex::Complex(float x , float y ) :real(x), imaginary(y) {};
Complex::Complex(float x , double y ) :real(x), imaginary(y) {};
Complex::Complex(double x , int y ) :real(x), imaginary(y) {};
Complex::Complex(double x , float y ) :real(x), imaginary(y) {};
Complex::Complex(double x , double y ) :real(x), imaginary(y) {};

Complex::Complex(const Complex& z):real(z.real),imaginary(z.imaginary){}

Complex::~Complex(){}

double Complex::Norm()
{
    return sqrt(real*real+imaginary*imaginary);
}

Complex Complex::operator-()const
{
    return Complex(-real, -imaginary);
}

Complex Complex::operator+(const Complex& z)const
{
    return Complex(real+z.real,imaginary+z.imaginary);
}

Complex Complex::operator-(const Complex& z)const
{
    return Complex(real-z.real,imaginary-z.imaginary);
}

Complex Complex::operator*(const Complex& z)const
{
    return Complex(real*z.real-imaginary*z.imaginary,real*z.imaginary+imaginary*z.real);
}

Complex Complex::operator/(const Complex& z)const
{
    try
    {
        if (z.real == 0 && z.imaginary == 0)
        {
            throw "Divided by zero!\n";
        }
    }
    catch (const char* s)
    {
        cout << s;
    }
    double norm2 = z.real * z.real + z.imaginary * z.imaginary;
    return Complex((real * z.real + imaginary * z.imaginary) / norm2, (imaginary * z.real - real * z.imaginary) / norm2);
}

Complex& Complex::operator+=(const Complex& z)
{
    real += z.real;
    imaginary += z.imaginary;
    return *this;
}

Complex& Complex::operator-=(const Complex& z)
{
    real -= z.real;
    imaginary -= z.imaginary;
    return *this;
}

Complex& Complex::operator*=(const Complex& z)
{
    double x = real * z.real - imaginary * z.imaginary;
    double y = real * z.imaginary + imaginary * z.real;
    real = x;
    imaginary = y;
    return *this;
}

Complex& Complex::operator/=(const Complex& z)
{
    try
    {
        if (z.real == 0 && z.imaginary == 0)
        {
            throw "Divided by zero!\n";
        }
    }
    catch (const char* s)
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

Complex& Complex::operator=(const Complex& z)
{
    real = z.real;
    imaginary = z.imaginary;
    return *this;
}

bool Complex::operator==(const Complex& z)const
{
    return (real == z.real) && (imaginary == z.imaginary);
}
bool Complex::operator!=(const Complex& z)const
{
    return (real != z.real) || (imaginary != z.imaginary);
}

template <typename T>
Complex Complex::operator+(const T& x)const
{
    return Complex(real + x, imaginary);
}

template <typename T>
Complex Complex::operator-(const T& x)const
{
    return Complex(real - x, imaginary);
}

template <typename T>
Complex Complex::operator*(const T& x)const
{
    return Complex(real * x, imaginary * x);
}

template <typename T>
Complex Complex::operator/(const T& x)const
{
    try
    {
        if (x== 0)
        {
            throw "Divided by zero!\n";
        }
    }
    catch (const char* s)
    {
        std::cout << s;
    }
    return Complex(real / x, imaginary / x);
}

template <typename T>
Complex& Complex::operator+=(const T& x)
{
    real += x;
    return *this;
}

template <typename T>
Complex& Complex::operator-=(const T& x)
{
    real -= x;
    return *this;
}

template <typename T>
Complex& Complex::operator*=(const T& x)
{
    real *= x;
    imaginary *= x;
    return *this;
}

template <typename T>
Complex& Complex::operator/=(const T& x)
{
    real /= x;
    imaginary /= x;
    return *this;
}

template <typename T>
Complex& Complex::operator=(const T& x)
{
    real = x;
    return *this;
}

template <typename T>
bool Complex::operator==(const T& x)const//�ж�ʵ��
{
    return (imaginary == 0) && (real == x);
}

template <typename T>
bool Complex::operator!=(const T& x)const//�жϲ�����ʵ��
{
    return (imaginary != 0) || (real != x);
}

template <typename T>
Complex operator+(const T& x, Complex z)
{
    return Complex(x + z.real, z.imaginary);
}

template <typename T>
Complex operator-(const T& x, Complex z)
{
    return Complex(x - z.real, z.imaginary);
}

template <typename T>
Complex operator*(const T& x, Complex z)
{
    return Complex(x * z.real, x * z.imaginary);
}

template <typename T>
Complex operator/(const T& x, Complex z)
{
    return Complex(x) / z;
}

Complex Complex::Conjugate()const
{
    return Complex(real, -imaginary);
}
Complex::operator double()
{
    return Norm();
}

ostream& operator<<(ostream& os, const Complex& z)//����<<���
{
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
        if (z.imaginary > 0 && z.real != 0) os << "+";
        if (z.imaginary == -1) os << "-";
        else if (z.imaginary != 1) os << z.imaginary;
        os << "i";
    }
    return os;
}
istream& operator>>(istream& is, Complex& z)//����>>���룬Ҫ���ʽΪ a+bi (a����Ϊ0��b����Ϊ1)
{
    string s;
    char c;
    c = is.get();
    while (c == '\n') c = is.get();
    while (c != ' ' && c != '\n')
    {
        s += c;
        c = is.get();
    }
    int pos = int(s.find('i'));//�ж�i�Ƿ����
    if (pos == s.npos)//������i,��Ϊ��ʵ��
    {
        z.real = stod(s);
        z.imaginary = 0;
        return is;
    }
    //����i�������д����鲿 (i����д���鲿�Ľ�β)
    int sgnpos = 0;//�ҵ��鲿ǰ����λ��
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
    if (sgnpos == 0)//�ڶ�λ��û�з��ַ��ţ���Ϊ��������
    {
        if (s[0] == '-')//��λΪ����
        {
            if (s[1] == 'i')//������-i
            {
                z.imaginary = -1;
            }
            else z.imaginary = -stod(s.substr(1));//�ӵ�1λ��ʼ��ȡ����
        }
        else if (s[0] == 'i')//�����ַ�i ��ֻ������i
        {
            z.imaginary = 1;
        }
        else
        {
            z.imaginary = stod(s);
        }
        z.real = 0;//ʵ��Ϊ0
    }
    else
    {
        z.real = stod(s.substr(0, sgnpos));
        if (s[sgnpos + 1] == 'i')//����λ����i�����鲿ϵ��Ϊ1�����Ϸ���
        {
            z.imaginary = sgn;
        }
        else z.imaginary = sgn * stod(s.substr(sgnpos + 1));
    }
    return is;
}

template Complex Complex::operator+(const int&)const;//��ʽʵ���� T=int
template Complex Complex::operator-(const int&)const;
template Complex Complex::operator*(const int&)const;
template Complex Complex::operator/(const int&)const;
template Complex& Complex::operator+=(const int&);
template Complex& Complex::operator-=(const int&);
template Complex& Complex::operator*=(const int&);
template Complex& Complex::operator/=(const int&);
template Complex& Complex::operator=(const int&);
template bool Complex::operator==(const int&)const;
template bool Complex::operator!=(const int&)const;

template Complex operator+(const int&, Complex);
template Complex operator-(const int&, Complex);
template Complex operator*(const int&, Complex);
template Complex operator/(const int&, Complex);

template Complex Complex::operator+(const float&)const;//��ʽʵ���� T=float
template Complex Complex::operator-(const float&)const;
template Complex Complex::operator*(const float&)const;
template Complex Complex::operator/(const float&)const;
template Complex& Complex::operator+=(const float&);
template Complex& Complex::operator-=(const float&);
template Complex& Complex::operator*=(const float&);
template Complex& Complex::operator/=(const float&);
template Complex& Complex::operator=(const float&);
template bool Complex::operator==(const float&)const;
template bool Complex::operator!=(const float&)const;
template Complex operator+(const float&, Complex);
template Complex operator-(const float&, Complex);
template Complex operator*(const float&, Complex);
template Complex operator/(const float&, Complex);

template Complex Complex::operator+(const double&)const;//��ʽʵ���� T=double
template Complex Complex::operator-(const double&)const;
template Complex Complex::operator*(const double&)const;
template Complex Complex::operator/(const double&)const;
template Complex& Complex::operator+=(const double&);
template Complex& Complex::operator-=(const double&);
template Complex& Complex::operator*=(const double&);
template Complex& Complex::operator/=(const double&);
template Complex& Complex::operator=(const double&);
template bool Complex::operator==(const double&)const;
template bool Complex::operator!=(const double&)const;

template Complex operator+(const double&, Complex);//��ʽʵ����
template Complex operator-(const double&, Complex);
template Complex operator*(const double&, Complex);
template Complex operator/(const double&, Complex);

#endif
