# Matrix-QR

## Matrix.hpp

Matrix template

**RealMatrix** for `Matrix <double>`

**ComplexMatrix** for `Matrix <Complex>`

function:

- tranpose
- trace
- rank
- determinate
- SolveEquation
- istream >> & ostream <<

## Complex.hpp

Complex number

- fundamental **algebraic operations** + - * /
- **conjugate**
- Norm()
- equal (==) & not equal (!=)
- sin() & cos()
- exp()
- pow(z,p)  or  z^p
- ntrt() ------**the nth root**
- log()
- sinh() & cosh()
- ostream <<  -------output form is “a+bi”
- istream >> --------inputform is “a+bi”

## EigenValue.hpp

- HessenBerg_Transformation (transform a matrix to a **HessenBerg matrix**)
- HessenMatQR (find **eigenvalues** of a HessenBerg matrix using QR decomposition)
- EigenValue (**find a real matrix’s eigenvalues** and return as a complex array)