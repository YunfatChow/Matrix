# Matrix-QR

## Matrix.hpp

Matrix template

**RealMatrix** for `Matrix <double>`

**ComplexMatrix** for `Matrix <Complex>`

function:

- tranpose
- trace
- determinate
- SolveEquation
- istream >> & ostream <<

## Complex.hpp

Complex number

- fundamental **algebraic operations** + - * /

- **conjugate**
- equal (==) & not equal (!=)

- istream >> & ostream <<

## EigenValue.hpp

- HessenBerg_Transformation (transform a matrix to a **HessenBerg matrix**)
- HessenMatQR (find **eigenvalues** of a HessenBerg matrix using QR decomposition)
- EigenValue (**find a real matrix’s eigenvalues** and return as a complex array)