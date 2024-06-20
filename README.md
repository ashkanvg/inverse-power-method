
# Inverse Power Method for Eigenvalue Calculation

This repository contains a Python implementation of the Inverse Power Method to calculate the eigenvalue of a given matrix. The code includes several helper functions for matrix operations, including forward and backward substitution, LU decomposition, Cholesky decomposition, and more.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
  - [lowerMatrix](#lowermatrix)
  - [upperMatrix](#uppermatrix)
  - [isSymmetric](#issymmetric)
  - [gauss](#gauss)
  - [cholesky](#cholesky)
  - [LU](#lu)
  - [solve](#solve)
  - [isDefinitePositive](#isdefinitepositive)
  - [generateSeries](#generateseries)
- [Example](#example)
- [Contributing](#contributing)

## Installation

Ensure you have Python and NumPy installed. You can install NumPy using pip:

```bash
pip install numpy
```

## Usage

To use the code, simply run the script and follow the prompts to enter the size of the matrix, the matrix elements, the value of rho, and the number of final steps.

```bash
python inverse_power_method.py
```

## Functions

### lowerMatrix
Performs forward substitution to solve the lower triangular matrix equation \( L \cdot y = b \).

### upperMatrix
Performs backward substitution to solve the upper triangular matrix equation \( U \cdot x = y \).

### isSymmetric
Checks if a given matrix is symmetric.

### gauss
Performs Gaussian elimination to decompose the matrix and return the transformed matrix and the pivot indices.

### cholesky
Performs Cholesky decomposition on a symmetric positive definite matrix.

### LU
Performs LU decomposition using Gaussian elimination.

### solve
Solves the system of linear equations given \( L \) and \( U \) matrices from LU decomposition.

### isDefinitePositive
Checks if a matrix is positive definite using Cholesky decomposition.

### generateSeries
Generates the series of vectors for the Inverse Power Method, printing intermediate vectors at each step.

## Example

Here is an example of input values:
```
Enter the size of matrix 'A': 2
Enter A[0][0]:
1
Enter A[0][1]:
2
Enter A[1][0]:
-1
Enter A[1][1]:
3
Enter the rho of the algorithm:(number) 1
Enter the final steps of the algorithm:(number) 5
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.
