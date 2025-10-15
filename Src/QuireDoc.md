# Quire Documentation
## What is quire?
Quire is a high-precision fixed-point accumulator used to prevent rounding errors during computations, especially in vector and matrix multiplications. In our case, we use it to implement rounding deferred summation. Rounding deferred summation is a method of summation that doesn't round any intermediate values, and instead defers rounding to the final answer, preventing rounding errors. Specifically, we use it in two functions: to compute the inner product of two vectors, and to compute the multiplication of a matrix by a vector.
## How it works
Ideally, the quire would be built-in and optimized for any summation automatically.
Right now, we simulate the behavior of the quire using a variable of high precision, specifically an `mpf_class`. Whenever we do a dot product of two vectors, we store the sum as an `mpf_class`, and we cast the operands to `mpf_class` before adding them to the sum.
Our inner product function looks something like this:
```
Initialize an mpf_class "sum" to zero
Loop through the elements of the two vectors and do the following:
    Upcast both operands to mpf_class
    Multiply the operands
    Add the product to the sum
Downcast sum back to the desired precision and return the result
```
Our matrix-vector multiplication algorithm works similarly.
For more details about how these two functions work, check the comments in Matrix.hh for innerProductQuire and matVecQuire.
## Bit width
The standard bit width for the quire is 16 * the current data type's bit width. For example, a float of size 32 bits would use a quire of size 512 bits. However, the bit width can be adjusted or optimized as needed.