# Posit-Research
To run the solvers, first navigate to the Src folder and compile using the following command:  
`make NewBenchmarks`  
Then, you can run the solvers using the following command:  
```./NewBenchmarks```  
The command line will then display a series of questions to know which solver to run and with which settings to use. Follow the prompts and answer with "y" or "n" for yes-no questions and any text for other questsion. Then, the program will solve the matrix (or matrices) and save the final results in the plots folder. For Conjugate Gradient (CG), the program also saves the residual at each iteration. Details about each question are explained below:  
### Enter test ID
This questions asks which method to use for solving the linear system. Currently, only the CG method and Trisolve method are available. Enter 0 for CG and 1 for Trisolve.  
### Test all matrices
This question asks whether to solve all matrices in the MatrixMarket folder, or to solve only a specific matrix. Note that solving all matrices could take a very long time, and a failure in any matrix could cause the program to terminate and not solve any more matrices. There is also no guarantee on the order in which matrices are solved. However, any matrices that are solved before the program is terminated (for any reason) will be saved. More on that later.
### Enter mtx file name
This questions asks for the name of the file that contains the matrix to be solved. The matrix must be in the MatrixMarket folder. Please enter the entire file name, including the .mtx extension.
### Enter matrix name, and Enter matrix identifier
These two questions are used to identify the results stored in the plots folder. The final results are stored in a file called "CG\<identifier\>" or "Trisolve\<identifier\>". If multiple matrices are given the same identifier, their results will be stored in the same file, and their given names will differentiate between them. If the solver used is CG, the per-iteration residual will be stored in a file called "\<name\>\<identifier\>.csv".
### Scale?
Whether to use scaling or not.
### Quire?
Whether to use rounding deferred summation or not. **Currently, this doesn't work for Trisolve.**
### Stochastic?
Whether to use stochastic rounding or not. Note that if both Quire and Stochastic are chosen, stochastic will be ignored, and Quire will be used. **Currently, this doesn't work for Trisolve.**
### Cholesky? (Trisolve Only)
Whether to use the Cholesky version or not.
### Enter tolerance (CG Only)
This question asks about the tolerance. When the residual at an iteration is less than the tolerance, the method will count it as converged. This question accepts typing values in scientific notation. 1e-7 is a good default value to use.
### Enter iterations (CG Only)
This question asks about the maximum iterations for CG. If the number of iterations exceeds this number, the algorithm will stop, even if convergence hasn't been reached yet. 10000 is a good default value to use.