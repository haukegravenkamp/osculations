# osculations

**Compute dispersion curves of elastic waveguides utilizing uniform block-diagonalization**

These Matlab codes compute the generalized eigenvalues of a matrix function (aka matrix flow) $E(k)$, which corresponds to solving  
$E(k)\,\phi = \lambda\,M\,\phi$  
with eigenvector \phi, eigenvalue \lambda, constant matrix $M$, and $E(k)$ is a matrix function depending on one parameter $k$.
Specifically, we focus on matrix flows of the polynomial form  
$E(k)=k^2  E_0 - k\,E_1 + E_2$  
and assume the matrix flow to be Hermitian and $M$ nonsingular.  

Other matrix flows can be incorporated by changing the function `matrixFlow.m` as long as these conditions are met.

The implementation was done for the application of computing dispersion curves of elastic guided waves. When discretizing the cross-section of such a waveguide by finite elements, we obtain a matrix flow of the form shown above.

The strategy for computing such eigencurves is to first test whether the matrix flow can be uniformly decomposed, i.e., whether there is a similarity transformation that leads to the matrix flow being block-diagonal with the same block structure for any $k$. If such a decomposition exists, the eigencurves are computed separately for each block. This is then more efficient than solving the complete system at once. Furthermore, it is known that the eigencurves of each non-decomposable block do not cross. Hence, we do not require mode-tracing to determine which solutions belong to which curve. As a consequence, crossing-points are trivially distinguished from points where two curves approach each other very closely without crossing. The latter situation is known as *osculation, veering, avoided crossing, mode repulsion*, ...

# literature

The approach used here for this very application is described in detail in the following paper. The examples included in the repository reproduce the results in the paper:  
> [1] Gravenkamp, H., Plestenjak, B., & Kiefer, D. A.. Notes on osculations and mode tracing in semi-analytical waveguide modeling (under review). Ultrasonics.

The general concept of block-diagonalizing matrix flows in this manner was previously presented by Frank Uhlig, e.g. in  
> [2] Uhlig, F. (2022). On the unitary block-decomposability of 1-parameter matrix flows and static matrices. Numerical Algorithms, 89(2), 529–549. https://doi.org/10.1007/s11075-021-01124-7

The details on how the matrices are computed for the case of elastic waveguides are found in, e.g.,  
> [3] Gravenkamp, H., Song, C., & Prager, J. (2012). A numerical approach for the computation of dispersion relations for plate structures using the scaled boundary finite element method. Journal of Sound and Vibration, 331, 2543–2557. https://doi.org/10.1016/j.jsv.2012.01.029   
> [4] Gravenkamp, H., Man, H., Song, C., & Prager, J. (2013). The computation of dispersion relations for three-dimensional elastic waveguides using the Scaled Boundary Finite Element Method. Journal of Sound and Vibration, 332, 3756–3771. https://doi.org/10.1016/j.jsv.2013.02.007

# usage

Simply download the folder and run any of the `example*.m` files in Matlab.

Then you can try it with your own matrices or modify the matrix flow to some other function.  
If you encounter any issues of have an interesting application to discuss, do not hesitate to get in touch.



