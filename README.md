# Htps 
### Efficient spatial interpolation using R

This repository contain a reproducible example of the article "Efficient estimation for a smoothing thin plate spline in a two-dimensional space". For this, do the following:

1. Download all the folders of this this repository.
2. Copya and paste the folder $\textbf{htps}$ in your desktop. 
3. Download the Microsoft SDK (to use MPI in Windows): https://developer.microsoft.com/en-us/windows/downloads/windows-sdk/
4. Make sure link the $\texttt{c++}$ header library "htool" and the MPI folder to your $\texttt{Makevars}$ files correctly. For example, in the $\textbf{src}$ folder of the "htps" library you can find two $\texttt{Makevars}$ files: $\texttt{Makevars.win}$ and $\texttt{Makevars}$. So, for the $\texttt{Makevars.win}$:
   
* PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) -I/Users/Usuario/Desktop/htps/htool/include -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
* PKG_LIBS = `$`(SHLIB_OPENMP_CXXFLAGS) `$`(LAPACK_LIBS) `$`(BLAS_LIBS) `$`(FLIBS) -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi
</br>

and the same for $\texttt{Makevars}$. 
   
5. Compile the "htps" library using "Rcpp" of R.
6. Run the "htps.R" code that you can find in the $\textbf{r_codes}" of this repository to reproduce the results. 
