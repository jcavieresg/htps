# Htps 
### Efficient spatial interpolation using R

This repository contain a reproducible example of the article "Efficient estimation for a smoothing thin plate spline in a two-dimensional space". For this, do the following:

1. Download all the folders of this repository.
2. Copya and paste the folder $\textbf{htps}$ in your desktop. 
3. Download the Microsoft SDK (to use MPI in Windows): https://developer.microsoft.com/en-us/windows/downloads/windows-sdk/
4. Make sure link the $\texttt{c++}$ header library "htool" and the MPI folder to your $\texttt{Makevars}$ files correctly. For example, in the $\textbf{src}$ folder of the "htps" library you can find two $\texttt{Makevars}$ files: $\texttt{Makevars.win}$ and $\texttt{Makevars}$. So, for the $\texttt{Makevars.win}$:
</br>

* CXX_STD = CXX11
* PKG_CXXFLAGS = `$`(SHLIB_OPENMP_CXXFLAGS) -I/Users/Usuario/Desktop/htps/htool/include -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
* PKG_LIBS = `$`(SHLIB_OPENMP_CXXFLAGS) `$`(LAPACK_LIBS) `$`(BLAS_LIBS) `$`(FLIBS) -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi
</br>

and the same for $\texttt{Makevars}$ file. 
   
5. Compile the "htps" library using "Rcpp" of R following the "build_package.R" code. 
6. Run the "htps.R" code to reproduce the results.

**IMPORTANT**
The "htool" library used in this experiment is a slight adaptation of the original library written in $\texttt{c++}$ ([https://htool-documentation.readthedocs.io/en/latest/](https://github.com/htool-ddm/htool)). For a more comprehensive and original version of "htool" follow this link: https://htool-documentation.readthedocs.io/en/latest/


If you have any suggestions or would like to improve the way "htool" is linked with R, please feel free to send an email to: j.cavieres.g@gmail.com
