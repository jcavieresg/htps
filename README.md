# Htps 
### Efficient spatial interpolation using R

This repository contain a reproducible example of the article "Efficient estimation for a smoothing thin plate spline in a two-dimensional space". For this, do the following:

1. Download all the folders of this repository.
2. Copy a and paste the folder $\textbf{htps}$ in your desktop.
3. Download MPI for Windows: https://www.microsoft.com/en-us/download/details.aspx?id=105289
4. Download the Microsoft SDK (to link the header libraries): https://developer.microsoft.com/en-us/windows/downloads/windows-sdk/
5. Ensure that you correctly link the "htool.hpp" $\texttt{c++}$ header library and the MPI folder in your $\texttt{Makevars}$ file. So, for the $\texttt{Makevars.win}$:
</br>

* CXX_STD = CXX11
* PKG_CXXFLAGS = `$`(SHLIB_OPENMP_CXXFLAGS) -Iinclude -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
* PKG_LIBS = `$`(SHLIB_OPENMP_CXXFLAGS) `$`(LAPACK_LIBS) `$`(BLAS_LIBS) `$`(FLIBS) -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi
</br>

6. Compile the $\texttt{htps.cpp}$ which is in the $\textbf{htps}$ folder using "Rcpp" by following the "build_package.R" code. Here you also have to replace "Usuario" with the username you are using on your PC or laptop.
7. Run the "htps.R" code to reproduce the results.

**IMPORTANT**
The "htool" library used in this experiment is a slight adaptation of the original library written in $\texttt{c++}$ ([https://htool-documentation.readthedocs.io/en/latest/](https://github.com/htool-ddm/htool)). For a more comprehensive and original version of "htool" follow this link: https://htool-documentation.readthedocs.io/en/latest/


If you have any suggestions or would like to improve the way "htool" is linked with R, please feel free to send an email to: j.cavieres.g@gmail.com

![plot_hmat](https://github.com/jcavieresg/htps/assets/55828236/b9c6bf46-9c62-44d8-8377-eb209324323b)


