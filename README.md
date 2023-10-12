# Htps 
### Efficient spatial interpolation using R

This repository contain a reproducible example of the article "Efficient estimation for a smoothing thin plate spline in a two-dimensional space". For this, do the following:

1. Download the "htps" library from this repository.
2. Download the Microsoft SDK for MPI (Windows): https://developer.microsoft.com/en-us/windows/downloads/windows-sdk/
3. Make sure link the $\texttt{c++}$ header library "htool" to your $\texttt{Makevars}$ files and to the MPI folder correctly. For example, in the "src" folder of the "htps" library you can find two $\texttt{Makevars}$ files: "$\texttt{Makevars.win}$ and $\texttt{Makevars}$. 
4. Compile the "htps" library using "Rcpp" of R. 
