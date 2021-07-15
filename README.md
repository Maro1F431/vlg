# C++ implementation of the BOUNDINGECCENTRICITIES algorithm (for very large graphs) using IGRAPH

## AUTHORS:

Mathieu TAMMARO
Paul ADAM

# Build the project

We provide you a simple makefile</br>

- If you want to debug the project (use gdb) just type:</br>
	**make debug**
- If just want to run the project a simple **make** will suffise</br>

We also provide a make clean in order to erase the executable file

# How to use it

To use our program you just need a .txt file containing a list of eges,
if you have doubt the example `facebook_combined.txt` is provided. Here is
an example of how to use our solution:</br>
-	**./a.out path/to/file.txt**

Furthermore you can use the --print optin to print the eccentricity list as such:\
	**./a.out --print path/to/file.txt**
