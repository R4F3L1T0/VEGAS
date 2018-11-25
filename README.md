# VEGAS
### An adaptive Monte Carlo technique to numerically calculate definite integrals.


VEGAS is an algorithm that numerically calculate, with any desired accuracy, definite integrals in any number of dimensions. Despite the deterministic methods, it does not slow down drastically as the dimensionality of the integration domain rises. It's based on importance sampling, it samples points from the probability distribution described by the function |f|, so that the points are concentrated in the regions that make the largest contribution to the integral. The most important feature of the VEGAS algorithm is that it can be used when the probability distribution is not known, in fact the algorithm calculate the distribution on its own step-by-step. For this reason it belongs to the category of adaptive MonteCarlo techniques.  

In this repository you will find:

  - my VEGAS 1D implementation (**vegas.cpp**), the function to be 
  integrated and the integration domain are hard coded but you can change
  them as appropriate.
  
  - some plotting tools (**plotter.py**).
  
  - a very simple BASH (**vegas.sh**) executable that runs the VEGAS algortihm and consequently
  plot the results.
  
This picture shows the peculiar features of the VEGAS algorithm.  
The point inside the integration domain are pulled out of the probability law p(x), that is calculated steb-by-step (importance sampling), while a "crude" Monte Carlo would extract the points uniformly (simple sampling).  
That's key feature of the VEGAS algorithm that allows it to reach higher accuracy in a significantly lower time and to work properly in any number of dimentions.
