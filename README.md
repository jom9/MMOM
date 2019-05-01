# Phys 490

The purpose of this repository is for documentation of numerical analysis of different materials and their optical properties

cauchyfitting.m
Reference: Fundemental of Optics by Jenkins, White pgs 479 - 485

This program takes in values of lamda and real index of refraction and fits them to the cauchy equation: n(λ)=Σai/λ^(2i)
This equation is an approximation of the following relations:

1.n^2-k^2=1+Σ(ai*λ^2)/[(λ^2-λi^2)^2+gi*λ^2]
2.2nk = ΣAi*(gi)^.5*λ^3/[(λ^2-λi^2)^2+gi*gi^2]

where k is the extiction factor and gi is the frictional force. These equations were derived using the helmhotlz model assumptions.

In this case, only the first three terms are taken for the fit. So the function will be in the form of n=a+b/λ^2+c/λ^4.The coeffiecents are 
stored in the coeff vector, and the values for n are stored in nfitted variables for later use.


createFit1.m
This is the helper function that creates the fit for cauchyfitting.m. It performs nonlinear regression using the method of least squares.
More info can be found here: [https://www.mathworks.com/help/curvefit/custom-nonlinear-models.html](url)

