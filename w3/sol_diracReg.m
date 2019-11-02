function y = sol_diracReg(x, eps)
%Dirac function of x
%   sol_diracReg( x, epsilon ) Computes the derivative of the heaviside
%   function of x with respect to x. Regularized based on epsilon.

y = eps ./ (pi * (eps^2 + x.^2)); %TODO 19: Line to complete
