function y = heavisideReg(x, eps)

y = 0.5 * (1. + (2. / pi) * atan(x / eps));
