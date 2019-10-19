function x = gauss_seidel(A, b, iters)
if nargin < 3
    iters = 10;
end
n = size(A, 2);
x = zeros(n, 1);
for k = 1:iters
    x_old = x;
    for i = 1:n
        x(i) = (b(i) - A(i,:)*x + A(i,i)*x(i)) / A(i,i);
    end
    if norm(x-x_old) < 1e-5
        break;
    end
    disp(norm(x-x_old));
end
