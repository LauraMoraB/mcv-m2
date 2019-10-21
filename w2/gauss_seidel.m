function x = gauss_seidel(A, b, x, args)
    if isfield(args, 'iters')
        iters = args.iters;
    else
        iters = 1000;
    end
    if isfield(args, 'tol')
        tol = args.tol;
    else
        tol = 1e-5;
    end
    if isfield(args, 'omega')
        omega = args.omega;
    else
        omega = 1;
    end
    
    n = size(A,1);
    tic
    for k = 1:iters
        x_old = x;
        for i = 1:n
            x(i) = (1-omega)*x(i) + (omega/A(i,i)) * (b(i) - A(i,:)*x + A(i,i)*x(i));
        end
        err = norm(x-x_old);
        fprintf('iter %03d, error: %.6f\n', k, err);
        if err < tol
            break;
        end
    end
    elapsedTime = toc;
    fprintf('Gauss-Seidel took %.3fs\n', elapsedTime);
end
