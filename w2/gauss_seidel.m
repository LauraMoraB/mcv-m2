function x = gauss_seidel(A, b, x, args)
    if isfield(args, 'iters')
        iters = args.iters;
    else
        iters = 100;
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

    for k = 1:iters
        x_old = x;
        for i = 1:n
            x(i) = (1-omega)*x(i) + (omega/A(i,i)) * (b(i) - A(i,:)*x + A(i,i)*x(i));
        end
        r = norm(x-x_old);
        if r < tol
            break;
        end
        disp(r);
    end
end
