function x = gradient_descent(A, b, x, args)
    if isfield(args, 'tol')
        tol = args.tol;
    else
        tol = 1e-2;
    end
    
    r = b - A*x; %Opposite direction to the gradient
    err = norm(r);

    tic
    i = 1;
    while err > tol
        q = A*r;
        alpha = (transpose(r)*r) /(transpose(q)*r);
        r = b - A*x; 
        x = x + alpha*r;
        
        err = norm(r);
        fprintf('iter %03d, error: %.6f\n', i, err);
        
        i = i+1;
    end
    elapsedTime = toc;
    fprintf('Gradient Descent took %.3fs\n', elapsedTime);
end
