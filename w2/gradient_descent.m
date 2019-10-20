function [x, time] = gradient_descent(A, b, x, args)

    if isfield(args, 'tol')
        tolerance = args.tol;
    else
        tolerance = 1e-2;
    end
    
    % INIT 
    r = b - A*x; % direction oposite to the gradient
    err = norm(r);
    k = 1;
    tic
    while err > tolerance
        q = A*r;
        alpha = (transpose(r)*r) /(transpose(q)*r);
        
        r = b - A*x; 

        x = x + alpha*r;
        
        err = norm(r);
        
        fprintf('iter %03d, error: %.6f\n', k, err);
        k = k+1;
    end
    time = toc;

end
