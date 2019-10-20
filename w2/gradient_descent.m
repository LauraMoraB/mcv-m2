function [x] = gradient_descent(A, b, x, args)

    if isfield(args, 'iters')
        iters = args.iters;
    else
        iters = 1000;
    end
    if isfield(args, 'tol')
        epsilon = args.tol;
    else
        epsilon = 1e-5;
    end
  
    % Start iteration
    for k = 1:iters
        %disp(k)
        r = b-A*x;
        e = transpose(r)*r;
        alpha = e /(transpose(r)*A*r);
        x = x + alpha*r;
        err = norm(r)
        
        fprintf('iter %03d, error: %.6f\n', k, err);
        % stopping criteria
        if  err < epsilon
            break
        end
        
    end


end
