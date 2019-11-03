function [seg] = sol_ChanVeseIpol_GDExp(I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni)
%Implementation of the Chan-Vese segmentation following the explicit
%gradient descent in the paper of Pascal Getreur "Chan-Vese Segmentation".
%It is the equation 19 from that paper

%I     : Gray color image to segment
%phi_0 : Initial phi
%mu    : mu lenght parameter (regularizer term)
%nu    : nu area parameter (regularizer term)
%eta   : epsilon for the total variation regularization
%lambda1, lambda2: data fidelity parameters
%tol   : tolerance for the sopping criterium
% epHeaviside: epsilon for the regularized heaviside. 
% dt     : time step
%iterMax : MAximum number of iterations
%reIni   : Iterations for reinitialization. 0 means no reinitializacion

[ni,nj]=size(I);
hi=1;
hj=1;

phi=phi_0;
for nIter = 1:iterMax
    phi_old=phi;
    
    %Fixed phi, Minimization w.r.t c1 and c2 (constant estimation)
    %H = heavisideReg(phi, epHeaviside);
    H = (phi > 0);
    c1 = sum(I(:) .* H(:)) / sum(H(:)); %TODO 1: Line to complete
    c2 = sum(I(:) .* ~H(:)) / sum(~H(:)); %TODO 2: Line to complete

    P = zeros(ni+2,nj+2);
    P(2:end-1,2:end-1) = phi;

    %Boundary conditions
    P(1,:)   = P(2,:); %TODO 3: Line to complete
    P(end,:) = P(end-1,:); %TODO 4: Line to complete

    P(:,1)   = P(:,2); %TODO 5: Line to complete
    P(:,end) = P(:,end-1); %TODO 6: Line to complete
    
    %Regularized Dirac's Delta computation
    delta_phi = sol_diracReg(phi, epHeaviside);   %notice delta_phi=H'(phi)	
    
    %derivatives estimation
    %i direction, forward finite differences
    phi_iFwd  = P(3:end,2:end-1) - P(2:end-1,2:end-1); %TODO 7: Line to complete
    phi_iBwd  = P(2:end-1,2:end-1) - P(1:end-2,2:end-1); %TODO 8: Line to complete
    
    %j direction, forward finitie differences
    phi_jFwd  = P(2:end-1,3:end) - P(2:end-1,2:end-1); %TODO 9: Line to complete
    phi_jBwd  = P(2:end-1,2:end-1) - P(2:end-1,1:end-2); %TODO 10: Line to complete
    
    %centered finite diferences
    phi_icent   = (P(3:end,2:end-1) - P(1:end-2,2:end-1)) / 2; %TODO 11: Line to complete
    phi_jcent   = (P(2:end-1,3:end) - P(2:end-1,1:end-2)) / 2; %TODO 12: Line to complete
    
    %A and B estimation (A y B from the Pascal Getreuer's IPOL paper "Chan
    %Vese segmentation
    A = mu ./ (sqrt(eta^2 + phi_iFwd.^2 + phi_jcent.^2)); %TODO 13: Line to complete
    B = mu ./ (sqrt(eta^2 + phi_icent.^2 + phi_jFwd.^2)); %TODO 14: Line to complete
    
    C1 = mu ./ (sqrt(eta^2 + phi_icent.^2 + phi_jFwd.^2));
    C2 = mu ./ (sqrt(eta^2 + phi_icent.^2 + phi_jBwd.^2));
    C3 = mu ./ (sqrt(eta^2 + phi_iFwd.^2 + phi_jcent.^2));
    C4 = mu ./ (sqrt(eta^2 + phi_iBwd.^2 + phi_jcent.^2));
    
    %%Equation 22, for inner points
    K = (P(2:end-1,3:end) .* C1 + P(2:end-1,1:end-2) .* C2 + ...
         P(3:end,2:end-1) .* C3 + P(1:end-2,2:end-1) .* C4);
    differenceFromAverage = (- lambda1*(I-c1).^2 + lambda2*(I-c2).^2);
    phi = (phi + (dt*delta_phi) .* (K - nu + differenceFromAverage)) ./ ...
           (1. + (dt*delta_phi) .* (C1 + C2 + C3 + C4)); %TODO 15: Line to complete
            
    %Reinitialization of phi
    if reIni>0 && mod(nIter, reIni)==0
        indGT = phi >= 0;
        indLT = phi < 0;
        
        phi=double(bwdist(indLT) - bwdist(indGT));
        
        %Normalization [-1 1]
        nor = min(abs(min(phi(:))), max(phi(:)));
        phi=phi/nor;
    end
  
    %Diference. This stopping criterium has the problem that phi can
    %change, but not the zero level set, that it really is what we are
    %looking for.
    dif = sqrt(mean((phi(:) - phi_old(:)).^2));

    fprintf('iter: %d, dif: %.6f\n',nIter,dif);

    if nIter>=3 && dif <= tol
        break;
    end
          
    if mod(nIter, 10)==0
        %Plot the level sets surface
        subplot(1,2,1) 
            %The level set function
            surfc(phi, 'LineStyle','none')  %TODO 16: Line to complete 
            hold on
            %The zero level set over the surface
            contour(phi, 'LineColor', 'blue'); %TODO 17: Line to complete
            hold off
            title('Phi Function');

        %Plot the curve evolution over the image
        subplot(1,2,2)
            imagesc(I);        
            colormap gray;
            hold on;
            contour(phi, 'LineColor', 'blue') %TODO 18: Line to complete
            title('Image and zero level set of Phi')

            axis off;
            hold off
        drawnow;
    end
end

seg = phi > 0;
