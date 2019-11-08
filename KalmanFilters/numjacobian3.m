%% NUMERICAL APPROXIMATION OF JACOBIAN
% 
% Numerical approximation of the Jacobian of a non-linear system
% in the form f(x(t),u(t),t) + G w(t), using five-point method
% or complex step arithmetic for improved accuracy. Note that this function 
% only works for function handles with two variables: x and u. Both can be vectors.
%
% Based on:
% Lai, Kok-Lam, et al. 
% "New complex-step derivative approximations with application to second-order kalman filtering." 
% AIAA Guidance, Navigation and Control Conference, San Francisco, California. 2005.
%
% M.A. van den Hoek
% Delft University of Technology
% Faculty of Aerospace Engineering
% Department of Control & Simulation
%

function [Jx,Ju,Jq] = numjacobian3(f,x,u,q)
    
    nx = size(x,1); % set size of state vector
    nu = size(u,1); % set size of input vector
    nq = size(q,1); % set size of input vector
    m = size(f(x,u,q),1); % number of functions in system
    step = 1e-8; % difference step size

    % preallocate matrix for speed
    Jx = zeros(m,nx);
    Ju = zeros(m,nu);
    Jq = zeros(m,nq);
    
%     for i=1:n
%        %xstep = x;
%        % compute the i-th component partial derivative 
%        % numerically using first order (Euler) forward difference approximation
%        %xstep(i) = x(i)+step;
%        %J(:,i) = (feval(f,xstep,0)-fx)/step;
%        
%        % Make step vector for i-th index by only setting variable xi
%        % non-zero to find df(x)/dxi
%        h = zeros(1,n); h(i) = step;
%        
%        % Five-point method, higher order first derivative approximation
%        xstep2 = x + 2*h;
%        xstep1 = x + h;
%        xstepm2 = x - 2*h;
%        xstepm1 = x - h;
%        J(:,i) = (-feval(f,xstep2,u) + 8*feval(f,xstep1,u) - 8*feval(f,xstepm1,u) + feval(f,xstepm2,u))/(12*step);
%     end
    
    
    % Alternatively, use complex step arithmetic for higher accuracy
    h = nx * eps;                            % differentiation step size, eps = machine epsilon
    for k = 1:nx                             % loop for each independent variable 
        x1 = x;                              % reference point
        x1(k) = x1(k) + h*1i;                % increment in kth independent variable
        Jx(:,k) = imag(f(x1,u,q))/h;            % complex step differentiation
    end
    
    h = nu * eps;                            % differentiation step size, eps = machine epsilon
    for k = 1:nu                             % loop for each independent variable 
        u1 = u;                              % reference point
        u1(k) = u1(k) + h*1i;                % increment in kth independent variable
        Ju(:,k) = imag(f(x,u1,q))/h;         % complex step differentiation
    end
    
    h = nq * eps;                            % differentiation step size, eps = machine epsilon
    for k = 1:nq                             % loop for each independent variable 
        q1 = q;                              % reference point
        q1(k) = q1(k) + h*1i;                % increment in kth independent variable
        Jq(:,k) = imag(f(x,u,q1))/h;         % complex step differentiation
    end
    
end
