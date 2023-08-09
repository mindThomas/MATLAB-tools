function y_hat = lm_func(t,p,c)
% y_hat = lm_func(t,p,c)
%  
% three example functions used for nonlinear least squares curve-fitting
% in the tutorial file lm_examp.m
% to demonstrate the Levenberg-Marquardt function, lm.m
%
% -------- INPUT VARIABLES ---------
%  t     = independent variable values (assumed to be error-free)        (m x 1)
%  p     = parameter values , n = 4 in these examples                    (n x 1)
%  c     = optional vector of other constants
% 
% ---------- OUTPUT VARIABLES -------
% y_hat  = curve-fit fctn evaluated at points t and with parameters p    (m x 1)
%
% H.P. Gavin, Dept. Civil & Environ. Eng'g, Duke Univ.

global example_number

%  example 1:  easy for LM ... a poor initial guess is ok
if example_number == 1
   y_hat = p(1)*exp(-t/p(2)) + p(3)*t.*exp(-t/p(4));
end

%  example 2: medium for LM ... local minima
if example_number == 2
   mt = max(t);
   y_hat = p(1)*(t/mt) + p(2)*(t/mt).^2 + p(3)*(t/mt).^3 + p(4)*(t/mt).^4;
end

% example 3: difficult for LM ... needs a very good initial guess
if example_number == 3
   y_hat = p(1)*exp(-t/p(2)) + p(3)*sin(t/p(4));
end

% LM_FUNC ------------------------------------ 13 Apr 2011, 27 Jun 2011
