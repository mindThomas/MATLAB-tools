function [p,redX2,sigma_p,sigma_y,corr_p,R_sq,cvg_hst] = lm(func,p,t,y_dat,weight,dp,p_min,p_max,c,opts)
% [p,redX2,sigma_p,sigma_y,corr_p,R_sq,cvg_hst] = lm(func,p,t,y_dat,weight,dp,p_min,p_max,c,opts)
%
% Levenberg Marquardt curve-fitting: minimize sum of weighted squared residuals
% ----------  INPUT  VARIABLES  -----------
% func   = function of n independent variables, 't', and m parameters, 'p', 
%          returning the simulated model: y_hat = func(t,p,c)
% p      = initial guess of parameter values                             (n x 1)
% t      = independent variables (used as arg to func)                   (m x 1)
% y_dat  = data to be fit by func(t,p)                                   (m x 1)
% weight = weights or a scalar weight value ( weight >= 0 ) ...          (m x 1)
%          inverse of the standard measurement errors
%          Default:  ( 1 / ( y_dat' * y_dat ))
% dp     = fractional increment of 'p' for numerical derivatives
%          dp(j)>0 central differences calculated
%          dp(j)<0 one sided 'backwards' differences calculated
%          dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
%          Default:  0.001;
% p_min  = lower bounds for parameter values                             (n x 1)
% p_max  = upper bounds for parameter values                             (n x 1)
% c      = an optional matrix of values passed to func(t,p,c)
% opts   = vector of algorithmic parameters
%             parameter    defaults    meaning
% opts(1)  =  prnt            3        >1 intermediate results; >2 plots
% opts(2)  =  MaxIter      10*Npar     maximum number of iterations
% opts(3)  =  epsilon_1       1e-3     convergence tolerance for gradient
% opts(4)  =  epsilon_2       1e-3     convergence tolerance for parameters
% opts(5)  =  epsilon_3       1e-1     convergence tolerance for red. Chi-square
% opts(6)  =  epsilon_4       1e-1     determines acceptance of a L-M step
% opts(7)  =  lambda_0        1e-2     initial value of L-M paramter
% opts(8)  =  lambda_UP_fac   11       factor for increasing lambda
% opts(9)  =  lambda_DN_fac    9       factor for decreasing lambda
% opts(10) =  Update_Type      1       1: Levenberg-Marquardt lambda update
%                                      2: Quadratic update 
%                                      3: Nielsen's lambda update equations
%
% ----------  OUTPUT  VARIABLES  -----------
% p       = least-squares optimal estimate of the parameter values
% redX2   = reduced Chi squared error criteria - should be close to 1
% sigma_p = asymptotic standard error of the parameters
% sigma_y = asymptotic standard error of the curve-fit
% corr_p  = correlation matrix of the parameters
% R_sq    = R-squared cofficient of multiple determination  
% cvg_hst = convergence history ... see lm_plots.m
 
%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 4 May 2016
%   modified from: http://octave.sourceforge.net/optim/function/leasqr.html
%   using references by
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.
%   Sam Roweis       http://www.cs.toronto.edu/~roweis/notes/lm.pdf
%   Manolis Lourakis http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf
%   Hans Nielson     http://www2.imm.dtu.dk/~hbn/publ/TR9905.ps
%   Mathworks        optimization toolbox reference manual
%   K. Madsen, H.B., Nielsen, and O. Tingleff
%   http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf

 global   iteration  func_calls

 tensor_parameter = 0;                  % set to 1 of parameter is a tensor

 iteration  = 0;                        % iteration counter
 func_calls = 0;                        % running count of function evaluations

 p = p(:); y_dat = y_dat(:); t = t(:);  % make column vectors
 Npar   = length(p);                    % number of parameters
 Npnt   = length(y_dat);                % number of data points
 p_old  = zeros(Npar,1);                % previous set of parameters
 y_old  = zeros(Npnt,1);                % previous model, y_old = y_hat(t;p_old)
 X2     = 1e-3/eps;                     % a really big initial Chi-sq value
 X2_old = 1e-3/eps;                     % a really big initial Chi-sq value
 J      = zeros(Npnt,Npar);             % Jacobian matrix
 DoF    = Npnt - Npar + 1;              % statistical degrees of freedom


 if length(t) ~= length(y_dat)
    disp('lm.m error: the length of t must equal the length of y_dat');
    length_t = length(t)
    length_y_dat = length(y_dat)
    X2 = 0; corr_p = 0; sigma_p = 0; sigma_y = 0; R_sq = 0; cvg_hist = 0;
    if ~tensor_parameter, 
         return;                
    end
 end

 if nargin <  5, weight = 1/(y_dat'*y_dat); end
 if nargin <  6, dp = 0.001; end
 if nargin <  7, p_min   = -100*abs(p); end
 if nargin <  8, p_max   =  100*abs(p); end
 if nargin <  9, c       =  1; end
 if nargin < 10,                % Algorithmic Paramters
%         prnt MaxIter  eps1  eps2  eps3  eps4  lam0  lamUP lamDN UpdateType 
   opts = [  3,10*Npar, 1e-3, 1e-3, 1e-1, 1e-1, 1e-2,    11,    9,        1 ];
 end
 prnt          = opts(1);        % >1 intermediate results; >2 plots
 MaxIter       = opts(2);        % maximum number of iterations
 epsilon_1     = opts(3);        % convergence tolerance for gradient
 epsilon_2     = opts(4);        % convergence tolerance for parameters
 epsilon_3     = opts(5);        % convergence tolerance for Chi-square
 epsilon_4     = opts(6);        % determines acceptance of a L-M step
 lambda_0      = opts(7);        % initial value of damping paramter, lambda
 lambda_UP_fac = opts(8);        % factor for increasing lambda
 lambda_DN_fac = opts(9);        % factor for decreasing lambda
 Update_Type   = opts(10);       % 1: Levenberg-Marquardt lambda update
                                 % 2: Quadratic update 
                                 % 3: Nielsen's lambda update equations

%if ( tensor_parameter && prnt == 3 ) prnt = 2; end

 plotcmd='figure(102); plot(t(:,1),y_init,''-k'',t(:,1),y_hat,''-b'',t(:,1),y_dat,''o'',''color'',[0,0.6,0],''MarkerSize'',4); title(sprintf(''\\chi^2_\\nu = %f'',X2/DoF)); drawnow ';

 p_min=p_min(:); p_max=p_max(:); % make column vectors

 if length(dp) == 1, dp = dp*ones(Npar,1); end

 idx   = find(dp ~= 0);                  % indices of the parameters to be fit
 Nfit = length(idx);                     % number of parameters to fit
 stop = 0;                               % termination flag

 y_init = feval(func,t,p,c);       % residual error using p_try

 if ( var(weight) == 0 )                 % identical weights vector 
        weight = abs(weight(1))*ones(Npnt,1);        
        disp('using uniform weights for error analysis')
 else
         weight = abs(weight(:));
 end

% initialize Jacobian with finite difference calculation
 [JtWJ,JtWdy,X2,y_hat,J] = lm_matx(func,t,p_old,y_old,1,J,p,y_dat,weight,dp,c);

 if ( max(abs(JtWdy)) < epsilon_1 )
        fprintf(' *** Your Initial Guess is Extremely Close to Optimal ***\n')
        fprintf(' *** epsilon_1 = %e\n', epsilon_1);
        stop = 1;
 end

 switch Update_Type
  case 1                                 % Marquardt: init'l lambda
        lambda  = lambda_0;
  otherwise                              % Quadratic and Nielsen
        lambda  = lambda_0 * max(diag(JtWJ)); nu=2;
 end

 X2_old = X2;                            % previous value of X2 

 cvg_hst = ones(MaxIter,Npar+3);         % initialize convergence history

 while ( ~stop && iteration <= MaxIter )        % --- Start Main Loop

   iteration = iteration + 1;
 
% incremental change in parameters
   switch Update_Type
     case 1                                       % Marquardt
       h = ( JtWJ + lambda*diag(diag(JtWJ)) ) \ JtWdy;
     otherwise                                    % Quadratic and Nielsen
       h = ( JtWJ + lambda*eye(Npar) ) \ JtWdy;
   end

%  big = max(abs(h./p)) > 2;                      % this is a big step

   % --- Are parameters [p+h] much better than [p] ?

   p_try = p + h(idx);                            % update the [idx] elements 
   p_try = min(max(p_min,p_try),p_max);           % apply constraints

   delta_y = y_dat - feval(func,t,p_try,c);       % residual error using p_try
   if ~all(isfinite(delta_y))                     % floating point error; break
     stop = 1;
     break     
   end
   func_calls = func_calls + 1;
   X2_try = delta_y' * ( delta_y .* weight );     % Chi-squared error criteria

   if ( Update_Type == 2 )                        % Quadratic
%    One step of quadratic line update in the h direction for minimum X2
     alpha =  JtWdy'*h / ( (X2_try - X2)/2 + 2*JtWdy'*h ) ;
     h = alpha * h;

     p_try = p + h(idx);                          % update only [idx] elements
     p_try = min(max(p_min,p_try),p_max);         % apply constraints

     delta_y = y_dat - feval(func,t,p_try,c);     % residual error using p_try
     func_calls = func_calls + 1;
     X2_try = delta_y' * ( delta_y .* weight );   % Chi-squared error criteria
   end

%  switch Update_Type                             % Nielsen
%    case 1
%      rho = (X2 - X2_try) / ( h' * (lambda*diag(diag(JtWJ))*h + JtWdy) );
%    otherwise
       rho = (X2 - X2_try) / ( h' * (lambda * h + JtWdy) );
%  end

   if ( rho > epsilon_4 )                         % it IS significantly better

     dX2 = X2 - X2_old;
     X2_old = X2;
     p_old = p;
     y_old = y_hat;
     p = p_try(:);                           % accept p_try

     [JtWJ,JtWdy,X2,y_hat,J] = ...
                       lm_matx(func,t,p_old,y_old,dX2,J,p,y_dat,weight,dp,c);

                                % decrease lambda ==> Gauss-Newton method

     switch Update_Type
       case 1                                   % Levenberg
         lambda = max(lambda/lambda_DN_fac,1.e-7);
       case 2                                   % Quadratic
         lambda = max( lambda/(1 + alpha) , 1.e-7 );
       case 3                                   % Nielsen
         lambda = lambda*max( 1/3, 1-(2*rho-1)^3 ); nu = 2;
     end

     if ( prnt > 2 )
       eval(plotcmd);
     end

   else                                           % it IS NOT better

     X2 = X2_old;                             % do not accept p_try

     if ( ~rem(iteration,2*Npar) )            % rank-1 update of Jacobian
       [JtWJ,JtWdy,dX2,y_hat,J] = ...
                       lm_matx(func,t,p_old,y_old,-1,J,p,y_dat,weight,dp,c);
     end

                                % increase lambda  ==> gradient descent method

     switch Update_Type
       case 1                                   % Levenberg
         lambda = min(lambda*lambda_UP_fac,1.e7);
       case 2                                   % Quadratic
         lambda = lambda + abs((X2_try - X2)/2/alpha);
       case 3                                   % Nielsen
         lambda = lambda * nu;   nu = 2*nu;
     end

   end

   if ( prnt > 1 )
     fprintf('>%3d:%3d | chi_sq=%10.3e | lambda=%8.1e \n', ...
                                       iteration,func_calls,X2/DoF,lambda );
     fprintf('    param:  ');
     for pn=1:Npar
       fprintf(' %10.3e', p(pn) );
     end
     fprintf('\n');
     fprintf('    dp/p :  ');
     for pn=1:Npar
       fprintf(' %10.3e', h(pn) / p(pn) );
     end
     fprintf('\n');
   end

% update convergence history ... save _reduced_ Chi-square
   cvg_hst(iteration,:) = [ func_calls  p'  X2/DoF lambda ];


   if ( max(abs(JtWdy)) < epsilon_1  &&  iteration > 2 ) 
     fprintf(' **** Convergence in r.h.s. ("JtWdy")  **** \n')
     fprintf(' **** epsilon_1 = %e\n', epsilon_1);
     stop = 1;
   end
   if ( max(abs(h)./(abs(p)+1e-12)) < epsilon_2  &&  iteration > 2 ) 
     fprintf(' **** Convergence in Parameters **** \n')
     fprintf(' **** epsilon_2 = %e\n', epsilon_2);
     stop = 1;
   end
   if ( X2/DoF < epsilon_3 &&  iteration > 2 ) 
     fprintf(' **** Convergence in reduced Chi-square  **** \n')
     fprintf(' **** epsilon_3 = %e\n', epsilon_3);
     stop = 1;
   end
   if ( iteration == MaxIter )
     disp(' !! Maximum Number of Iterations Reached Without Convergence !!')
     stop = 1;
   end

 end                                        % --- End of Main Loop

 % --- convergence achieved, find covariance and confidence intervals

% ---- Error Analysis ----

 if var(weight) == 0   % recompute equal weights for paramter error analysis 
   weight = DoF/(delta_y'*delta_y) * ones(Npnt,1);
 end

 if nargout > 1                                % reduced Chi-square
   redX2 = X2 / DoF;
 end

 [JtWJ,JtWdy,X2,y_hat,J] = lm_matx(func,t,p_old,y_old,-1,J,p,y_dat,weight,dp,c);

 if nargout > 2                                % standard error of parameters 
   covar_p = inv(JtWJ);
   sigma_p = sqrt(diag(covar_p));
 end

 if nargout > 3                                % standard error of the fit
%  sigma_y = sqrt(diag(J * covar_p * J'));        % slower version of below
   sigma_y = zeros(Npnt,1);
   for i=1:Npnt
     sigma_y(i) = J(i,:) * covar_p * J(i,:)';        
   end
   sigma_y = sqrt(sigma_y);
 end

 if nargout > 4                                % parameter correlation matrix
   corr_p = covar_p ./ [sigma_p*sigma_p'];        
 end

 if nargout > 5                                % coefficient of multiple determination
   R_sq = corr([y_dat y_hat]);
   R_sq = R_sq(1,2).^2;                
 end

 if nargout > 6                                % convergence history
   cvg_hst = cvg_hst(1:iteration,:);
 end

% endfunction  # ---------------------------------------------------------- LM


function J = lm_FD_J(func,t,p,y,dp,c)
% J = lm_FD_J(func,t,p,y,{dp},{c})
%
% partial derivatives (Jacobian) dy/dp for use with lm.m
% computed via Finite Differences
% Requires n or 2n function evaluations, n = number of nonzero values of dp
% -------- INPUT VARIABLES ---------
% func = function of independent variables, 't', and parameters, 'p',
%        returning the simulated model: y_hat = func(t,p,c)
% t  = independent variables (used as arg to func)                       (m x 1)
% p  = current parameter values                                          (n x 1)
% y  = func(t,p,c) initialised by user before each call to lm_FD_J       (m x 1)
% dp = fractional increment of p for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
%      Default:  0.001;
% c  = optional vector of constants passed to y_hat = func(t,p,c)
%---------- OUTPUT VARIABLES -------
% J  = Jacobian Matrix J(i,j)=dy(i)/dp(j)         i=1:n; j=1:m 

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.


 global  func_calls

 m=length(y);              % number of data points
 n=length(p);              % number of parameters

 if nargin < 5
        dp = 0.001*ones(1,n);
 end

 ps=p; J=zeros(m,n); del=zeros(n,1);         % initialize Jacobian to Zero

 for j=1:n                 % START --- loop over all parameters

   del(j) = dp(j) * (1+abs(p(j)));   % parameter perturbation
   p(j)   = ps(j) + del(j);          % perturb parameter p(j)

   if del(j) ~= 0
     y1=feval(func,t,p,c);
     func_calls = func_calls + 1;

     if (dp(j) < 0)                  % backwards difference
       J(:,j) = (y1-y)./del(j);
     else                            % central difference, additional func call
       p(j) = ps(j) - del(j);
       J(:,j) = (y1-feval(func,t,p,c)) ./ (2 .* del(j));
       func_calls = func_calls + 1;
     end
   end

   p(j)=ps(j);                       % restore p(j)

 end                       % END --- loop over all parameters

% endfunction # -------------------------------------------------- LM_FD_J



function J = lm_Broyden_J(p_old,y_old,J,p,y)
% J = lm_Broyden_J(p_old,y_old,J,p,y)
% carry out a rank-1 update to the Jacobian matrix using Broyden's equation
%---------- INPUT VARIABLES -------
% p_old = previous set of parameters                                     (n x 1)
% y_old = model evaluation at previous set of parameters, y_hat(t;p_old) (m x 1)
% J  = current version of the Jacobian matrix                            (m x n)
% p     = current  set of parameters                                     (n x 1)
% y     = model evaluation at current  set of parameters, y_hat(t;p)     (m x 1)
%---------- OUTPUT VARIABLES -------
% J = rank-1 update to Jacobian Matrix J(i,j)=dy(i)/dp(j)  i=1:n; j=1:m  (m x n)

 h  = p - p_old;

 J = J + ( y - y_old - J*h )*h' / (h'*h);       % Broyden rank-1 update eq'n

% endfunction # ---------------------------------------------- LM_Broyden_J



function [JtWJ,JtWdy,Chi_sq,y_hat,J] = lm_matx(func,t,p_old,y_old,dX2,J,p,y_dat,weight,dp,c)
% [JtWJ,JtWdy,Chi_sq,y_hat,J] = lm_matx(func,t,p_old,y_old,dX2,J,p,y_dat,weight,{da},{c})
%
% Evaluate the linearized fitting matrix, JtWJ, and vector JtWdy, 
% and calculate the Chi-squared error function, Chi_sq 
% Used by Levenberg-Marquard algorithm, lm.m   
% -------- INPUT VARIABLES ---------
% func   = function ofpn independent variables, p, and m parameters, p,
%         returning the simulated model: y_hat = func(t,p,c)
% t      = independent variables (used as arg to func)                   (m x 1)
% p_old  = previous parameter values                                     (n x 1)
% y_old  = previous model ... y_old = y_hat(t;p_old);                    (m x 1)
% dX2    = previous change in Chi-squared criteria                       (1 x 1)
% J      = Jacobian of model, y_hat, with respect to parameters, p       (m x n)
% p      = current  parameter values                                     (n x 1)
% y_dat  = data to be fit by func(t,p,c)                                 (m x 1)
% weight = the weighting vector for least squares fit ...
%          inverse of the squared standard measurement errors
% dp     = fractional increment of 'p' for numerical derivatives
%          dp(j)>0 central differences calculated
%          dp(j)<0 one sided differences calculated
%          dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
%          Default:  0.001;
% c      = optional vector of constants passed to y_hat = func(t,p,c)
%---------- OUTPUT VARIABLES -------
% JtWJ    = linearized Hessian matrix (inverse of covariance matrix)     (n x n)
% JtWdy   = linearized fitting vector                                    (n x m)
% Chi_sq = Chi-squared criteria: weighted sum of the squared residuals WSSR
% y_hat  = model evaluated with parameters 'p'                           (m x 1)
% J      = Jacobian of model, y_hat, with respect to parameters, p       (m x n)

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

 global   iteration  func_calls

 Npnt = length(y_dat);               % number of data points
 Npar = length(p);                   % number of parameters 

 if nargin < 6
   dp = 0.001;
 end

%JtWJ = zeros(Npar);
%JtWdy  = zeros(Npar,1);

 y_hat = feval(func,t,p,c);          % evaluate model using parameters 'p'
 func_calls = func_calls + 1;

 if ( ~rem(iteration,2*Npar) || dX2 > 0 ) 
   J = lm_FD_J(func,t,p,y_hat,dp,c);            % finite difference
 else
   J = lm_Broyden_J(p_old,y_old,J,p,y_hat);     % rank-1 update
 end

 delta_y = y_dat - y_hat;            % residual error between model and data

 Chi_sq = delta_y' * ( delta_y .* weight );     % Chi-squared error criteria

 JtWJ  = J' * ( J .* ( weight * ones(1,Npar) ) );  

 JtWdy = J' * ( weight .* delta_y );
 
% endfunction  # ------------------------------------------------------ LM_MATX
