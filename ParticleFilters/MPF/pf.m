% Standard particle filter 
% 
% A call to this function should look like
%
%       g = pf(y,m,N) 
%
% where
%          y:  Measurements
%          m:  Model structure
%          N:  Number of particles to use
%
%          g:  Returned structure with the following fields
%       g.xf:  Filtered states, i.e.  E[x(t) | y_1,..,y_{t}]
%     g.info:  Information about possible divergence (info = 0 - filter 
%              converged, info = 1 - filter diverged).
%
% Written by:
%              Thomas Sch�n (schon@isy.liu.se)
%              Division of Automatic Control
%              Link�ping University
%              www.control.isy.liu.se/~schon
%              Last revised on September 30, 2010
%

function g = pf(y,m,N)
  Tfinal = size(y,2);            % Number of samples (measurements)
  xf     = zeros(m.nx,Tfinal);   % Allocate room for the filtered estimates
  x      = repmat(m.x0,1,N) + sqrtm(m.P0)*randn(m.nx,N);
  for t=1:Tfinal
    yNow = y(:,t);
    e    = repmat(yNow,1,N) - feval(m.h,x);
    q    = exp(-(1/2)*(sum(e.*(inv(m.R)*e))));  % Compute the importance weights  
    if sum(q)>1e-12              % Divergence check
      q       = q/sum(q);        % Normalize the importance weights
      xf(:,t) = sum(repmat(q,m.nx,1).*x,2);
      index   = sysresample(q);  % Resample
      x       = x(:,index);
      info    = 0;
     else   % The filter has diverged
       info = 1;
       xf   = zeros(4,Tfinal);
       disp(['Weights close to zero at t=',num2str(t),' !!!']);
       return;
     end;
    x = feval(m.f,x) + m.T*sqrtm(m.Q)*randn(size(m.Q,1),N);  % Predict the particles one step forward
  end;
  g.xf   = xf;
end