% Rao-Blackwellized (Marginalized) particle filter 
% 
% A call to this function should look like
%
%       g = mpf(y,m,N) 
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
%              Thomas Schön (schon@isy.liu.se)
%              Division of Automatic Control
%              Linköping University
%              www.control.isy.liu.se/~schon
%              Last revised on August 19, 2011
%

function g = mpf(y,m,N)
  Tfinal = size(y,2);         % Number of samples (measurements)
  NrPart = N;                 % Number of particles
  xf     = zeros(4,Tfinal);   % Allocate room for the filtered estimates

  % (1) Initialize
  xnp = repmat(m.x0(1),1,NrPart) + chol(m.P0(1,1))*randn(1,NrPart);  % Nonlinear states
  xlp = repmat(m.x0(2:4),1,NrPart);          % Conditionally linear Gaussian states
  Pl  = m.P0(2:4,2:4);
  Pp  = repmat(Pl,[1,1,NrPart]);             % Initial covariance matrix
  xlf = zeros(size(xlp));                    % Allocate room for the filtered quantities
  Pf  = zeros(size(Pp));                     % Allocate room for the filtered quantities
  for t=1:Tfinal
    yNow = y(:,t);
    yhat = [(0.1*xnp.^2).*sign(xnp); xlp(1,:) - xlp(2,:) + xlp(3,:)];
    e    = repmat(yNow,1,NrPart) - yhat;
    % (2) Compute the importance weights according to eq. (25a)
    for i=1:NrPart
      M = m.C*Pp(:,:,i)*m.C' + m.R;
      q(i) = exp(-(1/2)*(e(:,i)'*inv(M)*e(:,i)));
    end;
    if(sum(q)>1e-12)                % Divergence check
      q       = q/sum(q);           % Normalize the importance weights
      xf(1,t) = sum(q.*xnp,2);      % Compute estimate for the nonlinear states
      index   = sysresample(q);     % (3) Resample
      xnp     = xnp(:,index);       % Resampled nonlinear particles
      xlp     = xlp(:,index);       % Resampled linear particles
      Pp      = Pp(:,:,index);      % Resampled covariance matrices
      xlf     = xlf(:,index);       % Resampled linear particles
      Pf      = Pf(:,:,index);      % Resampled covariance matrices     
      info    = 0;
    else   % The filter has diverged
      info = 1;
      xf   = zeros(4,Tfinal);
      disp(['Weights close to zero at t=',num2str(t),' !!!']);
      return;
    end;
    % (4a) KF MU
    for i = 1:NrPart
      M         = m.C*Pp(:,:,i)*m.C' + m.R;    % Eq. (22c)
      K         = Pp(:,:,i)*m.C'*inv(M);       % Eq. (22d)
      yhat      = [(0.1*xnp(:,i).^2).*sign(xnp(:,i)); xlp(1,i) - xlp(2,i) + xlp(3,i)];
      xlf(:,i)  = xlp(:,i) + K*(yNow - yhat);  % Eq. (22a)
      Pf(:,:,i) = Pp(:,:,i) - K*M*K';          % Eq. (22b)
    end;
    xf(2:4,t) = mean(xlf,2);    % Compute estimate for the linear states
    % (4b) PF prediction according to Eq. (25b)
    xnf = xnp;
    for i = 1:NrPart
      xnp(i) = atan(xnf(i)) + xlf(1,i) + sqrt(m.An*Pf(:,:,i)*m.An' + m.Q(1,1))*randn(1);
    end;
    % (4c) KF TU
    for i = 1:NrPart
      N         = m.An*Pf(:,:,i)*m.An' + m.Q(1,1);       % Eq. (23c)
      L         = m.Al*Pf(:,:,i)*m.An'*inv(N);           % Eq. (23d)
      z         = xnp(i) - atan(xnf(i));                 % Eq. (24a)
      xlp(:,i)  = m.Al*xlf(:,i) + L*(z - m.An*xlf(:,i)); % Eq. (23a)
      Pp(:,:,i) = m.Al*Pf(:,:,i)*m.Al' + m.Q(2:end,2:end) - L*N*L'; % Eq. (23b)
    end;
  end;
  g.xf   = xf;
  g.info = info;
end
