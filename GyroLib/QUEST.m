function q = QUEST(fb, mb, fn, mn, wf, wm)
% q = QUEST(Sensors, Parameters)
% Function implements QUEST algorithm using measurements
% from three-component accelerometer with orthogonal axes and vector
% magnetometer
%
%   Input arguments:
%   fb  - Acceleration vector in body frame [3x1]
%   mb  - Magnetic field vector in body frame [3x1]
%   fn  - Gravity vector in navigation frame [3x1]
%   mn  - Magetic field vector in navigation frame [3x1]
%
%   Output arguments:
%   Cbn - estimated Direction Cosines Matrix (DCM)


%% Wahba's problem 
%%
% We need to find $C_n^b$ matrix to minimize loss function:

%%
% $L(A) = \frac{1}{2}\sum_i w_i |b_i-C_n^b n_i|$
%%
% $L(A) = \sum_i w_i - tr\left(AB^T\right)$
%%
% where:
%%
% $B=\sum_i w_i b_i n_i^T$
%% 
% Attitude matrix that maximizes $tr(AB^T)$ minimizes the loss function $L(A)$
%%
% $tr\left(AB^T\right) = q^TKq$
%%
% where:
%%
% $K = \left[\begin{array}{cc}B+B^T-tr(B)I & \sum_i w_i b_i \times n_i \\
% \sum_i w_i \left(b_i \times n_i\right)^T & tr(B)\end{array}\right]$
%%
% q-method finds the optimal quaternion as the normalized eigenvector of K
% with the largest eigenvalue:
%%
% $Kq_{opt}=\lambda_{max} q_{opt}$

%K matrix
B = wf*fb*fn' + wm*mb*mn';
K11 = B+B'-eye(3)*trace(B);
K22 = trace(B);
K12 = wf*cross(fb,fn)+wm*cross(mb,mn);
K21 = K12';
K = [K11, K12; K21, K22];

%Find eigenvalues and eigenvectors
[V,D] = eig(K);
lambdas = diag(D);

%Look for the largest eigenvalue
index = lambdas == max(lambdas);
q_opt = V(:,index);

%Transform q_opt to the Matlab's standard format
q = [q_opt(4), q_opt(1:3,1)'];
%Normalize quaternion
q = q/sqrt(q*q');

%As q-method yelds Cnb matrix
q = quat_conj(q);
end