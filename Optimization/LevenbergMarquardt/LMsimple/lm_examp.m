% Levenberg-Marquardt example/test

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 2 May 2016
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

clear all
hold off
randn('seed',0);	% specify a particular random sequence for msmnt error

% *** For this demonstration example, simulate some artificial measurements by
% *** adding random errors to the curve-fit equation.  

global	example_number

example_number = 3;			  % which example to run: 1 2 or 3  

consts = [ ];                             % optional vector of constants

Npnt = 100;				  % number of data points

t = [1:Npnt]';				  % independent variable, column vector
% true value of parameters ...              column vectors
if example_number == 1, p_true  = [ 20   10   1  50 ]'; end	
if example_number == 2, p_true  = [ 20  -24  30 -40 ]'; end	 
if example_number == 3, p_true  = [  6   20   1   5 ]'; end	

Npar = length(p_true);

msmnt_err = 0.5;
y_dat = lm_func(t,p_true,consts);
y_dat = y_dat + msmnt_err*randn(Npnt,1);  % add random measurement errors


% range of values for gridded paramter search
p1 = 0.1*p_true(1):0.1*p_true(1):2*p_true(1);
p2 = 0.1*p_true(2):0.1*p_true(2):2*p_true(2);
p3 = 0.1*p_true(3):0.1*p_true(3):2*p_true(3);
p4 = 0.1*p_true(4):0.1*p_true(4):2*p_true(4);

% gridded parameter search
for ip2 = 1:length(p2);
   for ip4 = 1:length(p4);
	pt = [ p_true(1)  p2(ip2) p_true(3) p4(ip4) ];
	delta_y = ( y_dat - lm_func(t,pt,consts) );
	X2(ip2,ip4) = (delta_y' * delta_y)/(Npnt-Npar);  % reduced chi-sq.
   end
end

figure(1); % ------------ plot shape of Chi-squared objective function
 clf
 mesh(p2,p4,log10(X2))
  xlabel('p_2')
  ylabel('p_4')
  zlabel('log_{10}(\chi^2_\nu)','rotation',90)
  hold on
plotfile = ['lm_exampA',int2str(example_number),'.eps'];

% *** Replace the lines above with a read-in of some
% *** experimentally-measured data.

% initial guess parameters  ...
if example_number == 1, p_init  = [  5   2  0.2  10 ]';  end
if example_number == 2, p_init  = [  4  -5  6    10 ]';  end
if example_number == 3, p_init  = [ 10  50  5   5.7 ]';  end

% the proper value of the weight is 1/(squared measurement error)
 weight = 1/msmnt_err^2;

p_min = -10*abs(p_init);
p_max =  10*abs(p_init);

%{
figure(3)
 clf
 plot(t,y_dat,'o');
  xlabel('t')
  ylabel('y(t)')
%}

% Algorithmic Parameters
%         prnt MaxIter  eps1  eps2  eps3  eps4  lam0  lamUP lamDN UpdateType 
   opts = [  3,    100, 1e-3, 1e-3, 1e-1, 1e-1, 1e-2,    11,    9,        1 ];

[p_fit,Chi_sq,sigma_p,sigma_y,corr,R_sq,cvg_hst] = ...
	lm('lm_func',p_init,t,y_dat,weight,-0.01,p_min,p_max,consts,opts);

y_fit = lm_func(t,p_fit,consts);

disp('    initial    true       fit        sigma_p percent')
disp(' -------------------------------------------------------')
disp ([ p_init  p_true  p_fit sigma_p 100*abs(sigma_p./p_fit) ])

n = length(p_fit);

Chi_sq
R_sq
corr

fn = sprintf('lm_examp_%d',example_number);

lm_plots ( t, y_dat, y_fit, sigma_y, cvg_hst, fn );

% lm_examp.m  ---------------------------------------------------------------- 
