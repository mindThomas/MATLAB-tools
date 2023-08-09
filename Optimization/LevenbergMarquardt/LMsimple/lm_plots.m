function lm_plots ( t, y_dat, y_fit, sigma_y, cvg_hst, filename )
% lm_plots ( t, y_dat, y_fit, sigma_y, cvg_hst, filename )
% Plot statistics of the results of a Levenberg-Marquardt least squares
% analysis with lm.m

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 2 May 2016


dkgrn = [ 0.0 , 0.4 , 0.0 ];

y_dat = y_dat(:);
y_fit = y_fit(:);

[max_it,n] = size(cvg_hst); n = n-3;

figure(101);  % plot convergence history of parameters, reduced chi^2, lambda
 clf
subplot(211)
  plot( cvg_hst(:,1), cvg_hst(:,2:n+1), '-o','linewidth',4);
  for i=1:n
   text(1.02*cvg_hst(max_it,1),cvg_hst(max_it,1+i), sprintf('%d',i) );
  end
   ylabel('parameter values')
subplot(212)
  semilogy( cvg_hst(:,1) , [ cvg_hst(:,n+2) cvg_hst(:,n+3)], '-o','linewidth',4)
   text(cvg_hst(1,1),cvg_hst(1,n+2), '\chi^2_\nu','FontSize',16,'color','k');
   text(cvg_hst(1,1),cvg_hst(1,n+3), '\lambda', 'FontSize',16, 'color','k');
   text(cvg_hst(max_it,1),cvg_hst(max_it,n+2), '\chi^2_\nu','FontSize',16,'color','k');
   text(cvg_hst(max_it,1),cvg_hst(max_it,n+3), '\lambda', 'FontSize',16, 'color','k');
%  legend('\chi^2/(m-n+1)}', '\lambda', 3);
   ylabel('\chi^2_\nu and \lambda')
   xlabel('function calls')


figure(102); % ------------ plot data, fit, and confidence interval of fit
confidence_level = 0.99;   % confidence level for error confidence interval;
z = norminv( (1+confidence_level)/2);
 clf
   plot(t,y_dat,'og', t,y_fit,'-b', 
        t,y_fit+z*sigma_y,'.k', t,y_fit-z*sigma_y,'.k');
    legend('y_{data}','y_{fit}','99% c.i.','',0);
    ylabel('y(t)')
    xlabel('t')
% subplot(212)
%  semilogy(t,sigma_y,'-r','linewidth',4);
%    ylabel('\sigma_y(t)')

figure(103); % ------------ plot histogram of residuals, are they Gaussean?
 clf
 hist(real(y_dat - y_fit))
  title('histogram of residuals')
  axis('tight')
  xlabel('y_{data} - y_{fit}')
  ylabel('count')
