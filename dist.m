clc
close all
clear all

beta = 1;
alph = 0:0.5:4;
x = 0:.1:3;

mu = 0;  % mean zero
sig = 0.5:0.5:1.5;  % standard deviation 

for i = 1:length(x)
    for j = 1:length(alph)
        K = -(x(i)/beta);
        G(i,j) = (beta^-alph(j))*x(i)^(alph(j)-1)*exp(K)/(gamma(alph(j)));  % gamma distribution 
        KG = -(x(i)/beta)^(alph(j));
        W(i,j) = alph(j)*(beta^-alph(j))*(x(i)^(alph(j)-1))*exp(KG);        % Weibull density 
        Wd(i,j) = 1-exp(KG);                                                % Weibull distribution 
       
    end
end

for i = 1:length(x)
    for j = 1: length(alph)-1
        T = -(x(i)/beta);
        TF(i,j) = 1-T*((x(i)/beta)^(j)/factorial(j));
    end
    TFF(i) = sum(TF(i,:));
end 
figure(1)
subplot(221)
plot(x,G)
legend('\alpha = 0','\alpha = 0.5','\alpha = 1','\alpha = 2','\alpha = 3'); title('Gamma density function '); ylabel('f(x)');
subplot(222)
plot(x,W)
legend('\alpha = 0','\alpha = 0.5','\alpha = 1','\alpha = 2','\alpha = 3'); title('Weibull density function'); ylabel('f(x)');
subplot(223)
plot(x,Wd)
legend('\alpha = 0','\alpha = 0.5','\alpha = 1','\alpha = 2','\alpha = 3'); title('Weibull distribution function'); ylabel('F(x)');
subplot(224)
plot(x,TF)
legend('\alpha = 0','\alpha = 0.5','\alpha = 1','\alpha = 2','\alpha = 3'); title('Gamma distribution function');ylabel('F(x)');


for i = 1:length(x)
    for j = 1:length(sig)
        KK = -((log(x(i))- mu)^2)/(2*sig(j)^2);
        F(i,j) = (1/(x(i)*sqrt(2*pi*sig(j)^2)))*(exp(KK));  % log normal density function 
    end 
end 

figure(2)
plot(x,F)
legend('\sigma = 0.5','\sigma = 1','\sigma = 1.5');
