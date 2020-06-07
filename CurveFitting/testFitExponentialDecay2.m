a_ = single(-12.2);
b_ = single(30);
c_ = single(-0.2);

std_dev = 0.3;

n = 10;
x = single(zeros(n,1));
y = single(zeros(n,1));
for (i = 1:n)
    x(i) = 2 + (i-1);
    y(i) = a_ + b_ * exp(c_ * x(i)) + 0*std_dev*randn(1,1);
end
	
figure(1);
plot(x, y, '*');
hold on;
fplot(@(x) a_ + b_ .* exp(c_.*x), [min(x), max(x)]);

[a2,b2,c2] = fitExponentialDecay2(x, y);
fplot(@(x) a2 + b2 .* exp(c2.*x), [min(x), max(x)]);
hold off;