x = Gaussian(10, 2);
noise = Gaussian(0, 1);
deterministic_offset = Gaussian(1, 0); % 0 variance = deterministic
% y = A*x + b
A = 99/10;
b = deterministic_offset.add(noise);

% z = [x, y]
z = x.join_transform(A, b)
y = z.marginalize(1)

% p(x|y)
z.conditional(2, 0)
q = z.conditional(2)

% p(x,y) = p(x|y) * p(y)
joint = q.join(y)