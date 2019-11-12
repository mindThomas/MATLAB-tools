mu = [0
      1];
C = rot2(deg2rad(80)) * diag([1,2]) * rot2(deg2rad(80))';  

% Create bivariate joint distribution between x and y
% p(x,y)
joint = Gaussian(mu, C);

% Marginalize out variable x to get
% p(y)
py = joint.marginalize(1);

% Compute conditional distribution at a few places
% p(x|y) = p(x,y) / p(y)
y = 1;

conditional = joint.conditional(2, y);
conditional.plot(figure(1));

comparison = joint.pdf([-1, y]') / py.pdf(y)

