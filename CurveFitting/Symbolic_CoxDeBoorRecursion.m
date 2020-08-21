% Cox-de Boor recursion formula
% Used to evaluate the basis function at a given u associated with
% a particular control point, i, for a B-spline basis of order p

% Knots
u = sym('u');
ui = sym('u%d', [7, 1])
s = [];
for (i = 1:(length(ui)-1))
    s = [s; sym(sprintf('s%d%d', i, i+1))];
end
s

%%
% Cox-de Boor recursion formula
% Build the recursion tree.
layer = [];
% Starting with the first layer which is just the interval functions
% defined in 's'
p = 0;
layer{p +1} = s;

% Compute the next n layers
n = 3;

for (l = 1:n)
    % Next layer is then given as:
    p = p + 1;
    layer{p +1} = [];
    for (i = 1:(length(layer{p-1 +1})-1))
        layer{p +1} = [layer{p +1};
              (u - ui(i)) / (ui(i+p) - ui(i)) * layer{p-1 +1}(i) + ...
              (ui(i+p+1) - u) / (ui(i+p+1) - ui(i+1)) * layer{p-1 +1}(i+1) ...
            ];
    end
end

%%
% For a given order basis function, p, derive the activation from control point, j, 
% in the region u_i <= u < u_i+1
p = 3; % 1st order basis function
j = 0; % activation from the first control point
i = 0; % in the span u0 to u1

activation = simplify(diff(layer{p +1}(j +1), s(i +1)));

% Offset by p-1 to get a center-aligned activation function
%ss = sym('s');
%expand(subs(activation, u, ss+p-1))

%% Compute the derivate with respect to u
dlayer = {};
dlayer{p +1} = simplify(diff(layer{p +1}, u));
dactivation = simplify(diff(dlayer{p +1}(j +1), s(i +1)))