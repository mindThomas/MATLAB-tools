function SOC = segop(k,Q,x0)

% function SOC = segop(k,Q,xO). This function returns the
% parameters for the segmentally optimal composite curve.
% It receives the .IG curve parameters, x0, data points Q,
% and knot sequence k. It was written by E. J. Lane.

[P,ang,dt] = ktangdt(x0);   % Separates the vector xO
                            % into its subcomponents.

bdt = bstdst(dt,Q,P,ang,k); % Call to the function which finds
                            % the optimum distances for a segment.

bdt1 = [];                  % ECR

for i = 1 : 2 % Loop to assemble the "best" distances.
    bdt1 =  [bdt1 bdt(i,:)];
end

% Assemble the vector of parameters for the curve.
SOC = [P(1,:) P(2,:) ang bdt1];
