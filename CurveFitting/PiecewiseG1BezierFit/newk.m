function nk = newk(Q,P)

% function nk = newk(Q,P). This function takes data points, Q,and
% knot points, P, as input. The function finds the closest data
% point of the cubic segment that is associated with that knot
% point, and returns a new k-array, nk. It ensures the data points
% of Q are properly associated with the proper segment of the curve.
% "dpkpc" is a global variable that is initially equal to the old
% k. This was written by M. R. Holmes and revised by E. J. Lane.

global dpkpc

[r,m] = size(Q);
[s,n] = size(P);

nk(1) = 1; nk(n) = m; 	% Ensures the knot sequence starts and
			% ends with the 1st and last points in Q.
for i = 2:n-1
	js = dpkpc(i-1); 
	je = dpkpc(i+1); % variables to pick out
	jm = dpkpc(i); 	 % interior knot positions.
	z  = je - js + 1; 
	mm = jm - js + 1;
	R = Q(:,js:je) - P(:,i) * ones(1,z); % Finds differences
					     % between data points and
					     % knot point being checked.
	for jj = 1:z
		D(jj) = R(:,jj)'*R(:,jj);    % Ensures differences are
					     % positive for comparison.
	end
	if mm < z
		
		sd = sign(D(mm) - D(mm+1)) ; 	% Compares for smallest
						% difference to find
	elseif mm > 1 				% new dividing points.
		sd = sign(D(mm-1) - D(mm));
	else
		sd = 0;
	end

	while D(mm) - D(mm+sd) > 0
		if (mm == 2) && (sd < 0)
			break;
		end

		if (mm == m-1) && (sd > 0)
			break;
		end
		mm = mm + sd;
	end
	nk(i)= mm+ js - 1; % knot positions.
end

dpkpc = nk; 	% knot positions or sequence.