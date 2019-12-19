function [R,t]=registerSVD(P2,P1)
%Given two sets of points, it computes the least square minimization to
%find the best rotation R and translation t that would match the two
%shapes. It uses a singular value decomposition based method

% The article follow the algorithm described in: 
% Estimating 3-D rigid body transformations: a comparison
% of four major algorithms
% D.W. Eggert1 , A. Lorusso2 , R.B. Fisher3

% Find the mean of both point set, which will be their center of mass.

sizeN = size(P1,2);

mp1(1)=mean(P1(1,:));
mp1(2)=mean(P1(2,:));
mp1(3)=mean(P1(3,:));

mp2(1)=mean(P2(1,:));
mp2(2)=mean(P2(2,:));
mp2(3)=mean(P2(3,:));

cmp1 = P1-repmat(mp1',1,sizeN);
cmp2 = P2-repmat(mp2',1,sizeN);
H = zeros(3,3);
for i = 1:size(cmp1,2)
    H = H + cmp2(:,i)*cmp1(:,i)';
end

[U,S,V] = svd(H);

R = V*U';

t = mp1' - R*mp2';