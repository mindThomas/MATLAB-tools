% TEST  --  Find Closest Point
%
% Given a line that goes from point A to point B, find the point D on the
% line AB that is closest to another point C.
%

% Construct a random set of lines and points
nRow = 2;
nCol = 3;
n = nRow*nCol;
A = rand(2,n);
B = rand(2,n);
C = rand(2,n);

% Compute the closest points
[D,lambda,isConvex] = findClosestPoint(A,B,C);

% Plot the solution:
figure(1); clf;
for i=1:nRow
    for j=1:nCol
        idx = nCol*(i-1) + j;
        subplot(nRow,nCol,idx); hold on;
        plot([A(1,idx), B(1,idx)],[A(2,idx), B(2,idx)],'k-','LineWidth',2);
        plot([C(1,idx), D(1,idx)],[C(2,idx), D(2,idx)],'k--','LineWidth',1);
        plot(C(1,idx), C(2,idx),'go','MarkerSize',10,'LineWidth',3);
        plot(D(1,idx), D(2,idx),'rs','MarkerSize',10,'LineWidth',3);
        axis equal;
        title(['Test: ' num2str(idx)]);
        xlabel('x')
        ylabel('y')
    end
end



