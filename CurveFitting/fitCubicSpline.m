function [S C] = fitCubicSpline(u,x,y,dya,dyb)
    % vectors x and y contain n+1 points and the corresponding function values
    % vector u contains all discrete samples of the continuous argument of f(x)
    % dya and dyb are the derivatives f'(x_0) and f'(x_n), respectively 
    n=length(x);       % number of interpolating points
    k=length(u);       % number of discrete sample points
    C=zeros(n,k);      % the n-1 cubic interpolating polynomials
    A=2*eye(n);        % coefficient matrix on left-hand side
    A(1,2)=1;
    A(n,n-1)=1;   
    d=zeros(n,1);      % vector on right-hand side
    d(1)=((y(2)-y(1))/(x(2)-x(1))-dya)/h0;  % first element of d
    for i=2:n-1
        h0=x(i)-x(i-1);
        h1=x(i+1)-x(i);
        h2=x(i+1)-x(i-1);       
        A(i,i-1)=h0/h2;
        A(i,i+1)=h1/h2;
        d(i)=((y(i+1)-y(i))/h1-(y(i)-y(i-1))/h0)/h2; % 2nd divided difference
    end
    d(n)=(dyb-(y(n)-y(n-1))/h1)/h1;   % last element of d
    M=6*inv(A)*d;                     % solving linear equation system for M's
    for i=2:n
        h=x(i)-x(i-1);
        x0=u-x(i-1);
        x1=x(i)-u;
        C(i-1,:)=(x1.^3*M(i-1)+x0.^3*M(i))/6/h... % the ith cubic polynomial
                 -(M(i-1)*x1+M(i)*x0)*h/6+(y(i-1)*x1+y(i)*x0)/h;  
        idx=find(u>x(i-1) & u<=x(i));  % indices between x(i-1) and x(i)
        S(idx)=C(i-1,idx);             % constructing spline by cubic polynomials
    end
end