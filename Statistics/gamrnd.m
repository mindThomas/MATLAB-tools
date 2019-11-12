function rnd = gamrnd(a,b,method)
% This function produces independent random variates from the Gamma distribution.
%
%  INPUTS
%    a       [double]    n*1 vector of positive parameters.
%    b       [double]    n*1 vector of positive parameters.
%    method  [string]    'BawensLubranoRichard' or anything else (see below).
%
%  OUTPUT
%    rnd     [double]    n*1 vector of independent variates from the gamma(a,b) distribution.
%                        rnd(i) is gamma distributed with mean a(i)b(i) and variance a(i)b(i)^2.
%
%  ALGORITHMS
%    Described in Bauwens, Lubrano and Richard (1999, page 316) and Devroye (1986, chapter 9).

% Copyright (C) 2006-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if (nargin < 2)
    error('gamrnd:: Two input arguments are needed!');
end

if nargin==2
    method= 'BauwensLubranoRichard';
    if any(a<1)
        method = 'Devroye';
        Devroye.small = 'Best'; % 'Weibull' , 'Johnk' , 'Berman' , 'GS' , 'Best'
                                % REMARK: The first algorithm (Weibull) is producing too much extreme values.
    end
    if ~strcmpi(method,'BauwensLubranoRichard')
        Devroye.big = 'Best'; % 'Cheng' , 'Best'
                              % REMARK 1: The first algorithm (Cheng) is still producing obviously wrong simulations.
                              % REMARK 2: The second algorithm seems slightly slower than the algorithm advocated by Bauwens,
                              %           Lubrano and Richard, but the comparison depends on the value of a (this should be
                              %           investigated further).
    end
else
    error('gamrnd:: Selection of method not yet implemented')
end

[ma,na] = size(a);
[mb,nb] = size(b);

if ma~=mb || na~=nb
    error('gamrnd:: Input arguments must have the same size!');
end

if na~=1
    error('gamrnd:: Input arguments must be column vectors');
end

if (any(a<0)) || (any(b<0)) || (any(a==Inf)) || (any(b==Inf))
    error('gamrnd:: Input arguments must be finite and positive!');
end

integers = isint(a);
doubles = (isnum(a) ~= isint(a));

integer_idx = find(integers);
double_idx = find(doubles);

number_of_integer_a = length(integer_idx);
number_of_double_a = length(double_idx);

rnd = NaN(ma,1);

if number_of_integer_a
    small_idx = find(a(integer_idx)<30);
    big_idx = find(a(integer_idx)>=30);
    number_of_small_a = length(small_idx);
    number_of_big_a = length(big_idx);
    if number_of_small_a
        % Exact sampling.
        for i=1:number_of_small_a
            rnd(integer_idx(small_idx(i))) = sum(exprnd(ones(a(integer_idx(small_idx(i))),1)))*b(integer_idx(small_idx(i)));
        end
    end
    if number_of_big_a
        % Gaussian approximation.
        rnd(integer_idx(big_idx)) = sqrt(a(integer_idx(big_idx))).* b(integer_idx(big_idx)) .* randn(number_of_big_a, 1) + a(integer_idx(big_idx)) .* b(integer_idx(big_idx));
    end
end


if number_of_double_a
    if strcmpi(method,'BauwensLubranoRichard')
        % Algorithm given in Bauwens, Lubrano & Richard (1999) page 316.
        rnd(double_idx) = knuth_algorithm(a(double_idx),b(double_idx));
    else% Algorithm given in  Devroye (1986, chapter 9)
        small_idx = find(a(double_idx)<1);
        big_idx = find(a(double_idx)>1);
        number_of_small_a = length(small_idx);
        number_of_big_a = length(big_idx);
        if number_of_small_a
            if strcmpi(Devroye.small,'Weibull')
                % Algorithm given in Devroye (1986, page 415) [Rejection from the Weibull density]
                rnd(double_idx(small_idx)) = weibull_rejection_algorithm(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'Johnk')
                % Algorithm given in Devroye (1986, page 418) [Johnk's gamma generator]
                rnd(double_idx(small_idx)) = johnk_algorithm(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'Berman')
                % Algorithm given in Devroye (1986, page 418) [Berman's gamma generator]
                rnd(double_idx(small_idx)) = berman_algorithm(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'GS')
                % Algorithm given in Devroye (1986, page 425) [Ahrens and Dieter, 1974]
                rnd(double_idx(small_idx)) = ahrens_dieter_algorithm(a(double_idx(small_idx)),b(double_idx(small_idx)));
            elseif strcmpi(Devroye.small,'Best')
                % Algorithm given in Devroye (1986, page 426) [Best, 1983]
                rnd(double_idx(small_idx)) = best_1983_algorithm(a(double_idx(small_idx)),b(double_idx(small_idx)));
            end
        end
        if number_of_big_a
            if strcmpi(Devroye.big,'Cheng')
                % Algorithm given in Devroye (1986, page 413) [Cheng's rejection algorithm GB]
                rnd(double_idx(big_idx)) = cheng_algorithm(a(double_idx(big_idx)),b(double_idx(big_idx)));
            elseif strcmpi(Devroye.big,'Best')
                % Algorithm given in Devroye (1986, page 410) [Best's rejection algorithm XG]
                rnd(double_idx(big_idx)) = best_1978_algorithm(a(double_idx(big_idx)),b(double_idx(big_idx)));
            end
        end
    end
end


function gamma_variates = weibull_rejection_algorithm(a,b)
nn = length(a);
mm = nn;
cc = 1./a ;
dd = a.^(a./(1-a)).*(1-a);
ZE = NaN(nn,2);
X  = NaN(nn,1);
INDEX = 1:mm;
index = INDEX;
while mm
    ZE(index,:) = exprnd(ones(mm,2));
    X(index) = ZE(index,1).^cc(index);
    id = find( (ZE(:,1)+ZE(:,2) > dd + X) );
    if isempty(id)
        mm = 0;
    else
        index = INDEX(id);
        mm = length(index);
    end
end
gamma_variates = X.*b;


function gamma_variates = johnk_algorithm(a,b)
nn = length(a);
mm = nn;
aa = 1./a ;
bb = 1./b ;
INDEX = 1:mm;
index = INDEX;
UV = NaN(nn,2);
X  = NaN(nn,1);
Y  = NaN(nn,1);
while mm
    UV(index,:) = rand(mm,2);
    X(index) = UV(index,1).^aa(index);
    Y(index) = UV(index,2).^bb(index);
    id = find(X+Y>1);
    if isempty(id)
        mm = 0;
    else
        index = INDEX(id);
        mm = length(index);
    end
end
gamma_variates = exprnd(ones(nn,1)).*X./(X+Y);


function gamma_variates = berman_algorithm(a,b)
nn = length(a);
mm = nn;
aa = 1./a ;
cc = 1./(1-a) ;
INDEX = 1:mm;
index = INDEX;
UV = NaN(nn,2);
X  = NaN(nn,1);
Y  = NaN(nn,1);
while mm
    UV(index,:) = rand(mm,2);
    X(index) = UV(index,1).^aa(index);
    Y(index) = UV(index,2).^cc(index);
    id = find(X+Y>1);
    if isempty(id)
        mm = 0;
    else
        index = INDEX(id);
        mm = length(index);
    end
end
Z = gamrnd(2*ones(nn,1),ones(nn,1));
gamma_variates = Z.*X.*b ;


function gamma_variates = ahrens_dieter_algorithm(a,b)
nn = length(a);
mm = nn;
bb = (exp(1)+a)/exp(1);
cc = 1./a;
INDEX = 1:mm;
index = INDEX;
UW = NaN(nn,2);
V  = NaN(nn,1);
X  = NaN(nn,1);
while mm
    UW(index,:) = rand(mm,2);
    V(index) = UW(index,1).*bb(index);
    state1 = find(V(index)<=1);
    state2 = find(V(index)>1);%setdiff(index,index(state1));
    ID = [];
    if ~isempty(state1)
        X(index(state1)) = V(index(state1)).^cc(index(state1));
        test1 = UW(index(state1),2) <= exp(-X(index(state1))) ;
        id1 = find(~test1);
        ID = INDEX(index(state1(id1)));
    end
    if ~isempty(state2)
        X(index(state2)) = -log(cc(index(state2)).*(bb(index(state2))-V(index(state2))));
        test2 = UW(index(state2),2) <= X(index(state2)).^(a(index(state2))-1);
        id2 = find(~test2);
        if isempty(ID)
            ID = INDEX(index(state2(id2)));
        else
            ID = [ID,INDEX(index(state2(id2)))];
        end
    end
    mm = length(ID);
    if mm
        index = ID;
    end
end
gamma_variates = X.*b ;


function gamma_variates = best_1983_algorithm(a,b)
nn = length(a);
mm = nn;
tt = .07 + .75*sqrt(1-a);
bb = 1 + exp(-tt).*a./tt;
cc = 1./a;
INDEX = 1:mm;
index = INDEX;
UW = NaN(nn,2);
V  = NaN(nn,1);
X  = NaN(nn,1);
Y  = NaN(nn,1);
while mm
    UW(index,:) = rand(mm,2);
    V(index) = UW(index,1).*bb(index);
    state1 = find(V(index)<=1);
    state2 = find(V(index)>1);%setdiff(index,index(state1));
    ID = [];
    if ~isempty(state1)
        X(index(state1)) = tt(index(state1)).*V(index(state1)).^cc(index(state1));
        test11 = UW(index(state1),2) <= (2-X(index(state1)))./(2+X(index(state1))) ;
        id11 = find(~test11);
        if ~isempty(id11)
            test12 = UW(index(state1(id11)),2) <= exp(-X(index(state1(id11)))) ;
            id12 = find(~test12);
        else
            id12 = [];
        end
        ID = INDEX(index(state1(id11(id12))));
    end
    if ~isempty(state2)
        X(index(state2)) = -log(cc(index(state2)).*tt(index(state2)).*(bb(index(state2))-V(index(state2)))) ;
        Y(index(state2)) = X(index(state2))./tt(index(state2)) ;
        test21 = UW(index(state2),2).*(a(index(state2)) + Y(index(state2)) - a(index(state2)).*Y(index(state2)) ) <= 1 ;
        id21 = find(~test21);
        if ~isempty(id21)
            test22 = UW(index(state2(id21)),2) <= Y(index(state2(id21))).^(a(index(state2(id21)))-1) ;
            id22 = find(~test22);
        else
            id22 = [];
        end
        if isempty(ID)
            ID = INDEX(index(state2(id21(id22))));
        else
            ID = [ID,INDEX(index(state2(id21(id22))))];
        end
    end
    mm = length(ID);
    if mm
        index = ID;
    end
end
gamma_variates = X.*b ;


function  gamma_variates = knuth_algorithm(a,b)
nn = length(a);
mm = nn;
bb = sqrt(2*a-1);
dd = 1./(a-1);
Y = NaN(nn,1);
X = NaN(nn,1);
INDEX = 1:mm;
index = INDEX;
while mm
    Y(index) = tan(pi*rand(mm,1));
    X(index) = Y(index).*bb(index) + a(index) - 1 ;
    idy1 = find(X(index)>=0);
    idn1 = setdiff(index,index(idy1));
    if ~isempty(idy1)
        test = log(rand(length(idy1),1)) <= ...
               log(1+Y(index(idy1)).*Y(index(idy1))) + ...
               (a(index(idy1))-1).*log(X(index(idy1)).*dd(index(idy1))) - ...
               Y(index(idy1)).*bb(index(idy1)) ;
        idy2 = find(test);
        idn2 = setdiff(idy1,idy1(idy2));
    else
        idy2 = [];
        idn2 = [];
    end
    index = [ INDEX(idn1) , INDEX(index(idn2)) ] ;
    mm = length(index);
end
gamma_variates = X.*b;


function  gamma_variates = cheng_algorithm(a,b)
nn = length(a);
mm = nn;
bb = a-log(4);
cc = a+sqrt(2*a-1);
UV = NaN(nn,2);
Y  = NaN(nn,1);
X  = NaN(nn,1);
Z  = NaN(nn,1);
R  = NaN(nn,1);
index = 1:nn;
INDEX = index;
while mm
    UV(index,:) = rand(mm,2);
    Y(index) = a(index).*log(UV(index,2)./(1-UV(index,2)));
    X(index) = a(index).*exp(UV(index,2));
    Z(index) = UV(index,1).*UV(index,2).*UV(index,2);
    R(index) = bb(index) + cc(index).*Y(index)-X(index);
    test1 = (R(index) >= 4.5*Z(index)-(1+log(4.5)));
    jndex = index(find(test1));
    Jndex = setdiff(index,jndex);
    if ~isempty(Jndex)
        test2 = (R(Jndex) >= log(Z(Jndex)));
        lndex = Jndex(find(test2));
        Lndex = setdiff(Jndex,lndex);
    else
        Lndex = [];
    end
    index = INDEX(Lndex);
    mm = length(index);
end
gamma_variates = X.*b;


function  gamma_variates = best_1978_algorithm(a,b)
nn = length(a);
mm = nn;
bb = a-1;
cc = 3*a-.75;
UV = NaN(nn,2);
Y  = NaN(nn,1);
X  = NaN(nn,1);
Z  = NaN(nn,1);
W  = NaN(nn,1);
index = 1:nn;
INDEX = index;
while mm
    UV(index,:) = rand(mm,2);
    W(index) = UV(index,1).*(1-UV(index,1));
    Y(index) = sqrt(cc(index)./W(index)).*(UV(index,1)-.5);
    X(index) = bb(index)+Y(index);
    jndex = index(find(X(index)>=0));
    Jndex = setdiff(index,jndex);
    if ~isempty(jndex)
        Z(jndex) = 64*W(jndex).*W(jndex).*W(jndex).*UV(jndex,2).*UV(jndex,2);
        kndex = jndex(find(Z(jndex)<=1-2*Y(jndex).*Y(jndex)./X(jndex)));
        Kndex = setdiff(jndex,kndex);
        if ~isempty(Kndex)
            lndex = Kndex(find(log(Z(Kndex))<=2*(bb(Kndex).*log(X(Kndex)./bb(Kndex))-Y(Kndex))));
            Lndex = setdiff(Kndex,lndex);
        else
            Lndex = [];
        end
        new_index = INDEX(Lndex);
        %mm = length(index);
    end
    index = union(new_index,INDEX(Jndex));
    mm = length(index);
end
gamma_variates = X.*b;
