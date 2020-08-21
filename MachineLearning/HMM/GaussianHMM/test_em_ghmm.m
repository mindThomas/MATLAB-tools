%%%%%% Example of Training/Testing a 2d-mixture of 2 gaussians driven by
%%%%%% HMM



d                                   = 2;
m                                   = 2;
L                                   = 1;
R                                   = 1;
Ntrain                              = 3000;
Ntest                               = 10000;
options.nb_ite                      = 30;

PI                                  = [0.5 ; 0.5];
A                                   = [0.95 0.05 ; 0.05 0.95];
M                                   = cat(3 , [-1 ; -1] , [2 ; 2]);
S                                   = cat(3 , [1 0.3 ; 0.3 0.8] , [0.7 0.6; 0.6 1]);

[Ztrain , Xtrain]                   = sample_ghmm(Ntrain , PI , A , M , S , L);
Xtrain                              = Xtrain - 1;

%%%%% initial parameters %%%%

PI0                                 = rand(d , 1 , R);
sumPI                               = sum(PI0);
PI0                                 = PI0./sumPI(ones(d , 1) , : , :);

A0                                  = rand(d , d , R);
sumA                                = sum(A0);
A0                                  = A0./sumA(ones(d , 1) , : , :);

M0                                  = randn(m , 1 , d , R);
S0                                  = repmat(cat(3 , [2 0 ; 0 2] , [3 0; 0 2]) , [1 , 1 , 1, R]);

%%%%% EM algorithm %%%%

[logl , PIest , Aest , Mest , Sest] = em_ghmm(Ztrain , PI0 , A0 , M0 , S0 , options);


[x , y]                             = ndellipse(M , S);
[xest , yest]                       = ndellipse(Mest , Sest);

Ltrain_est                          = likelihood_mvgm(Ztrain , Mest , Sest);
Xtrain_est                          = forward_backward(PIest , Aest , Ltrain_est);
Xtrain_est                          = Xtrain_est - 1;

ind1                                = (Xtrain_est == 0);
ind2                                = (Xtrain_est == 1);

Err_train                           = min(sum(Xtrain ~= Xtrain_est , 2)/Ntrain , sum(Xtrain ~= ~Xtrain_est , 2)/Ntrain);

figure(1) ,
h                                   = plot(Ztrain(1 , ind1) , Ztrain(2 , ind1) , 'k+' , Ztrain(1 , ind2) , Ztrain(2 , ind2) , 'g+' , x , y , 'b' , xest  , yest ,'r', 'linewidth' , 2);
legend([h(1) ; h(3:m:end)] , 'Train data' , 'True'  , 'Estimated' , 'location' , 'best')
title(sprintf('Train data, Error rate = %4.2f%%' , Err_train*100))

%%%%% Test data  %%%%


[Ztest , Xtest]                     = sample_ghmm(Ntest , PI , A , M , S , L);
Xtest                               = Xtest - 1;


Ltest_est                           = likelihood_mvgm(Ztest , Mest , Sest);
Xtest_est                           = forward_backward(PIest , Aest , Ltest_est);
Xtest_est                           = Xtest_est - 1;


ind1                                = (Xtest_est == 0);
ind2                                = (Xtest_est == 1);

Err_test                            = min(sum(Xtest ~= Xtest_est , 2)/Ntest , sum(Xtest ~= ~Xtest_est , 2)/Ntest);

figure(2),
h                                   = plot(Ztest(1 , ind1) , Ztest(2 , ind1) , 'k+' , Ztest(1 , ind2) , Ztest(2 , ind2) , 'g+' , x , y , 'b' , xest  , yest ,'r', 'linewidth' , 2);
legend([h(1) ; h(3:m:end)] , 'Test data' , 'True'  , 'Estimated' , 'location' , 'best')
title(sprintf('Test data, Error rate = %4.2f%%' , Err_test*100))