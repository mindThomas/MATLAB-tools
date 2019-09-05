function [error,Reallignedsource]=ICPmanu_allign2(target,source)

[IDX1(:,1),IDX1(:,2)]=knnsearch(target,source);
[IDX2(:,1),IDX2(:,2)]=knnsearch(source,target);
IDX1(:,3)=1:length(source(:,1));
IDX2(:,3)=1:length(target(:,1));

m1=mean(IDX1(:,2));
s1=std(IDX1(:,2));
IDX1=IDX1(IDX1(:,2)<(m1+1.65*s1),:);

m2=mean(IDX2(:,2));
s2=std(IDX2(:,2));
IDX2=IDX2(IDX2(:,2)<(m2+1.65*s2),:);

Datasetsource=vertcat(source(IDX1(:,3),:),source(IDX2(:,1),:));
Datasettarget=vertcat(target(IDX1(:,1),:),target(IDX2(:,3),:));

[error,Reallignedsource,transform] = procrustes(Datasettarget,Datasetsource);
Reallignedsource=transform.b*source*transform.T+repmat(transform.c(1,1:3),size(source,1),1);


