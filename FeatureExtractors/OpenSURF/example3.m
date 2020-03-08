% Example 3, Affine registration
% Load images
I1=im2double(imread('TestImages/lena1.png'));
I2=im2double(imread('TestImages/lena2.png'));

% Get the Key Points
Options.upright=true;
Options.tresh=0.0001;
Ipts1=OpenSurf(I1,Options);
Ipts2=OpenSurf(I2,Options);

% Put the landmark descriptors in a matrix
D1 = reshape([Ipts1.descriptor],64,[]);
D2 = reshape([Ipts2.descriptor],64,[]);

% Find the best matches
err=zeros(1,length(Ipts1));
cor1=1:length(Ipts1);
cor2=zeros(1,length(Ipts1));
for i=1:length(Ipts1),
    distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
    [err(i),cor2(i)]=min(distance);
end

% Sort matches on vector distance
[err, ind]=sort(err);
cor1=cor1(ind);
cor2=cor2(ind);

% Make vectors with the coordinates of the best matches
Pos1=[[Ipts1(cor1).y]',[Ipts1(cor1).x]'];
Pos2=[[Ipts2(cor2).y]',[Ipts2(cor2).x]'];
Pos1=Pos1(1:30,:);
Pos2=Pos2(1:30,:);

% Show both images
I = zeros([size(I1,1) size(I1,2)*2 size(I1,3)]);
I(:,1:size(I1,2),:)=I1; I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=I2;
figure, imshow(I); hold on;

% Show the best matches
plot([Pos1(:,2) Pos2(:,2)+size(I1,2)]',[Pos1(:,1) Pos2(:,1)]','-');
plot([Pos1(:,2) Pos2(:,2)+size(I1,2)]',[Pos1(:,1) Pos2(:,1)]','o');

% Calculate affine matrix
Pos1(:,3)=1; Pos2(:,3)=1;
M=Pos1'/Pos2';

% Add subfunctions to Matlab Search path
functionname='OpenSurf.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir '/WarpFunctions'])

% Warp the image
I1_warped=affine_warp(I1,M,'bicubic');

% Show the result
figure,
subplot(1,3,1), imshow(I1);title('Figure 1');
subplot(1,3,2), imshow(I2);title('Figure 2');
subplot(1,3,3), imshow(I1_warped);title('Warped Figure 1');


