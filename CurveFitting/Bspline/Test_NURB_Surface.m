%% Test basis function
surface = NURB_PixelSurface(0);

%u = min(spline.knots):0.02:max(spline.knots);
u = 0:0.02:10;
B = zeros(size(u));
dB = zeros(size(u));
ddB = zeros(size(u));
dddB = zeros(size(u));
for (i = 1:length(u))    
    B(i) = surface.basis(5, u(i));    
    dB(i) = surface.dbasis(5, u(i));    
    ddB(i) = surface.ddbasis(5, u(i));  
    dddB(i) = surface.dddbasis(5, u(i));  
end

figure(1);
plot(u, B);
hold on;
plot(u, dB);
plot(u, ddB);
plot(u, dddB);
hold off;

legend('3rd order basis', '1st derivative', '2nd derivative', '3rd derivative');

%%
image = imread('map.png');
image = image(:,:,1);
%image = image(230:260, 30:70);
image = image(1:1:end, 1:1:end);

imshow(image);

%%
surface = NURB_PixelSurface(image);

[X,Y] = meshgrid(1:0.1:size(image,2), 1:0.1:size(image,1));
S = zeros(size(X));

for (i = 1:size(X,1))
    for (j = 1:size(X,2))
        x = X(i,j);
        y = Y(i,j);
        S(i,j) = surface.eval(x,y);
    end
end

figure(2);
surf(X, Y, S);
shading flat