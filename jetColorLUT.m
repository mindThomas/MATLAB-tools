lutSize = 1024;
staticLut = zeros(lutSize, 3);
jet = [0, 0, 1;
       0, 0.5, 0.5;
       0.8, 0.8, 0;
       1, 0, 0];

jet2 = [0, 0, 0.5;
   0, 0, 1;
   0, 1, 1;
   1, 1, 0;
   1, 0, 0
   0.5, 0, 0];
   
colorIndexMax = size(jet,1) - 2;

for (i = 0:(lutSize-1))
    znorm = i / (lutSize - 1);
    colorIndex = uint16(min(floor((size(jet,1) - 1) * znorm), colorIndexMax));
    alpha = (size(jet,1) - 1) * znorm - double(colorIndex);
    staticLut(i+1, :) = jet(colorIndex + 1, :) * (1-alpha) + jet(colorIndex + 2, :) * alpha;
end
staticLut = min(max(staticLut, 0), 1);

im = reshape(staticLut, [1, lutSize, 3]);
im = repmat(im, 100, 1);

figure(1);
imshow(im);

figure(2);
plot(staticLut(:,1), 'r');
hold on;
plot(staticLut(:,2), 'g');
plot(staticLut(:,3), 'b');