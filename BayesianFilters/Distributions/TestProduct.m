dist1 = Gaussian(0, 3);
dist2 = Gaussian(3, 3);

fig = figure(1);
dist1.plot(fig);
hold on;
dist2.plot(fig);

dist1.product(dist2).normalize().plot(fig)
hold off;