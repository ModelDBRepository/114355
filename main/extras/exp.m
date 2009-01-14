lambda = 50;

x = 0:0.0001:100;
y1 = lambda/100 * exp(-lambda/100*x);
y2 = lambda * exp(-lambda*x);
y2 = y2 / 100;

figure
hold on
plot(x,y1,'r');
plot(x,y2,'b');