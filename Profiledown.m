clear all
clc

t = zeros(1,11);
y(1) = 108e3;
y(2) = 98e3; 
x(1) = 0;
vy(1) = 0;
vx(1) = 1.63e3;
g_m = -9.81/6;
tvect = zeros(1,11);
vxx(1) = vx(1);
vyy(1) = vy(1);
xx(1) = x(1);
yy(1) = y(1);
for i=1:10
    tvect(i+1) = -vy(i)/g_m + sqrt((vy(i)/g_m)^2-2/g_m*(y(i)-y(i+1)));
    vy(i+1) = vy(i) + g_m*tvect(i+1);
    x(i+1) = x(i) + vx(i)*tvect(i+1);
    vx(i+1) = vx(i);
    %Before thrust
    vx(i+1) = vx(i+1) - 0.163e3;
    vy(i+1) = vy(i+1) + 0.0566e3;
    %After thrust
    vxx(i+1) = vx(i+1);
    vyy(i+1) = vy(i+1);
    yy(i+1) = y(i+1);
    xx(i+1) = x(i+1);
    t(i+1) = t(i) + tvect(i+1);
    %
    y(i+2) = y(i+1)-10e3;
    y(i+1) = y(i)-10e3;
    x(i+2) = x(i+1);
   
end

figure(1)
plot(xx, yy)

figure(2)
plot(t,xx)

figure(3)
plot(t,yy)

figure(4)
plot(t,vxx)

figure(5)
plot(t,vyy)