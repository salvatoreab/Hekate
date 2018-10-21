close all
clear all
clc

%% LEO

%Data
mu = 398600; %[km^3/s^2]

%Initial Conditions
y0 = [-6564, 0, 0, 0, -7.793, 0]; %[km, km, km, km/s, km/s, km/s]

%Integration time vector
tspan = linspace(0, 10000, 10000); %[s]

%Integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Integration
[~, Y] = ode113(@(t,y) Orbital_Lab1_fun(y,mu), tspan, y0, options);

%Position and velocity
r = Y(:,1:3);  %[km]
v = Y(:,4:6);  %[km/s]

%% Lunar injection

%Initial Conditions
y0 = [-6564, 0, 0, 0, -10.93, 0]; %[km, km, km, km/s, km/s, km/s]

%Integration time vector
tspan = linspace(0, 233266, 10000); %[s]

%Integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Integration
[~, P] = ode113(@(t,y) Orbital_Lab1_fun(y,mu), tspan, y0, options);

%Position and velocity
r1 = P(:,1:3);  %[km]
v1 = P(:,4:6);  %[km/s]

%% Hyperbola
%Data
mu = 4904; %[km^3/s^2]

%Initial Conditions
y0 = [-61241, 6082, 0, 0.5386, 0.0540, 0]; %[km, km, km, km/s, km/s, km/s]

%Integration time vector
tspan = linspace(0, 91419, 10000); %[s]

%Integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Integration
[~, Q] = ode113(@(t,y) Orbital_Lab1_fun(y,mu), tspan, y0, options);

%Position and velocity
r2 = Q(:,1:3);  %[km]
v2 = Q(:,4:6);  %[km/s]

for i = 1:length(tspan)
    rr2(i,:) = norm(r2(i,:));    
end
    [rp,rpi] = min(rr2);
    tspan(rpi);
    v2(rpi);

r2(:,1) = r2(:,1) + ones(length(r2(:,1)),1)*(382518.07+1.1e4);
r2(:,2) = r2(:,2) + ones(length(r2(:,2)),1)*(-37990.66-4.9e3);

%% Ellipse

%Data
mu = 4904; %[km^3/s^2]

%Initial Conditions
y0 = [3460.5, 2344.7, 0, 0.7461, -1.1022, 0]; %[km, km, km, km/s, km/s, km/s]

%Integration time vector
tspan = linspace(0, 35347, 10000); %[s]

%Integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Integration
[~, R] = ode113(@(t,y) Orbital_Lab1_fun(y,mu), tspan, y0, options);

%Position and velocity
r3 = R(:,1:3);  %[km]
v3 = R(:,4:6);  %[km/s]

for i = 1:length(tspan)
    rr3(i,:) = norm(r3(i,:));    
end
    [ra,rai] = max(rr3);
    tspan(rai);
    v3(rai)

r3(:,1) = r3(:,1) + ones(length(r3(:,1)),1)*(382518.07+1.1e4);
r3(:,2) = r3(:,2) + ones(length(r3(:,2)),1)*(-37990.66-4.9e3);

%% Plane change

%Data
mu = 4904; %[km^3/s^2]

%Initial Conditions
y0 = [-10672, -7212.2, 0, 0, 0, 0.308]; %[km, km, km, km/s, km/s, km/s]

%Integration time vector
tspan = linspace(0, 28309, 10000); %[s]

%Integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Integration
[~, S] = ode113(@(t,y) Orbital_Lab1_fun(y,mu), tspan, y0, options);

%Position and velocity
r4 = S(:,1:3);  %[km]
v4 = S(:,4:6);  %[km/s]

for i = 1:length(tspan)
    rr4(i,:) = norm(r4(i,:));    
end
    [rp4,rpi4] = min(rr4);
    tspan(rpi4);
    v4(rpi4);
    rp4

r4(:,1) = r4(:,1) + ones(length(r4(:,1)),1)*(382518.07+1.1e4);
r4(:,2) = r4(:,2) + ones(length(r4(:,2)),1)*(-37990.66-4.9e3);

%% Circularize

%Data
mu = 4904; %[km^3/s^2]

%Initial Conditions
y0 = [1518.7, 1026.4, -2.6874, 0, 0, -1.634]; %[km, km, km, km/s, km/s, km/s]

%Integration time vector
tspan = linspace(0, 8309, 10000); %[s]

%Integration options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Integration
[~, T] = ode113(@(t,y) Orbital_Lab1_fun(y,mu), tspan, y0, options);

%Position and velocity
r5 = T(:,1:3);  %[km]
v5 = T(:,4:6);  %[km/s]

for i = 1:length(tspan)
    rr5(i,:) = norm(r5(i,:));    
end
    [ra5,rai5] = max(rr5);
    tspan(rai5);
    v5(rai5);
    ra5

r5(:,1) = r5(:,1) + ones(length(r5(:,1)),1)*(382518.07+1.1e4);
r5(:,2) = r5(:,2) + ones(length(r5(:,2)),1)*(-37990.66-4.9e3);

%%
%Orbit representation
figure (1)
hold on
Terra3d
Luna3d
plot3(r(:,1),r(:,2),r(:,3), 'r', 'linewidth', 1.5)
plot3(r1(:,1),r1(:,2),r1(:,3), 'r', 'linewidth', 1.5)
plot3(r2(:,1),r2(:,2),r2(:,3), 'r', 'linewidth', 1.5)
plot3(r3(:,1),r3(:,2),r3(:,3), 'r', 'linewidth', 1.5)
plot3(r4(:,1),r4(:,2),r4(:,3), 'r', 'linewidth', 1.5)
plot3(r5(:,1),r5(:,2),r5(:,3), 'r', 'linewidth', 1.5)
grid on
hold off
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
title('ORBIT TRAJECTORY')
set(gca,'Color', [0.05, 0.05, 0.05])
ax = gca;
ax.GridColor = 'w';
ax.GridAlpha = 0.2;
