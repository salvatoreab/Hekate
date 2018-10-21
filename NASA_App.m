clear all
clc

m_earth = 5.97e24;
m_moon = 7.35e22;
R = 384400;
mu = 4.904e12/1e9;
r = R*(m_moon/(3*m_earth))^(1/3);
alpha0 = 0;
alpha1 = pi/2;
i = alpha1 - alpha0;

v = sqrt(mu/r);
dV = 2*v*sin(i/2);
