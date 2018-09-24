clear all; close all; clc;

% Simulate
'Simulate'
type vanderpoldemo;
dt = 0.01;
n = 2;
Mu = 1.0;

x0 = [4; 4];
tspan=[dt:dt:200];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,xdat]=ode45(@(t,x) vanderpoldemo(t,x,Mu),tspan,x0,options);

r = 28;
stackmax = round(8/dt);
Npoly = 6;

compareEdmdHAVOK(xdat',dt,stackmax,Npoly,r);