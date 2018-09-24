%% Simulation
clear all; close all; clc;

% Simulate
'Simulate'
type vanderpoldemo;
dt = 0.001;
n = 2;
Muvals = [.1, .5, 1.0, 3.0];

traja = [];
trajb = [];

for j = 1:4
    clear('t');
    Mu = Muvals(j);
    x0 = [2; 0];
    tspan=[dt:dt:200];
    N = length(tspan);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
    [t,xdat]=ode45(@(t,x) vanderpoldemo(t,x,Mu),tspan,x0,options);

    rmax = 8;
    lambda = 0;
    delta = 1;
    stackmax = 3000;

    t = t(1:end/2);
    xdat = xdat(end/2:end,:);
    traja(j,:) = xdat(:,1)';
    trajb(j,:) = xdat(:,2)';
end

%% Visualization
n = length(t);
omega = (2*pi/t(end))*[0:n/2 -n/2:-1];
y = traja(3,:);
stackmaxes = [200 1000 1500 3000];
thing = log(abs(fft(y)));
thing = thing/max(thing);

figure;
for j = 1:4
    stackmax = stackmaxes(j);
    H = Hankel(y,stackmax-1,1);
    [U,S,V] = svd(H,'econ');
    r = 24;
    A = buildbasismodel(U(:,1:r),dt,stackmax);

    subplot(2,2,j);
    plot(omega,thing,'b', 'Linewidth',1)
    hold on
    for lambda=eig(A)
        plot([imag(lambda),imag(lambda)],[0,1],'r', 'Linewidth',2);
    end
    m = max(abs(eig(A)));
    xlim([-m-1, m+1])
    ylim([-.2, 1.2])
end
for j = 1:4
    subplot(2,2,j);
    xlim([-m-1,m+1]);
    ylim([-.2,1.2]);
    xlabel('\omega');
    ylabel('Normalized FFT','Fontname','Palatino');
end