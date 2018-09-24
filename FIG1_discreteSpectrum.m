clear all; close all; clc;

% Experiment parameters
tau = 1;
D = 100;
dt = 2*tau/D;
t = 0:dt:100;
N = 5;

% Build signal with random frequencies
frequencies = 10*rand(1,N);
N = length(frequencies);
amplitudes = rand(1,N) + 1;
phases = 2*pi*rand(1,N);

y = zeros(1,length(t));
for j = 1:N
    y = y + amplitudes(j)*sin((frequencies(j)*t+phases(j)));
end

% experiment!
[Ac,Ad, U, S, V, r, s, time] = fullstack_timedelay_model(y,D,dt,0,1i*[-frequencies, frequencies]);