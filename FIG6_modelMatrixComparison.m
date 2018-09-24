clear all; close all; clc;

load('Alorenz.mat');
load('Slorenz.mat');
dt = .001;
stackmax = 100;

figure;

subplot(2,2,1);
surf(A);
view(2);
axis('ij');
axis('equal');
axis('off');
title('SVD Coord. Operator', 'FontSize', 14);

subplot(2,2,2);
Arenorm = renormalize(A,S);
surf(Arenorm);
view(2);
axis('ij');
axis('equal');
axis('off');
title('Normalized SVD Coord. Operator', 'FontSize', 14);

QQ = buildLegendreOperator(length(A),stackmax*dt/2);

subplot(2,2,3);
surf(QQ);
view(2);
axis('ij');
axis('equal');
axis('off');
title('Legendre Coord. Operator', 'FontSize', 14);

subplot(2,2,4);
surf(renormalize(QQ,S));
view(2);
axis('ij');
axis('equal');
axis('off');
title('Normalized Legendre Coord. Operator', 'FontSize', 14);
caxis([min(min(Arenorm)), max(max(Arenorm))]);
