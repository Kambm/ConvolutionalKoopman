clear all; close all; clc;

load('Vlorenz.mat');
load('xdatLorenz.mat');

figure;
coords = [1 2 3 5 7 10];
rows = 2;
for j = 1:length(coords)
    subplot(rows,length(coords)/rows,j);
    set(gca,'Fontsize',16);
    set(gca,'Fontname','Palatino');
    axis off;
    plotcolor(xdat(1050:end-50,:),10*V(1000:end,coords(j)));
    view([0 0]);
    str = '$$ w _{'+string(coords(j))+'} $$';
    title(str,'Interpreter','latex');
end