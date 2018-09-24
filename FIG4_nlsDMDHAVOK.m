%% Params
clear all; close all; clc;

T = 16*pi;
slices = 2000;
rmin = 1;
rmax = 20;
rstep = 1;
rplot = 14; % plots the reconstructions at r = 14

%% SOlVE NLS Eqn

[usol,dt] = nls(slices,T);

X = usol.';
Xtrain = X(:,1:end/2);
Xtest = X(:,end/2+1:end);

%% Iterate over reconstructions
dmdtrainerrors = [];
dmdtesterrors = [];
edmdtrainerrors = [];
edmdtesterrors = [];
havoktrainerrors = [];
havoktesterrors = [];

for r = rmin:rstep:rmax
    r
    % DMD
    [Phi, omega, lambda, b, Xdmd] = DMDfull(Xtrain,{'dt',dt,'r',r});
    Xdmd = dmdReconstruct(Phi,omega,X(:,1),slices+1,dt);
    dmdtrainerrors(end+1) = norm(abs(Xdmd(:,1:end/2)-Xtrain))/(size(Xtrain,1)*size(Xtrain,2));
    dmdtesterrors(end+1) = norm(abs(Xdmd(:,end/2+1:end)-Xtest))/(size(Xtest,1)*size(Xtest,2));
    
    % eDMD
    Xaug = [Xtrain; Xtrain.*abs(Xtrain).^2];
    [PhiE, omegaE, lambdaE, bE, XEdmd] = DMDfull(Xaug,{'dt',dt,'r',r});
    XEdmd = dmdReconstruct(PhiE,omegaE,Xaug(:,1),slices+1,dt);
    XEdmd = XEdmd(1:size(Xtrain,1),:);
    edmdtrainerrors(end+1) = norm(abs(XEdmd(:,1:end/2)-Xtrain))/(size(Xtrain,1)*size(Xtrain,2));
    edmdtesterrors(end+1) = norm(abs(XEdmd(:,end/2+1:end)-Xtest))/(size(Xtest,1)*size(Xtest,2));
    
    % HAVOK
    delay = 10;
    H = Hankel(Xtrain,delay-1,1);
    [U,S,V] = svd(H,'econ');
    A = buildbasismodel(U(:,1:r),dt,delay);
    [PhiV,omegav] = eig(A);

    xhavok = dmdReconstruct(PhiV,diag(omegav),S(1:r,1:r)*V(1,1:r)',slices+1-delay,dt);
    Hhavok = U(:,1:r)*xhavok;
    Xhavok = unHankel(Hhavok,delay);
    
    havoktrainerrors(end+1) = norm(abs(Xhavok(:,1:end/2)-Xtrain))/(size(Xtrain,1)*size(Xtrain,2));
    havoktesterrors(end+1) = norm(abs(Xhavok(:,end/2+1:end)-Xtest))/(size(Xtest,1)*size(Xtest,2));
    if r == rplot
        figure(1);
        subplot(3,3,1);
        ylabel('$|u|$');
        surf(abs(Xtest));
        title('Truth','Fontname','Palatino','Fontsize',16);
        shading interp;
        
        
        subplot(3,3,2);
        surf(abs(Xdmd(:,end/2+1:end)));
        title('DMD Prediction','Fontname','Palatino','Fontsize',16);
        shading interp;
        
        subplot(3,3,3);
        surf(abs(Xhavok(:,end/2+1:end)));
        title('Havok Prediction','Fontname','Palatino','Fontsize',16);
        shading interp;
    end
end

figure(1);
subplot(3,3,[4:6]);
set(gca,'Fontsize',16);
hold on;
plot(rmin:rstep:rmax,log(dmdtrainerrors),'k','Linewidth',1);
plot(rmin:rstep:rmax,log(edmdtrainerrors),'--b','Linewidth',1);
plot(rmin:rstep:rmax,log(havoktrainerrors),'.-.r','Linewidth',1);
xlabel('Truncation rank','interpreter','latex');
ylabel('Log RMS Error','interpreter','latex');
title('Training Errors','Fontname','Palatino');
legend('DMD training error','eDMD training error','HAVOK training error');

subplot(3,3,[7:9]);
set(gca,'Fontsize',16);
hold on;
plot(rmin:rstep:rmax,log(dmdtesterrors),'k','Linewidth',1);
plot(rmin:rstep:rmax,log(edmdtesterrors),'--b','Linewidth',1);
plot(rmin:rstep:rmax,log(havoktesterrors),'.-.r','Linewidth',1);
xlabel('Truncation rank','interpreter','latex');
ylabel('Log RMS Error','interpreter','latex');
title('Test Errors','Fontname','Palatino');
legend('DMD test error','eDMD test error','HAVOK test error');


