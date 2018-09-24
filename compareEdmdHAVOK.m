function [Xh, Xdmd, U,S,V,r,omega,AH,yH] = compareEdmdHAVOK(xdat,dt,Ndelay,Npoly,r)
    close all; clc;
    %% Build HAVOK model
    
    % build model
    H = Hankel(xdat,Ndelay-1,1);
    [U,S,V] = svd(H,'econ');
    AH = V(1:end-1,1:r)\V(2:end,1:r);
    AH = renormalize(AH,S);
    AH = logm(AH)/dt;
    AH = AH';
    'HAVOK model built'
    
    % Havok reconstruction
    Vn = U'*H;
    size(AH)
    size(Vn(1:r,1))
    L = 1:length(V);
    time = dt*(L-1);
    sys = ss(AH,0*AH(:,r),eye(size(AH)),0*AH(:,r));
    [yH,t] = lsim(sys,0*Vn(r,L),time,Vn(1:r,1));
    size(yH)
    
    Xh = U(:,1:r)*yH';
    figure(1);
    subplot(2,1,1);
    set(gca,'Fontname','Palatino');
    set(gca,'Fontsize',16);
    hold on;
    plot(t,xdat(2,L),'k','Linewidth',2);
    plot(t,Xh(Ndelay+1,L),'r','Linewidth',2);
    xlabel('t');
    ylabel('x(t)');
    title('Time Delay Observables (\tau = 8)')
    
    'HAVOK Reconstruction'
    
    %% Edmd
    
    X = poolData(xdat',size(xdat,1),Npoly,0)';
    size(X)
    [Phi, omega, lambda, b, Xdmd] = DMDfull(X, {'dt',dt});
    'Edmd model built'
    size(Xdmd)
    L = 1:size(Xdmd,2);
    time = (L-1)*dt;
    figure(1);
    subplot(2,1,2);
    set(gca,'Fontname','Palatino');
    set(gca,'Fontsize',16);
    hold on;
    plot(time,xdat(2,L),'k','Linewidth',2);
    plot(time,Xdmd(3,L),'r','Linewidth',2);
    xlabel('t');
    ylabel('x(t)');
    title('Polynomial Observables (n = 6)');
    
    'Edmd reconstruction'

    
end