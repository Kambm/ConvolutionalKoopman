function [Ac,Ad, U, S, V, r, s, time] = fullstack_timedelay_model(x,N,dt,forcing,varargin)

    % Build Hankel Matrix and compute SVD coordinates
    H = Hankel(x,N-1,1);
    [U,S,V] = svd(H,'econ');  
    
    % Compute Rank Truncation
    r = truncationRank(S);
    
    Vn = U'*H;

    % Compute Linear Models
    
    % continuous time basis model
    Ac = buildbasismodel(U(:,1:r),dt,N);

    % discrete time basis model
    Ad = expm(Ac*dt);
    % discrete time DMD model
%     ADMD = (Vn(1:r,1:end-1)')\(Vn(1:r,2:end)');
%     % continuous time DMD model
%     ADMDc = logm(ADMD)/dt;
    [Phi, omega, lambda, b, Xdmd] = DMDfull(Vn(1:r,:),{'dt',dt});
    size(Xdmd)
    
    % Compare Spectra
    figure;
    set(gca,'FontSize',16);
    set(gca,'FontName','Times');
    hold on;
    scatter(real(eig(Ac)),imag(eig(Ac)),200,'*','Linewidth',1);
    nargin
    if nargin > 4
        length(varargin{1})
        length(eig(Ac))
        scatter(real(varargin{1}),imag(varargin{1}),200,'d','Linewidth',1);
        legend('Estimated Spectrum', 'True Spectrum');
    else
        legend('Estimated Spectrum');
    end
    title('Comparison of Spectra','Fontname','Times');
    ylabel('\Im\{\lambda\}');
    xlabel('\Re\{\lambda\}');
    axis('equal');
    box on;
    
    % Simulation and reconstruction
    
    maxrange = 5e+4;
    L = 1:min(min(length(V),length(Xdmd)),maxrange)
    time = dt*(L-1);
    if forcing
        Af = Ac(1:r-1,1:r-1);
        Bf = Ac(1:r-1,r);
        sys = ss(Af,Bf,eye(r-1),0*Bf);
        [s,t] = lsim(sys,Vn(r,L),time,Vn(1:r-1,1));
    else
        sys = ss(Ac,0*Ac(:,r),eye(size(Ac)),0*Ac(:,r));
        [s,t] = lsim(sys,0*Vn(r,L),time,Vn(1:r,1)');
    end
        
    figure;
    set(gca,'FontSize',16);
    set(gca,'FontName','Times');
    hold on;
    plot(time,Vn(1,L),'--k','Linewidth',2);
    plot(t,s(:,1),'.r','Linewidth',1);
    title('First Principal Component','Fontname','Times');    
    xlabel('Time');
    ylabel('\sigma_1 v_1');
    legend('Truth','Simulated (Basis Model)');
    box on;
    
    figure;
    hold on;
    set(gca,'FontSize',16);
    set(gca,'FontName','Times');
    
    [VV,WW] = eig(Ac);
    eigvecs = U(:,1:length(VV))*VV;
    size(eigvecs)
    p = ((1:N) - N/2)*dt;
    size(p)
    for j = 1:size(eigvecs,2)
        plot(p,real(eigvecs(:,j)),'linewidth',2);
    end
    xlabel('p');
    ylabel('Re\{\phi(p)\}');
    title('Real Part of Time Delay Eigenvectors','Fontname','Times');
end