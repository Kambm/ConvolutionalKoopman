function [usol,dt] = nls(slices,T)
  % space
  L=30; n=256;
  x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
  k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
  % time
  t=linspace(0,T,slices+1); dt=t(2)-t(1); 

  nf=1;
  slicesf=slices*nf;
  tf=linspace(0,2*pi*nf,slicesf+1);
  
 % initial conditions
  N=2;
  u=N*(sech(x)).';
  ut=fft(u);
  [t,utsol]=ode45('dmd_soliton_rhs',t,ut,[],k);
  for j=1:length(t)
      usol(j,:)=ifft(utsol(j,:));  % bring back to space
  end
end
