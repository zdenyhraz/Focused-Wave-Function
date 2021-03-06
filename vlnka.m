clear all;
close all;
clc;

%============================ SET ============================
k=1*10^7;%1*10^7
f=0.25*10^(-1);%1*10^(-1)
a=1*10^(-2);%1*10^(-2)
dz=0.3*10^(-4);%3*10^(-4)
drho=1.0*10^(-5);%1.5*10^(-5)

N=201;
PSI=zeros(N,N);

N1D=1001;
PSI1D=zeros(N1D);

x = 1000.*linspace(-dz,+dz,N);%[mm]
y = 1000.*linspace(-drho,+drho,N);%[mm]
[X,Y] = meshgrid(x,y);

for idx=1:N    
    if (x(idx)==0)
        zfidx=idx;
    end
    
    if (y(idx)==0)
        rho0idx=idx;
    end
end

%============================ 1D ============================

for it=1:N1D
    z=f-dz+(it-1)/(N1D-1)*2*dz;
    zi(it)=1000*(z-f);%[mm]
    
    S1D1(it)=1/(-i*k/2*(a/f)^2*(1-3/4*(a/f)^2+5/8*(a/f)^4-35/64*(a/f)^6));
    S1D2(it)=z/(sqrt(z^2+a^2));
    S1D3(it)=exp(i*k*(sqrt(z^2+a^2)-sqrt(f^2+a^2)))/(sqrt(z^2+a^2)-sqrt(f^2+a^2));
    S1D4(it)=exp(i*k*(z-f))/(z-f);
    
    PSI1D(it)=S1D1(it)*(S1D2(it)*S1D3(it)-S1D4(it));
end

%============================ 2D ============================

parfor iz=1:N
    fprintf('progress=%.2f%%\n',iz/N*100)
    for irho=1:N
        z=f-dz+(iz-1)/(N-1)*2*dz;
        rho=0-drho+(irho-1)/(N-1)*2*drho;
        phi=0;
        
        S1=-i*k*z/(2*pi);
        S2=@(r) r.*(exp(-i*k*sqrt(f^2+r.^2)))/(sqrt(f^2+r.^2));
        S3=@(p,r) (exp(i*k*sqrt(z^2+rho^2+r.^2-2*rho*r.*cos(p-phi)))/(z^2+rho^2+r.^2-2*rho*r.*cos(p-phi)))*(1+i./(k*sqrt(z^2+rho^2+r.^2-2*rho*r.*cos(p-phi))));
        I1=@(r) integral(@(p) S3(p,r),0,2*pi,'AbsTol',1e-0,'ArrayValued',true);
        I2=integral(@(r) S2(r)*I1(r),0,a,'AbsTol',1e-0,'ArrayValued',true);
        
        PSI(irho,iz)=S1*I2;
    end
end

INTENSITY1D=abs(PSI1D).^2;
INTENSITY=abs(PSI).^2;
MAXINTENSITY=max(max(INTENSITY));
INTENSITY=INTENSITY./MAXINTENSITY;

INTENSITYAXZ=INTENSITY(rho0idx,:);
INTENSITYAXR=INTENSITY(:,zfidx);

fprintf('maxintens=%d\n',MAXINTENSITY);
fprintf('rho0idx=%d\n',rho0idx);
fprintf('zfidx=%d\n',zfidx);

%$\psi(rho,phi,z)$

%============================ FIGURES ============================
fntsz=14;

%1D inteenzita
figure('Position',[100 300 800 600])
plot(zi,INTENSITY1D);
title('$\psi(0,-,z)$ Osova relativni intenzita pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
xlim([-1000*dz,1000*dz])
ylim([0,1])

%rez podel z
figure('Position',[100 300 800 600])
plot(x,INTENSITYAXZ);
title('$\psi(0,-,z)$ Rez relativni intenzity ve smeru osy $z$ pro $\rho=0$ pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
xlim([-1000*dz,1000*dz])
ylim([0,1])

%rez podel z zoom
figure('Position',[100 300 800 600])
plot(x,INTENSITYAXZ);
title('$\psi(0,-,z)$ Rez relativni intenzity ve smeru osy $z$ pro $\rho=0$ pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
xlim([-1000*dz,1000*dz])
ylim([0,0.1])

%rez podel rho
figure('Position',[100 300 800 600])
plot(y,INTENSITYAXR);
title('$\psi(\rho,-,f)$ Rez relativni intenzity ve smeru osy $\rho$ pro $z=f$ pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$\rho$[mm]','Interpreter','latex')
xlim([-1000*drho,1000*drho])
ylim([0,1])

%rez podel rho zoom
figure('Position',[100 300 800 600])
plot(y,INTENSITYAXR);
title('$\psi(\rho,-,f)$ Rez relativni intenzity ve smeru osy $\rho$ pro $z=f$ pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$\rho$[mm]','Interpreter','latex')
xlim([-1000*drho,1000*drho])
ylim([0,0.05])

%contour intensity
figure('Position',[100 300 800 600])
contour(X,Y,INTENSITY,[0,0.0001,0.001,0.002,0.003,0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],'ShowText','on');
title('$\psi(\rho,-,z)$ Relativni intenzita pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet

%{

figure('Position',[100 300 800 600])
contourf(X,Y,abs(PSI).^2,100,'LineColor','none');
title('$\psi(\rho,-,z)$ Relativni intenzita pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar

figure('Position',[100 300 800 600])
surfc(X,Y,abs(PSI).^2)
title('$\psi(\rho,-,z)$ Relativni intenzita pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
shading interp
colormap jet
axis tight

figure('Position',[100 300 800 600])
plot(zi,angle(PSI1D));
title('$\psi(0,-,z)$ Osova faze pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(x,angle(PSIAXZ));
title('$\psi(0,-,z)$ Rez faze ve smeru osy $z$ pro $\rho=0$ pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(y,angle(PSIAXR));
title('$\psi(\rho,-,f)$ Rez faze ve smeru osy $\rho$ pro $z=f$ pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$\rho$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
contour(X,Y,angle(PSI),20);
title('$\psi(\rho,-,z)$ Faze pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar

figure('Position',[100 300 800 600])
contourf(X,Y,angle(PSI),100,'LineColor','none');
title('$\psi(\rho,-,z)$ Faze pro f='+compose("%.2g",f)+', a='+compose("%.2g",a)+', k='+compose("%.2g",k),'fontsize',fntsz,'Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar
 
%}
