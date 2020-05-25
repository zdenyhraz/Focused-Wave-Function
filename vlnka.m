clear all;
close all;
clc;

%============================ SET ============================
k=1*10^7;
f=1*10^(-1);
a=1*10^(-2);
dz=3*10^(-4);
drho=1.5*10^(-5);

N1D=1001;
N=1001;

PSI1D=zeros(N1D);
PSI=zeros(N,N);

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
    
    S1D1(it)=(2*i)/(k*((a/f)^2));
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
        
        S1=-(i*k*z)/(2*pi);
        S2=@(r) r.*(exp(-i*k*sqrt(f^2+r.^2)))/(sqrt(f^2+r.^2));
        S3=@(p,r) (exp(i*k*sqrt(z^2+rho^2+r.^2-2*rho*r.*cos(p-phi)))/(z^2+rho^2+r.^2-2*rho*r.*cos(p-phi)))*(1+i./(k*sqrt(z^2+rho^2+r.^2-2*rho*r.*cos(p-phi))));
        I1=@(r) integral(@(p) S3(p,r),0,2*pi,'AbsTol',1e-0,'ArrayValued',true);
        I2=integral(@(r) S2(r)*I1(r),0,a,'AbsTol',1e-0,'ArrayValued',true);
        
        PSI(irho,iz)=S1*I2;
    end
end

maxim=max(max(PSI));
PSI=PSI./maxim;
PSIAXZ=PSI(rho0idx,:);
PSIAXR=PSI(:,zfidx);

fprintf('maxim=%d\n',maxim);
fprintf('rho0idx=%d\n',rho0idx);
fprintf('zfidx=%d\n',zfidx);

figure('Position',[100 300 800 600])
plot(zi,abs(PSI1D).^2);
title('$\psi(0,-,z)$ Osov? Intenzita','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(zi,angle(PSI1D));
title('$\psi(0,-,z)$ Osov? F?ze','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(x,abs(PSIAXZ).^2);
title('$\psi(0,-,z)$ Intenzita z-Osov? cut','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(x,angle(PSIAXZ));
title('$\psi(0,-,z)$ F?ze z-Osov? cut','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(y,abs(PSIAXR).^2);
title('$\psi(\rho,-,f)$ Intenzita $\rho$-Osov? cut','Interpreter','latex')
xlabel('$\rho$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
plot(y,angle(PSIAXR));
title('$\psi(\rho,-,f)$ F?ze $\rho$-Osov? cut','Interpreter','latex')
xlabel('$\rho$[mm]','Interpreter','latex')

figure('Position',[100 300 800 600])
contour(X,Y,abs(PSI).^2,[0,0.001,0.002,0.003,0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],'ShowText','on');
title('$\psi(\rho,-,z)$ Intenzita','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar

figure('Position',[100 300 800 600])
contour(X,Y,angle(PSI),20);
title('$\psi(\rho,-,z)$ F?ze','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar

figure('Position',[100 300 800 600])
contourf(X,Y,abs(PSI).^2,100,'LineColor','none');
title('$\psi(\rho,-,z)$ Intenzita','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar

figure('Position',[100 300 800 600])
surfc(X,Y,abs(PSI).^2)
title('$\psi(\rho,-,z)$ Intenzita','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
shading interp
colormap jet
axis tight

figure('Position',[100 300 800 600])
contourf(X,Y,angle(PSI),100,'LineColor','none');
title('$\psi(\rho,-,z)$ F?ze','Interpreter','latex')
xlabel('$z-f$[mm]','Interpreter','latex')
ylabel('$\rho$[mm]','Interpreter','latex')
colormap jet
colorbar











