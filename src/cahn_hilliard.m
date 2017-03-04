clear all; close all; clc;
imax=128; jmax=128;
hx=1; hy=1;
dx=hx; dy=hy;
a=10; b=1;
dt=2*10^-3;
tmax=500;
t=tmax/dt;
m=0:dt:tmax-dt;
%Initialization
phi=[0.49,0.51];
phi=(phi(2)-phi(1)).*rand(128,128)+phi(1);
% mesh(phi); contour(phi);

%Running the simulation
for k=1:t;
    
    %Finding out the Laplacians of phi
    for j=1:jmax
        for i=2:imax-1
            d2phi_x(i,j)=(phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/dx^2;
        end
        d2phi_x(1,j)=(phi(2,j)-2*phi(1,j)+phi(imax,j))/(dx^2);
        d2phi_x(imax,j)=(phi(1,j)-2*phi(imax,j)+phi(imax-1,j))/dx^2;
    end
    
    for i=1:imax
        for j=2:jmax-1
            d2phi_y(i,j)=(phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/dy^2;
        end
        d2phi_y(i,1)=(phi(i,2)-2*phi(i,1)+phi(i,jmax))/dy^2;
        d2phi_y(i,jmax)=(phi(i,1)-2*phi(i,jmax)+phi(i,jmax-1))/dy^2;
    end
    
    %Finding out the g'(phi) values
    for j=1:jmax
        for i=1:imax
            gdashphi(i,j)=32*b*(2*phi(i,j).^3-3*phi(i,j).^2+phi(i,j));
        end
    end
    
    %Finding value of mu
    for j=1:jmax
        for i=1:imax
            mu(i,j)=gdashphi(i,j)-a.*(d2phi_x(i,j)+d2phi_y(i,j));
        end
    end
    
    %Finding laplacian of mu
    for j=1:jmax
        for i=2:imax-1
            d2mu_x(i,j)=(mu(i+1,j)-2*mu(i,j)+mu(i-1,j))/dx^2;
        end
        d2mu_x(1,j)=(mu(2,j)-2*mu(1,j)+mu(imax,j))/dx^2;
        d2mu_x(imax,j)=(mu(1,j)-2*mu(imax,j)+mu(imax-1,j))/dx^2;
    end
    
    for i=1:imax
        for j=2:jmax-1
            d2mu_y(i,j)=(mu(i,j+1)-2*mu(i,j)+mu(i,j-1))/dy^2;
        end
        d2mu_y(i,1)=(mu(i,2)-2*mu(i,1)+mu(i,jmax))/dy^2;
        d2mu_y(i,jmax)=(mu(i,1)-2*mu(i,jmax)+mu(i,jmax-1))/dy^2;
    end
    
    %Finding out the new phi values
    for j=1:jmax
        for i=1:imax
            new_phi(i,j)=phi(i,j)+dt.*(d2mu_x(i,j)+d2mu_y(i,j));
        end
    end
    
    %Calculation of Interfacial energy, Mixing Energy and Total Energy
    
    %Finding out the divergence of phi
    for j=1:jmax
        for i=2:imax-1
            dphi_x(i,j)=(phi(i+1,j)-phi(i-1,j))/(2*dx);
        end
        dphi_x(1,j)=(phi(2,j)-phi(imax,j))/(2*dx);
        dphi_x(imax,j)=(phi(1,j)-phi(imax-1,j))/(2*dx);
    end
    
    for i=1:imax
        for j=2:jmax-1
            dphi_y(i,j)=(phi(i,j+1)-phi(i,j-1))/(2*dy);
        end
        dphi_y(i,1)=(phi(i,2)-phi(i,jmax))/(2*dy);
        dphi_y(i,jmax)=(phi(i,1)-phi(i,jmax-1))/(2*dy);
    end
    
    for j=1:jmax
        for i=1:imax
            div_phi(i,j)=dphi_x(i,j)+dphi_y(i,j);
        end
    end
    IE(k)=0; ME(k)=0; TE(k)=0; int_area(k)=0; int_length(k)=0;
    
    for j=1:jmax
        for i=1:imax
            IE(k)=IE(k)+0.5*a.*(div_phi(i,j)*div_phi(i,j));
            ME(k)=ME(k)+16*b.*((phi(i,j)-1).^2).*(phi(i,j)).^2;
            TE(k)=TE(k)+0.5*a.*(div_phi(i,j)*div_phi(i,j))+...
                16*b.*((phi(i,j)-1).^2).*(phi(i,j)).^2;
            int_area(k)=int_area(k)+phi(i,j);
            int_length(k)=int_length(k)+abs(sqrt(div_phi(i,j)));
        end
    end
    
    
    %Characteristic Length
    rc(k)=int_area(k)/int_length(k);
 
    phi=new_phi;
    %Plotting the results for phase field 
    figure(1);
    if k*dt==0.002
        subplot(2,3,1)
        axis square
        [hC hC]=contourf(phi);
        set(hC,'LineStyle','none');
        daspect([1 1 1]);
        xlabel('(a) t=0');
        title(k*dt);
    end
    if k*dt==4
        subplot(2,3,2)
        axis square
        [hC hC]=contourf(phi);
        set(hC,'LineStyle','none');
        daspect([1 1 1]);
        xlabel('(b) t=4');
        title(k*dt);
    end
    if k*dt==10
        subplot(2,3,3)
        axis square
        [hC hC]=contourf(phi);
        set(hC,'LineStyle','none');
        daspect([1 1 1]);
        xlabel('(c) t=10');
        title(k*dt);
    end
    if k*dt==50
        subplot(2,3,4)
        axis square
        [hC hC]=contourf(phi);
        set(hC,'LineStyle','none');
        daspect([1 1 1]);
        xlabel('(d) t=50');
        title(k*dt);
    end
    if k*dt==200
        subplot(2,3,5)
        axis square
        [hC hC]=contourf(phi);
        set(hC,'LineStyle','none');
        daspect([1 1 1]);
        xlabel('(e) t=200');
        title(k*dt);
    end
    if k*dt==500
        subplot(2,3,6)
        axis square
        [hC hC]=contourf(phi);
        set(hC,'LineStyle','none');
        daspect([1 1 1]);
        xlabel('(f) t=500');
        title(k*dt);
    end
end

figure(2);
loglog(m,IE,'-k','LineWidth',2);
hold on
loglog(m,ME,'-.b');
loglog(m,TE,'-.r');
%Including t^(-1/3) and t^(-1/4) lines in the plot
q=4:dt:tmax;
t_4=(q).^(-1/4)+10000;
t_3=(q).^(-1/3)+1000;
loglog(q,t_4,'-c');
loglog(q,t_3,'m');
title('Plot for Energy and Energy Rate decrease');
hold off;

figure(3);
loglog(m,rc,'-k','LineWidth',2);
hold on
w=4:dt:tmax;
tc_4=(w).^(1/4);
tc_3=(w).^(1/3);
loglog(w,tc_4,'-b');
loglog(w,tc_3,'-r');
title('Plot for Characteristic Length and growth');
hold off;