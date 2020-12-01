%Piotr Liebersbach
%ME 573 - CFD
%Homework 12
%12/13/19
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lid Driven Cavity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
fontsize = 18;

y_p=[1.0000 0.9766  0.9688 0.9609  0.9531 0.8516  0.7344 0.6172 0.5000 ...
     0.4531 0.2813 0.1719  0.1016 0.0703 0.0625 0.0547 0.0000];
%y coordinate
u_re100=[1.0000 0.8412 0.7887 0.7372 0.68717 0.2315 0.0033  -0.1364  -0.2058...
        -0.2109  -0.1566  -0.1015  -0.0643  -0.04775  -0.0419  -0.0371 0.0000];

x_p=[1.0000 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5000 ...
     0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0.0000];
v_re100=[0.0000 -0.05906  -0.0739 -0.0886 -0.10313 -0.16914 -0.22445 ...
  -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.1009 0.0923 0.0000];

  %initialize parameters
nu = 0.01;
L = 1.0;
rho = 1.0;
U_lid = 1.0;
dx = 0.05;
dy = dx;
gamma = 0;
C = 1;
dt = C*min((0.25*(dx^2)/nu) , (dx/U_lid));
w = 1.6;
Tolfac = 10^-6;
t_end = 15
%t_end = 80*dt
% t_end = 2*dt

%x - grid nodes
xu = (0:dx:L);
yu = (-0.5*dy:dy:L+0.5*dy);
[XU,YU] = meshgrid(xu,yu);
XU = XU'; YU = YU';
imax = length(xu);
Nyu = length(yu);

%y - grid nodes
xv = (-0.5*dx:dx:L+0.5*dx);
yv = (0:dy:L);
[XV, YV] = meshgrid(xv,yv);
XV = XV'; YV = YV';
Nxv = length(xv);
jmax = length(yv);

%P - grid nodes (center nodes)
xp = (-0.5*dx:dx:L+0.5*dx);
yp = (-0.5*dy:dy:L+0.5*dy);
[XP, YP] =  meshgrid(xp,yp);
XP = XP'; YP = YP';


%initialize  u and v fields, initialize derivatives
UU = zeros(size(XU));
VV = zeros(size(XV));
PP = zeros(size(XP));
UU_P = zeros(imax+1,jmax+1);
VV_P = UU_P;
Velocity_mag = UU_P;

F = UU;
G = VV;
g = PP;

AA = zeros(size(XU));
BB = AA;
%CC = AA;
DD = AA;
%EE = AA;

FF = zeros(size(XV));
GG = FF;
%HH = FF;
II = FF;
%JJ = FF;

%Set bc's on boundary
UU(1,:) = 0;
UU(imax,:) = 0;
VV(:,jmax) = 0;
VV(:,1) = 0;

nn = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Time Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = dt:dt:t_end
  fprintf('%f\n', n)
  VV(1,:) = -VV(2,:);
  VV(imax+1,:) = -VV(imax, :);
  UU(:,jmax+1) = 2*U_lid - UU(:,jmax);
  UU(:,1) = -UU(:,2);
  
  %calculate intermediate u velocity
  
  for i = 2:imax-1
    
    for j = 2:jmax
      
      AA(i,j) = nu*( (UU(i+1,j) - 2*UU(i,j) + UU(i-1,j))/dx^2 + ...
      (UU(i,j+1) - 2*UU(i,j) + UU(i,j-1))/dy^2 );
      
      %    BB(i,j) = (1/dx)*( ( (UU(i,j)+UU(i+1,j))/2 )^2 - ...
      %    ( (UU(i-1,j)+UU(i,j))/2 )^2 );
      
      BB(i,j) = (1/dx)*( ( (UU(i,j)+UU(i+1,j))/2 )^2 - ...
      ( (UU(i-1,j)+UU(i,j))/2 )^2 ) + (gamma/dx) * ...
      ( (0.25*abs(UU(i,j)+UU(i+1,j))*(UU(i,j)-UU(i+1,j))) - ...
      (0.25*abs(UU(i-1,j)+UU(i,j))*(UU(i-1,j)-UU(i,j))) );
      
      %    DD(i,j) = (1/dy) * ( (0.25*(VV(i,j)+VV(i+1,j))*(UU(i,j+1)+UU(i,j))) - ...
      %      (0.25*(VV(i,j-1)+VV(i+1,j-1))*(UU(i,j-1)+UU(i,j))) );
      
      DD(i,j) = (1/dy) * ( (0.25*(VV(i,j)+VV(i+1,j))*(UU(i,j+1)+UU(i,j))) - ...
      (0.25*(VV(i,j-1)+VV(i+1,j-1))*(UU(i,j-1)+UU(i,j))) ) + (gamma/dy) * ...
      ( (0.25*abs(VV(i,j)+VV(i+1,j))*(UU(i,j)-UU(i,j+1))) - ...
      (0.25*abs(VV(i,j-1)+VV(i+1,j-1))*(UU(i,j-1)-UU(i,j))) );
      
      F(i,j) = UU(i,j) + dt*( AA(i,j) - BB(i,j) - DD(i,j));
      
    end
    
  end
  
  %Calculate intermediate v velocity
  for i = 2:imax
    
    for j = 2:jmax-1
      
      FF(i,j) = nu* ( (VV(i+1,j) - 2*VV(i,j) + VV(i-1,j))/dx^2 + ...
      (VV(i,j+1) - 2*VV(i,j) + VV(i,j-1))/dy^2 );
      
      %    GG(i,j) = (1/dy)*( ( (VV(i,j)+VV(i,j+1))/2 )^2 - ...
      %    ( (VV(i,j-1)+VV(i,j))/2 )^2 );
      
      GG(i,j) = (1/dy)*( ( (VV(i,j)+VV(i,j+1))/2 )^2 - ...
      ( (VV(i,j-1)+VV(i,j))/2 )^2 ) + (gamma/dy) * ...
      ( (0.25*abs(VV(i,j)+VV(i,j+1))*(VV(i,j)-VV(i,j+1))) - ...
      (0.25*abs(VV(i,j-1)+VV(i,j))*(VV(i,j-1)-VV(i,j))) );
      
      %    II(i,j) = (1/dx) * ( (0.25*(VV(i+1,j)+VV(i,j))*(UU(i,j+1)+UU(i,j))) - ...
      %      (0.25*(VV(i-1,j)+VV(i,j))*(UU(i-1,j+1)+UU(i-1,j))) );
      
      II(i,j) = (1/dx) * ( (0.25*(VV(i+1,j)+VV(i,j))*(UU(i,j+1)+UU(i,j))) - ...
      (0.25*(VV(i-1,j)+VV(i,j))*(UU(i-1,j+1)+UU(i-1,j))) ) + (gamma/dx) * ...
      ( (0.25*abs(UU(i,j+1)+UU(i,j))*(VV(i,j)-VV(i+1,j))) - ...
      (0.25*abs(UU(i-1,j+1)+UU(i-1,j))*(VV(i-1,j)-VV(i,j))) );
      
      G(i,j) = VV(i,j) + dt*( FF(i,j) - II(i,j) - GG(i,j));
      
    end
    
  end
  
  %Compute PPE Source Term
  
  for i = 2:imax
    
    for j = 2:jmax
      
      g(i,j) = (rho/dt)*( (F(i,j) - F(i-1,j))/dx + (G(i,j) - G(i,j-1))/dy );
      
    end
    
  end
  
  p = zeros(imax+1,jmax+1);
  p2 = zeros(imax+1,jmax+1);
  Res = zeros(imax+1,jmax+1);
  RES = 0.5;
  iteration = 1;
  Tol = max(max(abs(g)))*Tolfac;
  
  %%%%%%%% PPE - SOR%%%%%%%%
%  Tol
%  RES
  while (RES > Tol) && (iteration < 3500)
%     fprintf('%f %d %f\n', n, iteration, RES);
    
    for i = 2:imax
      
      for j = 2:jmax
        
        p(i,j) = p2(i,j)*(1-w) + ( w/( (i>2) + (j<jmax) + (i<imax) + 1 ) )*...
        ( (i<imax)*p2(i+1,j) + (i>2)*p(i-1,j) + (j<jmax)*p2(i,j+1) + ...
        p(i,j-1) - g(i,j)*dx*dx );
        
        %      Res(i,j) = ( (i<imax)*(p(i+1,j)-p(i,j)) + (i>2)*(p(i-1,j)-p(i,j)) + ...
        %                 (j<jmax)*(p(i,j+1)-p(i,j)) + (p(i,j-1) - p(i,j)) )/(dx^2) - ...
        %                 g(i,j);
        
      end
      
    end
    
    for i = 2:imax
      
      for j = 2:jmax
        
        Res(i,j) = ( (i<imax)*(p(i+1,j)-p(i,j)) + (i>2)*(p(i-1,j)-p(i,j)) + ...
        (j<jmax)*(p(i,j+1)-p(i,j)) + (p(i,j-1) - p(i,j)) )/(dx^2) - ...
        g(i,j);
        
      end
      
    end
    
    
    p2 = p;
    RES = max(max(abs(Res)));
    iteration = iteration +1;       
    
  end
%   fprintf('%d %f\n', iteration, RES)
  %Correct u velocity
  for i = 2:imax-1
    
    for j = 2:jmax
      
      UU(i,j) = F(i,j) - (dt/dx)*(p(i+1,j) - p(i,j));
      
    end
    
  end
  %Correct v velocity
  for i = 2:imax
    
    for j = 2:jmax-1
      
      VV(i,j) = G(i,j) - (dt/dx)*( p(i,j+1) - p(i,j) );
      
    end
    
  end
  
  ke(nn) = sum(sum(0.5*rho*UU(1:imax,2:jmax).^2))+sum(sum(0.5*rho*VV(2:imax,i:jmax).^2));
  
 for i = 2:imax
   
   for j = 2:jmax
     
     UU_P(i,j) = ( UU(i,j) + UU(i-1,j))/2;
     VV_P(i,j) = ( VV(i,j) + VV(i,j-1))/2;
     Velocity_mag(i,j) = sqrt(UU_P(i,j)^2+VV_P(i,j)^2);
     
   end
   
 end
 
 UU_P = UU_P./Velocity_mag;
 VV_P = VV_P./Velocity_mag;
 
%  figure(1)
%  contourf(XP,YP,Velocity_mag)
%  hold on
%  quiver(XP,YP,UU_P,VV_P)
%  hold off
%  grid on
%  xlim([0 1])
%  ylim([0 1])
%  colorbar
%  title('Velocity Magnitude')
%  drawnow
%  
%  figure(2)
%  contourf(XP,YP,p)
%  title('pressure')
%  grid on
%  xlim([0 1])
%  ylim([0 1])
%  colorbar
%  drawnow
  idx = find(XU(:,1)==0.5);
  idy = find(YV(1,:)==0.5);

  figure(3)
  plot(y_p,u_re100,'rd',YU(1,:),UU(idx,:),'m-')
  xlim([0 1])
  title('U Velocity Profile')
  xlabel('Y Distance')
  ylabel('U Velocity')
  legend('Validation Data','Simulation')

  grid on
  drawnow

  figure(4)
  plot(x_p,v_re100,'rd',XV(:,1),VV(:,idy),'m-')
  xlim([0 1])
  title('V Velocity Profile')
  grid on
  xlabel('X distance')
  ylabel('V Velocity')
  legend('Validation Data','Simulation')
  drawnow
  

  
  nn = nn+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%% End Time Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(5)
  plot(dt:dt:t_end,ke)
  title('Kinetic Energy')
  grid on
  xlabel('Time')
  ylabel('Kinetic Energy')
  drawnow
  for i = 2:imax
    
    for j = 2:jmax
      
      UU_P(i,j) = ( UU(i,j) + UU(i-1,j))/2;
      VV_P(i,j) = ( VV(i,j) + VV(i,j-1))/2;
      Velocity_mag(i,j) = sqrt(UU_P(i,j)^2+VV_P(i,j)^2);
      
    end
    
  end
  
  UU_P = UU_P./Velocity_mag;
  VV_P = VV_P./Velocity_mag;
  
%   figure
%   contourf(XP,YP,Velocity_mag)
%   hold on
%   quiver(XP,YP,UU_P,VV_P)
%   hold off
%   grid on
%   xlim([0 1])
%   ylim([0 1])
%   colorbar
%   title('Velocity Magnitude')
%  drawnow
  figure
%   figure('units','normalized','position',[0 0.33 .3 .3])
  contourf(XP,YP,p)
  title('Pressure')
  grid on
  xlim([0 1])
  ylim([0 1])
  colorbar
%  drawnow
figure
% figure('units','normalized','position',[0 0.33 .3 .3])
contourf(XU,YU,UU)
grid on
title('U velocity')
xlim([0 1])
ylim([0 1])
colorbar
figure
% figure('units','normalized','position',[0 0.33 .3 .3])
contourf(XV,YV,VV)
title('V velocity')
grid on
xlim([0 1])
ylim([0 1])
colorbar

% figure
% contourf(XP,YP,g)
% title('source')
% grid on
% xlim([0 1])
% ylim([0 1])
% colorbar

% figure
% contourf(XP,YP,Res)
% title('Residual')
% grid on
% xlim([0 1])
% ylim([0 1])
% colorbar