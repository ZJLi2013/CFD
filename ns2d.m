% cfd midterm project

nx = 40;  % number of cells in x direction
ny = 40;  % number of cells in y direction

x =linspace(0,1,nx+1);
y =linspace(0,1,ny+1); 
hx = x(2) - x(1);
hy = y(2) - y(1);
[X,Y]=meshgrid(x,y);
% exten u v 
Ue = zeros(nx+1,ny+2);  
Ve = zeros(nx+2,ny+1);
% interior u v p
% U = zeros(nx-1,ny);
% V = zeros(nx,ny-1);
% P = zeros(nx,ny);
% % boundary contions
% vN = 0; uN = 1; % in physical domain
% vS = 0; uS = 0;
% vW = 0; uW = 0;
% vE = 0; uE = 0;

%-------------------------------------

tf = 1000;
dt = 0.001;
Re = 10; 
for t = 1: tf
%    % physical North
% Ve(:,end) = vN; % size(vN) = nx+2
% Ue(:,end) = 2*uN - Ue(:,end-1); % size(uN)=nx+2;
% % South
% Ve(:,1) = vS;
% Ue(:,1) = 2*uS - Ue(:,2);
% % West
% Ve(1,:)= 2*vW - Ve(2,:);
% Ue(1,:)= uW;
% % East
% Ve(end,:) = 2*vE - Ve(end-1,:);
% Ue(end,:) = uE; 

% modidfy on Saturday
% refer, bc for 2d lid driven cavity using ghost cells
% Ue(:,end) = 2 - Ue(:,end-1);
% Ve(:,end) = -Ve(:,end-1);
% 
% Ue(:,1) = - Ue(:,2);
% Ve(:,1) = - Ve(:,2); 
% 
% Ue(1,:) = -Ue(2,:);
% Ve(1,:) = -Ve(2,:);
% 
% Ue(end,:) = -Ue(end-1,:);
% Ve(end,:) = -Ve(end-1,:);

% north
Ue(:,end) = 2 - Ue(:,end-1);
Ve(:,end) = 0;
% south
Ue(:,1) = -Ue(:,2);
Ve(:,1) = 0; 
%west
Ue(1,:) = 0;
Ve(1,:) = -Ve(2,:);
%east
Ue(end,:) = 0;
Ve(end,:) = -Ve(end-1,:);
    
   visu = zeros(nx-1, ny);
   du2dx = zeros(nx-1,ny);
   duvdy = zeros(nx-1,ny); 
  for i = 2: nx
        for j = 2: ny+1       
 visu(i,j)=((Ue(i-1,j) - 2*Ue(i,j) + Ue(i+1,j))/hx^2 + ...
            (Ue(i,j+1)- 2*Ue(i,j) + Ue(i,j-1))/hy^2)/Re ;    
 
 du2dx(i,j) = (Ue(i+1,j)^2 - Ue(i-1,j)^2 + 2*Ue(i,j)*(Ue(i+1,j) - ...
             Ue(i-1,j)))/hx/4;
         
 duvdy(i,j) = (Ue(i,j) + Ue(i,j+1))*(Ve(i,j) + Ve(i+1,j))/4/hy -  ...
              (Ue(i,j) + Ue(i,j-1))*(Ve(i,j-1) + Ve(i+1,j-1))/4/hy;
       
 Ue(i,j) = Ue(i,j) + visu(i,j)*dt - (duvdy(i,j) + du2dx(i,j))*dt;
        end
  end
%  U = Ue(2:end-1,2:end-1);
    
visv = zeros(nx,ny-1);
dv2dy = zeros(nx,ny-1);
duvdx = zeros(nx,ny-1);
    for i = 2: nx+1
        for j=2:ny
 duvdx(i,j) = (Ue(i,j)+Ue(i,j+1))*(Ve(i,j)+Ve(i+1,j))/hx/4 - ...
              (Ue(i-1,j)+Ue(i-1,j+1))*(Ve(i-1,j)+Ve(i,j))/hx/4;

 dv2dy(i,j) = (Ve(i,j+1)^2 - Ve(i,j-1)^2 + 2*Ve(i,j)*(Ve(i,j+1) - ...
              Ve(i,j-1)))/hy/4;
          
 visv(i,j) = ((Ve(i-1,j) - 2*Ve(i,j) + Ve(i+1,j))/hx^2 + ...
             (Ve(i,j+1) - 2*Ve(i,j) + Ve(i,j-1))/hy^2)/Re;
         
 Ve(i,j) = Ve(i,j) + visv(i,j)*dt - (duvdx(i,j) + dv2dy(i,j))*dt;
        end
    end
%  V = Ve(2:end-1,2:end-1);
    
    dudx =  diff(Ue(:,2:end-1))/hx;
    dvdy =  diff(Ve(2:end-1,:),1,2)/hy;  
    rhs = reshape((dudx + dvdy)/dt, nx*ny, 1);
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx,1:nx-1,1,nx,nx);
Ax=Ex+Ex'-2*speye(nx);        %Dirichlet B.Cs
% axx = full(Ax);
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny,1:ny-1,1,ny,ny);
Ay=Ey+Ey'-2*speye(ny);        %Dirichlet B.Cs
Ay(1,1)=-1; Ay(ny,ny)=-1;  %Neumann B.Cs
A=kron(Ay/hy^2,speye(nx))+kron(speye(ny),Ax/hx^2);
P = A \ rhs;
P= reshape(P, nx, ny); 

   dpdx = diff(P,1,1);
   dpdy = diff(P,1,2);
   Ue(2:end-1,2:end-1) = Ue(2:end-1,2:end-1) - dt*dpdx;
   Ve(2:end-1,2:end-1) = Ve(2:end-1,2:end-1) - dt*dpdy;
   
   
% % visulization 
%  physu = Ue(:,2:end-1);
%  physv = Ve(2:end-1,:);
%  fieldu = [zeros(nx+1,1) avg(physu')' ones(nx+1,1)];
%  fieldv = [zeros(1,nx+1);  avg(physv); zeros(1,nx+1)];
%   quiver(X,Y,fieldu',fieldv');
%   axis([0 1 0 1]);

end
 


 physu = Ue(:,2:end-1);
 physv = Ve(2:end-1,:);
 fieldu = [zeros(nx+1,1) avg(physu')' ones(nx+1,1)];
 fieldv = [zeros(1,nx+1);  avg(physv); zeros(1,nx+1)];
 
%   % vector plot
%   figure(1)
%   quiver(X,Y,fieldu',fieldv');
%   axis([0 1 0 1]);
  
%plot velocity components at x=0.5
 figure(2)
 subplot(2,1,1)
plot(x,fieldu(:,nx/2),'r','Marker','.')
hold on
plot(x,fieldv(:,nx/2),'b','Marker','.')
hold on
title('Velocity Components - Vertical Dircetion');
xlabel('Vertical length y');
ylabel('u,v @x=0.5');
%plot velocity components at y=0.5
subplot(2,1,2)
plot(y,fieldu(ny/2,:),'r','Marker','.')
hold on
plot(y,fieldv(ny/2,:),'b','Marker','.')
hold on
title('Velocity Components - Horizontal Dircetion');
xlabel('Vertical length x');
ylabel('u,v @y=0.5');

% %
% Pe = zeros(nx+2, ny+2);
% Pe(2:end-1,2:end-1)= P;
% Pe(1,2:end-1) = P(1,:) ; % dp/dn = 0 @ west
% Pe(end,2:end-1) = P(end,:); % @east
% Pe(2:end-1,1) = P(:,1); % @south
% Pe(2:end-1,end)=P(:,end); % north
% fieldp = avg(avg(Pe)');
% 
% 
% % Plot pressure at y=0.5
% figure(2);
% subplot(2,1,1)
% plot(y,fieldp(nx/2,:),'r','Marker','.')
% hold on
% title('Pressure Along Horizontal Direction');
% xlabel('Horizontal Length x');
% ylabel('Pressure at y = 0.5');
% 
% % Plot pressure at x=0.5
% %pff=fliplr(pf);
% subplot(2,1,2)
% plot(x,fieldp(:,nx/2),'r','Marker','.')
% hold on
% view(0,-90)
% title('Pressure Along Vertical Direction');
% xlabel('Vertical Length y');
% ylabel('Pressure at x = 0.5');
% 
% 
% %Plot velocity u
% figure(3)
% surf(x, y, fieldu)
% view(-90,90)
% xlabel('x-axis'), ylabel('y-axis')
% %Plot velocity v
% figure(4);
% surf(x, y, fieldv)
% view(-90,90)
% xlabel('x-axis'), ylabel('y-axis')

   




 