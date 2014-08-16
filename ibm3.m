% ibm on staggered grid with frational step for 2D flow prolbem with 
% predefined movement equation of solid body
% auther zhengjiang li
% ---------------IBM algorithms describtion--------------------------------
% step1:: descritize the whole domain with Cartesian grid, descritize the
% solid body with only boundary marker points 
%
% in each time step do the following :
% claculate estimated velocity field in the whole domain by equation
% Du/Dt = -grad(P) + miu*laplace(u)
% this velocity can't satisfy the interface boundary
%
% step2:: identify the grid points into three class: fluid points, solid 
% points,forcing points
% fluid points: whose 4 neighbor points should all in the fluid domain
% solid points: inside the interface boundary
% forcing points: points in fluid domain with at least one neighbor point
% in the solid domain
%
% geometric detection algorithm( to identify the three points)
% for each marker points, search for the closest grid point, and mark the 9
% neighbors as vicinity points.
%
% for each marker point, we know its normal direction vector based on
% the "arc length" parameter-representation x, y, and these normal vector
% "vec_n" is pointing outside, so basically sign of the dot product of
% vec_n and (xi -int_x) or sign of cross product of vec_tan and (xi - int_x)
% can show whether this point(xi) is inside the interface or outside.
% 
% step2.2, identify forcing points
% by definiation, frist, forcing points should be in the fluid domain
% then test if at least one of the four diagonal points is in solid domain
%
% for predefined movment equation of the body, this part should have
% simpler implementation
%
% step3 :: velocity linear interpolation on forcing points(also higher 
% accuate interpolation is good but need more points information)
% we use two fluid points and one interface point to construct the linear
% dual-variable equation :  phi = b0 + x*b1 + y*b2
% so we get the coefficient and then get the velocity value "v@fp" at 
% forcing point
%
%step4:: calcualting forcing term at forcing points
% f@tn+1 = v@fp - v@tn - RHS
% RHS= miu*laplace(v) - udiv(u) - gradp/rho, all at tn
% there are many fractional ways here in the presentation papers
%
%step5:: recalculate the momentum equation with forcing term and come to 
% projection method to calculate pressure at tn+1, and update v@tn+1
% reference:
% 1 E. Balaras Computers & Fluid 33 (2004) 375-404
% 2 J. Yang Balaras Journal of Computational Physics 215 (2006) 12-40
% 3 Haoxiang Luo, Computers & Fluids 56 (2012) 61-76
% 4 Fotis Sotiropoulos, Progress in Aerospace Sciences 65 (2014) 1-21
% 5 Y. H. Teng. Journal of Computational physics 192 (2003) 593-623
% 6 C.C. Liao, Computers & Fluids 39 (2010) 152-167
% 7 J. Yang, Journal of Computational Physics 231 (2012) 5029-5061
%----------------end of algorithm description---------------------------
%--------------interface mesh with N points-----------------------------
N = 10;  % interface pints #
s = linspace(0, 2*pi, N+1); 
s(N+1) = [];
ds = 2*pi/N; % arclength, evenly distributed
R = 0.2; % radius of the circle in the flow
xx = 0.5 + R*cos(s);  % initial coord of interface points
yy = 0.5 + R*sin(s);  
% later in time iteration, as center of circle move,
% y0 need add v0*dt, and x0 need add u0*dt
% 
% normal vector of interface points
list_min = zeros(1,N); 
list_sub = zeros(1,N);
list_min(1)=N;
list_sub(N)=1;
for i=2:1:N
	list_min(i)=i-1;
end
for j=1:1:(N-1)
	list_sub(j)=j+1;
end
dxds(1:N)=(xx(list_sub) - xx(list_min))/ds/2;
dyds(1:N)=(yy(list_sub) - yy(list_min))/ds/2;
normalx = dyds; % x component of normal vector
normaly = -dxds; % y component of normal vector
%----------------prescribed body move equation-------------------------
% yy = yy + R*sin(2*pi*t); xx = xx; 
%-------------------------background mesh -----------------------------

Lx = 0;
Rx = 1;
By = 0;
Ty = 1; 
hx = 0.04; % approx with marker point spacing
hy = 0.04; 
nx = 1/hx;
ny = 1/hy; 
x =linspace(0,1,nx+1);
y =linspace(0,1,ny+1); 
[X,Y]=meshgrid(x,y);

% figure (1)
%  plot(X,Y), hold on
%  plot(Y,X), hold on
% plot([xx,xx(1)], [yy,yy(1)], '-')


incellx = zeros(1,N);
incelly = zeros(1,N);
valuep = zeros(nx, ny); % notation of different kind of points
mp2gpx = zeros(N,5);
mp2gpy = zeros(N,5);
% exten u v ( physical u, v plus dummy boundary u,v) 
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
tf = 100;
dt = 0.001;
Re = 10; 

for t = 1: tf

%------------------------------------------------------------------------
% step one
% claculate estimated velocity field in the whole domain by equation
% Du/Dt = -grad(P) + miu*laplace(u)   
% fixed boundary condition     
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
   rhsu = zeros(nx-1,ny);
  for i = 2: nx
        for j = 2: ny+1    
            ii = i-1;
            jj = j-1;
 visu(ii,jj)=((Ue(i-1,j) - 2*Ue(i,j) + Ue(i+1,j))/hx^2 + ...
            (Ue(i,j+1)- 2*Ue(i,j) + Ue(i,j-1))/hy^2)/Re ;    
 
 du2dx(ii,jj) = (Ue(i+1,j)^2 - Ue(i-1,j)^2 + 2*Ue(i,j)*(Ue(i+1,j) - ...
             Ue(i-1,j)))/hx/4;
         
 duvdy(ii,jj) = (Ue(i,j) + Ue(i,j+1))*(Ve(i,j) + Ve(i+1,j))/4/hy -  ...
              (Ue(i,j) + Ue(i,j-1))*(Ve(i,j-1) + Ve(i+1,j-1))/4/hy;
       
 rhsu(ii,jj) = visu(ii,jj) - duvdy(ii,jj) - du2dx(ii,jj); 
 Ue(i,j) = Ue(i,j) + rhsu(ii,jj)*dt;
        end
  end
    
visv = zeros(nx,ny-1);
dv2dy = zeros(nx,ny-1);
duvdx = zeros(nx,ny-1);
rhsv = zeros(nx, ny-1);
    for i = 2: nx+1
        for j=2:ny
            ii = i-1;
            jj = j-1;          
 duvdx(ii,jj) = (Ue(i,j)+Ue(i,j+1))*(Ve(i,j)+Ve(i+1,j))/hx/4 - ...
              (Ue(i-1,j)+Ue(i-1,j+1))*(Ve(i-1,j)+Ve(i,j))/hx/4;

 dv2dy(ii,jj) = (Ve(i,j+1)^2 - Ve(i,j-1)^2 + 2*Ve(i,j)*(Ve(i,j+1) - ...
              Ve(i,j-1)))/hy/4;
          
 visv(ii,jj) = ((Ve(i-1,j) - 2*Ve(i,j) + Ve(i+1,j))/hx^2 + ...
             (Ve(i,j+1) - 2*Ve(i,j) + Ve(i,j-1))/hy^2)/Re;
         
 rhsv(ii,jj) = visv(ii,jj) - duvdx(ii,jj) - dv2dy(ii,jj);
 Ve(i,j) = Ve(i,j) + rhsv(ii,jj)*dt;
        end
    end
%    
%-----------obtain estimated velocity at all grid points above---------

%---------------------- fluid, solid, forcing points determine------------
%  step two
for i = 1: N
    incellx(i)=(ceil((xx(i)-Lx)/hx)-1); %interface point in which cell
    incelly(i)=(ceil((yy(i)-By)/hy)-1);
% % mark vicinity background points with valuep = 100;  
% valuep(incellx(i):incellx(i)+1,incelly(i):incelly(i)+1) = 100;
% set a maker point - grid point matrix
mp2gpx(i,:)= [i,incellx(i),incellx(i)+1,incellx(i)+1,incellx(i)];
mp2gpy(i,:)= [i,incelly(i),incelly(i),incelly(i)+1,incelly(i)+1];
end
    boxxmax = max(max(mp2gpx(:,2:5)));
    boxxmin = min(min(mp2gpx(:,2:5)));
    boxymax = max(max(mp2gpy(:,2:5)));
    boxymin = min(min(mp2gpy(:,2:5)));
    
    
    for i = 1: nx+1
        for j=1:ny+1
           if X(i,j) < boxxmin|| X(i,j) > boxxmax || ...
              Y(i,j) < boxymin|| Y(i,j) > boxymax
              valuep(i,j) = 1; % fluid points outside of the bounding box
%            else
%               valuep(i,j) = -1; % solid and all forcing points
           end   
        end
    end
    % all forcing points mark as valuep=5
    valuep(boxxmin:boxxmax,boxymin) = 5;
    valuep(boxxmin:boxxmax,boxymax)=5;
    valuep(boxxmin,boxymin:boxymax)=5;
    valuep(boxxmax,boxymin:boxymax)=5;
    
    % interpolat points in fluid domain, mark as valuep=20
    valuep(boxxmin-1,boxymin-1:boxymax+1)=20;
    valuep(boxxmax+1,boxymin-1:boxymax+1)=20;
    valuep(boxxmin-1:boxxmax+1,boxymin-1)=20;
    valuep(boxxmin-1:boxxmax+1,boxymax+1)=20;
    
    valuep(boxxmin+1:boxxmax-1, boxymin+1:boxymax-1) = -1; %  solid points
    
    % ?? detect a point inside a triangular??
    
    % care about index, PHYSICAL velocity  Ue(:,2:end-1), Ve(2:end-1,:)
    % size(Ue(:,2:end-1)) = (26,25);
    % size(Ve(2:end-1,:)) = (25,26);
    
    % physical u, v
    
    phyu = rot90(Ue(:,2:end-1));
    phyv = rot90(Ve(2:end-1,:));
    bp2fp = zeros(1,N);
    force.x = zeros(nx-1, ny);
    force.y = zeros(nx, ny-1);
    for i=boxxmin:boxxmax
        for j=boxymin:boxymax
            if valuep(i,j) == 5
                for ii = 1:N
                    bp2fp(ii) = ((x(i) - xx(ii))^2 + (y(j) - yy(ii))^2)^0.5;
                end
                [mindis, mindp] = min(bp2fp);
                interp1.x = xx(mindp);
                interp1.y = yy(mindp);
                interp1.u = 0;
                interp1.v = 0;               
                
                   k=(interp1.y - y(j))/(interp1.x-x(i));
                   determ1 = (interp1.y - y(j))/k + x(i);
                   determ2 = (interp1.x - x(i))*k + y(j);  
                 
                   
                if i == boxxmin  % bottom
                   interp2.x = x(i);
                   interp2.y = y(j-1);
                   interp2.u = phyu(i,j-1);
                   interp2.v = phyv(i,j-1);
                   interp3.y = interp2.y;
                   if determ1 > x(i)
                      if determ1 > x(i+1)
                          interp3.x = x(i+2);
                          interp3.u = phyu(i+2,j-1);
                          interp3.v = phyv(i+2,j-1);
                      else
                       interp3.x = x(i+1); 
                       interp3.u = phyu(i+1,j-1);
                       interp3.v = phyv(i+1,j-1);
                      end
                   else
                       if determ1 < x(i-1)
                           interp3.x = x(i-2);
                           interp3.u = phyu(i-2,j-1);
                           interp3.v = phyv(i-2,j-1);
                       else
                       interp3.x = x(i-1);
                       interp3.u = phyu(i-1,j-1);
                       interp3.v = phyv(i-1,j-1);
                       end
                   end
                end
                
                if i == boxxmax % top
                    interp2.x =x(i);
                    interp2.y = y(j+1);
                    interp2.u = phyu(i,j+1);
                    interp2.v = phyv(i,j+1);
                    interp3.y = interp2.y;
                    if determ1 > x(i)
                       if determ1 > x(i+1)
                          interp3.x = x(i+2);
                          interp3.u = phyu(i+2,j+1);
                          interp3.v = phyv(i+2,j+1);
                       else
                       interp3.x = x(i+1);
                       interp3.u = phyu(i+1,j+1);
                       interp3.v = phyv(i+1,j+1);
                       end
                    else
                        if determ1 < x(i-1)
                            interp3.x = x(i-2);
                            interp3.u = phyu(i-2,j+1);
                            interp3.v = phyv(i-2,j+1);
                        else
                            interp3.x = x(i-1);
                            interp3.u = phyu(i-1,j+1);
                            interp3.v = phyv(i-1,j+1);
                        end
                    end
                end
                
                if j == boxymin % left
                   interp2.x = x(i-1);
                   interp2.y = y(j);
                   interp2.u = phyu(i-1,j);
                   interp2.v = phyv(i-1,j);
                   interp3.x = interp2.x;
                   if determ2 > y(j)
                      if determ2 > y(j+1)
                          interp3.y = y(j+2);
                          interp3.u = phyu(i-1,j+2);
                          interp3.v = phyv(i-1,j+2);
                      else
                      interp3.y = y(j+1);
                      interp3.u = phyu(i-1,j+1);
                      interp3.v = phyv(i-1,j+1);
                      end
                   else
                       if determ2 < y(j-1)
                           interp3.y = y(j-2);
                           interp3.u = phyu(i-1,j-2);
                           interp3.v = phyv(i-1,j-2);
                       else
                           interp3.y = y(j-1);
                           interp3.u = phyu(i-1,j-1);
                           interp3.v = phyv(i-1,j-1);
                       end
                   end
                end
                
                if j == boxymax % right
                   interp2.x = x(i+1);
                   interp2.y = y(j);
                   interp2.u = phyu(i+1,j);
                   interp2.v = phyv(i+1,j);
                   interp3.x = interp2.x;
                   if determ2 > y(j) 
                       if determ2 > y(j+1)
                             interp3.y = y(j+2); 
                             interp3.u = phyu(i+1,j+2);
                             interp3.v = phyv(i+1,j+2);
                       else
                           interp3.y = y(j+1);
                           interp3.u = phyu(i+1,j+1);
                           interp3.v = phyv(i+1,j+1);                          
                       end
                   else
                       if determ2 < y(j-1) 
                           interp3.y = y(j-2);
                           interp3.u = phyu(i+1,j-2);
                           interp3.v = phyv(i+1,j-2);
                       else
                           interp3.y = y(j-1);
                           interp3.u = phyu(i+1,j-1);
                           interp3.v = phyv(i+1,j-1);
                       end
                   end
                end 
                
                                         
% here we get all three interpolation points      
   matphi = [1 interp1.x interp1.y; 1 interp2.x interp2.y; ...
             1 interp3.x interp3.y];
   vetu = [interp1.u ; interp2.u; interp3.u];
   vetv = [interp1.v; interp2.v; interp3.v];
   coeffu = matphi\vetu;
   coeffv = matphi\vetv;
   % get the interpolation velocity at forcing point
   fp.u(i,j) = coeffu(1) + x(i)*coeffu(2) + y(j)*coeffu(3);
   fp.v(i,j) = coeffv(1) + x(i)*coeffv(2) + y(j)*coeffv(3);
%-------------------------------------------------------------------------
% step 4
%  calcualting forcing term at forcing points
% fu@tn+1 = (fp.u - u)@fp/dt - rhsu
% fv@tn+1 (fp.y -u)@fp/dt - rhv
% update u component  
% one prolbem, fpx, fpy here are on grid points, not the half points used
% in the staggered grids. 
 force.x(i,j) = (fp.u(i,j) - Ue(i+1,j+1))/dt - rhsu(i,j);         
%update v components
 force.y(i,j) = (fp.v(i,j) - Ve(i+1,j+1))/dt - rhsv(i,j);
 
        end
    end
end
% ------------------------------------------------------------------------
% step 5
% recalculate the momentum equation with forcing term and come to 
% projection method to calculate pressure at tn+1, and update v@tn+1
% 
%  this part use the same schedule as projection method on staggered grid
% only one difference, in RHS there is a forcing term
% and for stability, at first calculate the estimate velocity and here
% reculate the velocity again. be carefule to choose a good fractional
% factor actually. 
 
        for i = 2: nx
        for j = 2: ny+1   
            ii = i -1;
            jj = j -1; 
 visu(ii,jj)=((Ue(i-1,j) - 2*Ue(i,j) + Ue(i+1,j))/hx^2 + ...
            (Ue(i,j+1)- 2*Ue(i,j) + Ue(i,j-1))/hy^2)/Re ;    
 
 du2dx(ii,jj) = (Ue(i+1,j)^2 - Ue(i-1,j)^2 + 2*Ue(i,j)*(Ue(i+1,j) - ...
             Ue(i-1,j)))/hx/4;
         
 duvdy(ii,jj) = (Ue(i,j) + Ue(i,j+1))*(Ve(i,j) + Ve(i+1,j))/4/hy -  ...
              (Ue(i,j) + Ue(i,j-1))*(Ve(i,j-1) + Ve(i+1,j-1))/4/hy;
       
 rhsu(ii,jj) = visu(ii,jj) - duvdy(ii,jj) - du2dx(ii,jj); 
 Ue(i,j) = Ue(i,j) + (rhsu(ii,jj)+force.x(ii,jj))*dt;
        end
        end
    
        for i = 2: nx+1
        for j=2:ny
            ii = i-1;
            jj = j-1;
 duvdx(ii,jj) = (Ue(i,j)+Ue(i,j+1))*(Ve(i,j)+Ve(i+1,j))/hx/4 - ...
              (Ue(i-1,j)+Ue(i-1,j+1))*(Ve(i-1,j)+Ve(i,j))/hx/4;

 dv2dy(ii,jj) = (Ve(i,j+1)^2 - Ve(i,j-1)^2 + 2*Ve(i,j)*(Ve(i,j+1) - ...
              Ve(i,j-1)))/hy/4;
          
 visv(ii,jj) = ((Ve(i-1,j) - 2*Ve(i,j) + Ve(i+1,j))/hx^2 + ...
             (Ve(i,j+1) - 2*Ve(i,j) + Ve(i,j-1))/hy^2)/Re;
         
 rhsv(ii,jj) = visv(ii,jj) - duvdx(ii,jj) - dv2dy(ii,jj);
 Ve(i,j) = Ve(i,j) + (rhsv(ii,jj)+force.y(ii,jj))*dt;
        end
        end
        
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
               
end
            
 figure (1)
 plot(X,Y), hold on
 plot(Y,X), hold on
 plot([xx,xx(1)], [yy,yy(1)], '-'),hold on
 
 physu = Ue(:,2:end-1);
 physv = Ve(2:end-1,:);
 fieldu = [zeros(nx+1,1) avg(physu')' ones(nx+1,1)];
 fieldv = [zeros(1,nx+1);  avg(physv); zeros(1,nx+1)];
%Plot velocity u
figure(3)
surf(x, y, fieldu)
view(-90,90)
xlabel('x-axis'), ylabel('y-axis')
%Plot velocity v
figure(4);
surf(x, y, fieldv)
view(-90,90)
xlabel('x-axis'), ylabel('y-axis')












