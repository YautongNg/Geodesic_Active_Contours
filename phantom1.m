% Geodesic Active Contours


clearvars
close all
clc

[I0, ~] = LoadData(1);
[y,x] = size(I0); % Get the dimensions important for differencing
I = I0(:); % Reshape the image into a vector: L-R columns stack T-B

% Pick some parameters here
maxIter = 600; % Number of times to repeat gf + reinit steps
tol = 1e-6;
a = 5;
exitflag = 0; % Iteration flag, 0 for incomplete, 1 for complete
dt = 0.01;
gfIter = 10; % Number of iterations to take before reinitialization
reinitIter = 10; % Number of iterations to take in the reinitialization phase
epsin = 1e-6;
% edge detector g calculated below


%% Construct Operators
% y = 3; x = 4;
n = y*x;
e = ones(n,1);
e1 = ones(n,1);
for i = y:y:n
    e1(i) = 0;
end
e2 = ones(n,1);
for i = y+1:y:n
    e2(i) = 0;
end
e3 = ones(n,1);
for i = y:y:n
    e3(i) = 0;
end

% Lapalcian
L = spdiags([e e1 -4*e e2 e],[-y -1 0 1 y],n,n);

% x-Gradient
Dpx = spdiags([-e e],[0 y],n,n);
Dmx = spdiags([-e e],[-y 1],n,n);
Dcx = (1/2)*spdiags([-e e],[-y y],n,n);

% y-Gradient
Dpy = spdiags([e1 -e],[-1 0],n,n);
Dmy = spdiags([e -e1],[0 1],n,n);
Dcy = (1/2)*spdiags([e3 -e1],[-1 1],n,n);


%% Edge Detector
% Using g = 1./(1+||grad(I)||^2), here I is the initial image
% g should be a vector same length as u

g = 1./(1 + ((Dcx*I).^2 + (Dcy*I).^2));

g1 = max(Dcx*g,0);
g2 = min(Dcx*g,0);
g3 = max(Dcy*g,0);
g4 = min(Dcy*g,0);


%% Construct function u
% Recall that C is the zero level set of some function u
[X,Y] = meshgrid(1:x,1:y);
u = sqrt((X-1/2*x).^2 + (Y-1/2*y).^2) - 200;
% u = sin(Y./1.5) + cos(X./1.5) - 0.5;
u0 = u(:);

surf(u), shading interp


%% Main iterations
u_nxt = u0;
figure(2), clf
for k = 1:maxIter
    disp(['k is ' num2str(k)])
    for n = 1:gfIter
        u = u_nxt;
        
        
        
%         ngu = ((Dcx*u).^2 + (Dcy*u).^2).^(1/2); % Norm of the gradient of u
        ngu = ((Dcx*u).^2 + (Dcy*u).^2).^(1/2) + epsin; % Norm of the gradient of u
        gpxu = ((max(Dpx*u,0)).^2 + (min(Dmx*u,0)).^2 + (max(Dpy*u,0)).^2 + (min(Dmy*u,0)).^2).^(1/2); % Weird positive gradient operators
        K = Dpx*(Dmx*u./ngu) + Dpy*(Dmy*u./ngu);
        hypr = g1.*(Dpx*u) + g2.*(Dmx*u) + g3.*(Dpy*u) + g4.*(Dmy*u);
%         RHS = g.*(K.*ngu + a*gpxu + Dpx*u + Dpy*u); % Wrong
        RHS = g.*K.*(ngu-epsin) + a*g.*gpxu + hypr;
        u_nxt = u + dt*RHS;
        
        C((k-1)*gfIter + n) = norm(u-u_nxt);
        if (C((k-1)*gfIter + n) < tol)
            exitflag = 1;
            break;
        end
        
    end
    
    
    u_nxt = Reinitial2D(reshape(u_nxt,y,x),reinitIter);
    
    if (mod(k,12) == 1)
        clf
        pcolor(flipud(I0)), shading interp, colormap gray, hold on
        contour(flipud(u_nxt),[0 0],'r','linewidth',2);
        drawnow
    end
    
    
    u_nxt = u_nxt(:);

end

%% 
figure(1), clf
u = reshape(u_nxt,y,x);
surf(u), shading interp

figure(2), clf
pcolor(flipud(I0)), shading interp, colormap gray, hold on
contour(flipud(u),[0 0],'r','linewidth',5);
drawnow




