function out = MFPTGrowingDisk(params)

if (nargin == 0)
    
    % PDE parameters
    M = 10;
    k1 = 1; k2 = 0.1;
    D1 = 0.05; D2 = 0.05;
    
    % Initial radius and time
    
    %t0 = 0.5*(0.01/2)^2;
    t0 = 0.0005;
    R0 = 0.22;
    
    % Initial particle location
    
    epsilon = 0.1;
    x0 = (1- epsilon);
    y0 = 0;
    
    % Rate of area increase
    
    dAdt = 0.0;
    
    params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];
    
end

k1 = params(1);
k2 = params(2);
D1 = params(3);
D2 = params(4);
dAdt = params(5);
R0 = params(6);
t0 = params(7);
x0 = params(8);
y0 = params(9);
M = params(10);

if (k1*M==0)
    out = VolumeCheck(params);
    return
end

R = @(s) sqrt(R0^2 + dAdt.*s/pi);

tf = 2;

doPlot = 1;

tc1 = 2*t0;
tc2 = 50*t0;

t1 = linspace(t0,tc1,400);
t2 = linspace(tc1,tc2,400);
t3 = linspace(tc2,tf,1500);
t = [ t1 t2(2:end) t3(2:end) ];

%t = [t0 tf];

P = zeros(size(t));

model = createpde(2);
geometryFromEdges(model,@b);

setInitialConditions(model,@(s) initfun(s,params));

applyBoundaryCondition(model,'edge',[1:4],'q',@(r,s) QFun(r,s,params),'g',[0;0]);
applyBoundaryCondition(model,'edge',[1:4],'u',0,'EquationIndex',[1]);

%applyBoundaryCondition(model,'mixed','Edge',[1:4],'u',0,'EquationIndex',[1],'q',@QFun,'g',0);

specifyCoefficients(model,'m',0,'d',1,'c',@(r,s) ccoeff(r,s,params),'a',-[-M*k1 M*k1 k2 -k2]','f',@(r,s) fcoeff(r,s,params));
generateMesh(model,'Hmax',0.015);

%[p,e,~] = meshToPet(model.Mesh);
%indexBdy = find(e(7,:)==0);

sol = solvepde(model,t);
u = sol.NodalSolution;

% Get index of boundary points.

% Solve the PDE.
%
% xq = p(1,indexBdy); yq = p(2,indexBdy);
% [th, I] = sort(pi - atan2(yq,xq));

%d  = (p(1,:)-x0).^2 + (p(2,:)-y0).^2;
%u0 =  (1/R0^2) * exp(-d/(4*t0))/(4*t0 * pi);
%u0 = u0/GaussIntNew(model, R0^2 *u0');

%GaussIntNew(model, R0^2 *u0')


% uBdy = u(indexBdy,1,end);
% uBdy = uBdy(I);

% Evaluate the solution and its gradient on the boundary.
%
% [gradx,grady] = evaluateGradient(sol,[xq;yq],1,1:length(t));
%
% C = zeros(size(t));
%
% for n = 1:length(t)
%
%     ker = xq.*reshape(gradx(:,n),size(xq)) + yq.*reshape(grady(:,n),size(xq));
%
%     C(n) = -D1*trapz([th (2*pi+th(1))],[ker ker(1)]);
%
% end
% if (doPlot)
%
%     figure(1);
%     hold on
%     plot(t,C);
%     trapz(t,C)
%     hold off
%
% end
%close all;
%figure('position', [0, 0, 1000, 400]);

if (doPlot)
   close all; 
end

for j = 1:length(t)
    P(j) = GaussIntNew(model, R(t(j))^2 *(u(:,1,j)+u(:,2,j)));
      
%     if (doPlot)
%         
%         subplot(1,2,1), pdeplot(model,'xydata',u(:,1,j),'zdata',u(:,1,j));
%         title(['Free mass = ', P(j)]);
%         
%         subplot(1,2,2), pdeplot(model,'xydata',u(:,2,j),'zdata',u(:,2,j));
%         colormap('jet');
%         
%         drawnow;
%     end
end

P = [ P];
t = [ t];

out = P(length(t));

if (doPlot)
    
    figure('color','w')
    hold on
    plot(t,P,'-.k')
    box on;
    set(gca,'Fontsize',24);
    xlabel('$t$','Fontsize',28,'Interpreter','latex');
    ylabel('$P(t)$','Fontsize',28,'Interpreter','latex','rotation',0);
    
    xs = spline(t,-P);
    C = ppval(fnder(xs,1),t);
    xlim([0 1e-2])
    
    figure('color','w')
    plot(t,C,'-k');
    box on;
    set(gca,'Fontsize',24);
    xlabel('$t$','Fontsize',28,'Interpreter','latex');
    ylabel('$C(t)$','Fontsize',28,'Interpreter','latex','rotation',0);
    hold off
    CaptureMass = trapz(t,C);
    xlim([0 1e-2])
    
end

save('CaptureData.mat','t','C');

end

function bcMatrix = QFun(region,state,params)

%params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

t = state.time;
nr = length(region.x);

dAdt = params(5);
R0 = params(6);

R = (R0^2 + dAdt*t/pi)^(1/2);
Rp = (0.5*dAdt)/(pi*R);
ft = Rp/R;

bcMatrix = repmat([0 0;0 ft],[1,nr]);

end

function u0 = initfun(s,params)

%params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

x = s.x; y = s.y;

k1 = params(1);
k2 = params(2);
M = params(10);

R0 = params(6);
t0 = params(7);
x0 = params(8);
y0 = params(9);

M = length(x);
u0 = zeros(2,M);

d  = (x-x0).^2 + (y-y0).^2;

%u0(1,:) =  (1/R0^2) *(k2/(M*k1+k2)) *exp( -( d )/(4*t0))/(4*t0 * pi);
%u0(2,:) = (1/R0^2) *(M*k1/(M*k1+k2))*exp( - ( d )/(4*t0))/(4*t0 * pi);

u0(1,:) =  (1/R0^2) * exp(-d/(4*t0))/(4*t0 * pi);

end


function [x,y] = b(bs,s)

% Define the boundary

%CARDG Geometry File defining disc with a hole.
nbs= 4;

if nargin==0
    x=nbs;
    return
end

dl = zeros(4,nbs);
dl(1,1:end) = linspace(pi/2,2*pi,nbs) - pi/2;
dl(2,1:end) = linspace(pi/2,2*pi,nbs);
dl(3,1:4) = 1;

if nargin==1
    x=dl(:,bs);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 & n==1,
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) | n~=size(s,2),
    error('bs must be scalar or of same size as s');
end

nth=200;
th=linspace(0,2*pi,nth);

ii=find(bs<=4);

if ~isempty(ii)
    xt=cos(th);
    yt=sin(th);
    th=pdearcl(th,[xt;yt],s(ii),0,2*pi);
    x(ii) = cos(th);
    y(ii) = sin(th);
    
end

end

function c = ccoeff(region,state,params)

%params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

D1 = params(3);
D2 = params(4);
dAdt = params(5);
R0 = params(6);

nr = length(region.x);

t = state.time;

R = (R0^2 + dAdt*t/pi)^(1/2);

ft = 1/(R)^2;

%c = repmat([ft*D1; ft*D2],[1,nr]);

c = zeros(2,nr);
c(1,:) = ft*D1;
c(2,:) = ft*D2;

end

function f = fcoeff(region,state,params)

% params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

dAdt = params(5);
R0 = params(6);

N = 2;
nr = length(region.x);
f = zeros(N,nr);

if (dAdt ~= 0)
    
    ux = state.ux(1,:);
    uy = state.uy(1,:);
    vx = state.ux(2,:);
    vy = state.uy(2,:);
    
    t = state.time;
    
    x = region.x;
    y = region.y;
    
    R = (R0^2 + dAdt*t/pi)^(1/2);
    Rp = (0.5*dAdt)/(pi*R);
    
    ft = Rp/R;
    
    f(1,:) = ft*(x.*ux + y.*uy);
    f(2,:) = ft*(x.*vx + y.*vy);
    
end

end

function out = GaussIntNew(model,u)

p = model.Mesh.Nodes;
t = model.Mesh.Elements;

u = u';

it1=t(1,:);
it2=t(2,:);
it3=t(3,:);
Areas = 0.5* abs (p(1,it1).*(p(2,it2)-p(2,it3)) +...
    p(1,it2).*(p(2,it3)-p(2,it1))+ ...
    p(1,it3).*(p(2,it1)-p(2,it2)));
fun = @(s) s;

out = (1/3) * sum(Areas.*(fun(u(it1))+fun(u(it2))+fun(u(it3))));

end

