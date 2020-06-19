function out = VolumeCheck(params)

if (nargin == 0)
    
    % PDE parameters
    M = 30;
    k1 = 0.1; k2 = 1;
    D1 = 0.05; D2 = 0.05;
    
    % Initial radius and time
    
    %t0 = 0.5*(0.01/2)^2;
    t0 = 0.0001;
    R0 = 0.22;
    
    % Initial particle location
    
    epsilon = 0.09;
    x0 = (1- epsilon);
    y0 = 0;
    
    % Rate of area increase
    
    dAdt = 0.00;
    
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

R = @(s) sqrt(R0^2 + dAdt.*s/pi);

tf = 2;

doPlot = 0;

% tc1 = 2.5e-4;
% tc2 = 0.004;
% 
% t1 = linspace(t0,tc1,500);
% t2 = linspace(tc1,tc2,2000);
% t3 = linspace(tc2,tf,2000);
% t = [ t1 t2(2:end) t3(2:end) ];

t = [t0 tf];

model = createpde(1);
geometryFromEdges(model,@b);

setInitialConditions(model,@(s) initfun(s,params));

generateMesh(model,'Hmax',0.02);

specifyCoefficients(model,'m',0,'d',1,'c',@(r,s) ccoeff(r,s,params),'a',0,'f',@(r,s) fcoeff(r,s,params));

t = [ t0 tf ];

P = zeros(size(t));

%applyBoundaryCondition(model,'edge',[1:4],'q',@QFun);
applyBoundaryCondition(model,'edge',[1:4],'u',0);

sol = solvepde(model,t);

u = sol.NodalSolution;

for j = length(t)
    
    if (doPlot)
        
        pdeplot(model,'xydata',u(:,j),'zdata',u(:,j));
        colormap('jet');
        
    end
    
    P(j) = GaussIntNew(model,R(t(j)).^2 * u(:,j));
    %drawnow;    
    
end

out = P(end);
%plot(t,P)

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

u0 = zeros(1,length(x));

u0(1,:) =  exp( -( (x-x0).^2 + (y-y0).^2 )/(4*t0))/(4*t0 * pi);
u0 = u0/R0^2;


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

%params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0];

%params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

D1 = params(3);
dAdt = params(5);
R0 = params(6);

nr = length(region.x);

t = state.time;

R = (R0^2 + dAdt*t/pi)^(1/2);

ft = 1/(R)^2;

c = repmat([ft*D1],[1,nr]);

end

function f = fcoeff(region,state,params)

% params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

dAdt = params(5);
R0 = params(6);

N = 1;
nr = length(region.x);
f = zeros(N,nr);

ux = state.ux(1,:);
uy = state.uy(1,:);

t = state.time;

x = region.x;
y = region.y;

R = (R0^2 + dAdt*t/pi)^(1/2);
Rp = (0.5*dAdt)/(pi*R);

ft = Rp/R;

f(1,:) = ft*(x.*ux + y.*uy);

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


function bcMatrix = QFun(region,state)

global params

%params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0];

t = state.time;
nr = length(region.x);

dAdt = params(5);
R0 = params(6);

R = (R0^2 + dAdt*t/pi)^(1/2);
Rp = (0.5*dAdt)/(pi*R);
ft = Rp/R;

bcMatrix = repmat([ft],[1,nr]);

end