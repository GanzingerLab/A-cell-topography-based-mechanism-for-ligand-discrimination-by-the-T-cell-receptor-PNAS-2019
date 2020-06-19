function out = VolumeCheckMFPT(params)

% Output the Mean.

if (nargin == 0)
    
    % PDE parameters
    M = 0;
    k1 = 0.0; k2 = 0.25;
    D1 = 0.05; D2 = 0.05;
    
    % Initial radius and time
    
    %t0 = 0.5*(0.01/2)^2;
    t0 = 0.005;
    R0 = 0.22;
    
    % Initial particle location
    
    epsilon = 0.49;
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

R = @(s) sqrt(R0^2 + dAdt.*s/pi);

tf = 2;

doPlot = 0;
make_movie = 0;

if ( make_movie )
    vidObj = VideoWriter('FullDensity.avi');
    vidObj.FrameRate=30;
    open(vidObj);
end

tc1 = 2.5e-3;
tc2 = 0.004;

t1 = linspace(t0,tc1,500);
t2 = linspace(tc1,tc2,1000);
t3 = linspace(tc2,tf,1000);

t = [ t1 t2(2:end) t3(2:end) ];

t = linspace(t0,tf,400);

model = createpde(1);
geometryFromEdges(model,@b);

setInitialConditions(model,@(s) initfun(s,params));

generateMesh(model,'Hmax',0.04);

specifyCoefficients(model,'m',0,'d',1,'c',@(r,s) ccoeff(r,s,params),'a',0,'f',@(r,s) fcoeff(r,s,params));

%t = [ t0 tf ];

P = zeros(size(t));

[p0,e,T] = meshToPet(model.Mesh);

%applyBoundaryCondition(model,'edge',[1:4],'q',@QFun);
applyBoundaryCondition(model,'edge',[1:4],'u',0);

sol = solvepde(model,t);

u = sol.NodalSolution;

for j = 1:length(t)
    p = sqrt(R0^2 + dAdt*t(j)/pi) *p0;
    P(j) = GaussIntNew(model,R(t(j)).^2 * u(:,j));
    if (doPlot)
        %pdeplot(model,'xydata',u(:,j),'zdata',u(:,j));
        pdesurf(p,T,u(:,j));
        colormap('jet'); box on;
        
        xlim(R0*[-2 2]); ylim(R0*[-2 2]);
        zlim([0 50]); box on;
        title(['$t = $ ',num2str(t(j),'%4.3f'), ',  Remaining Mass = ',num2str(P(j),'%6.5f')],'Interpreter','latex','fontsize',18);
        xlabel('$x$','Interpreter','latex','fontsize',18);
        ylabel('$y$','Interpreter','latex','fontsize',18);
        zlabel('$u$','Interpreter','latex','rotation',0,'fontsize',18);
        set(gca,'Xtick',R0*[-2:2],'XTickLabel',{'-2R_0','-R_0','0','R_0','2R_0'})
        set(gca,'Ytick',R0*[-2:2],'YTickLabel',{'-2R_0','-R_0','0','R_0','2R_0'})
        
        %view((j/10),40);
        drawnow;
        
    end
    
    if(make_movie)
        %fr = getframe;
        fr = export_fig(gcf,'-transparent');
        writeVideo(vidObj,fr);
    end
    
end

close all

figure('color','w');

hold on;
xs = spline(t,-P);
C = ppval(fnder(xs,1),t);
mean = simps(t,t.*C);
out = P(end);
plot(t,C,'b','linewidth',2);

plot(mean,ppval(fnder(xs,1),mean),'ks','markerfacecolor','k')
box on;
xlabel('$t$','Interpreter','latex','fontsize',18);
ylabel('Occupancy time distribution','Interpreter','latex','fontsize',18);


hold off;

if(make_movie)
    close(vidObj);
end

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