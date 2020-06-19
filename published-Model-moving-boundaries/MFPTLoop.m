function MFPTLoop(n)

if (nargin==0)
    n = 1;
end

% The default parameters. To be adjusted for each particular Case

% PDE parameters
M = 30;
k1 = 0.1; k2 = 1.0; %kon = k1 (?m2/s), koff = k2 (s?1)
D1 = 0.05; D2 = 0.05; % diffusion in and outside the contact (?m2/s)

% Initial radius (?m) and time (s)

%t0 = 0.5*(0.01/2)^2;
t0 = 0.0001;
R0 = 0.22;

% Initial particle location

epsilon = 0.09;
x0 = (1- epsilon);
y0 = 0;

% Rate of area increase

dAdt = 0.0;

params = [k1, k2, D1, D2, dAdt, R0, t0, x0, y0, M];

para = 0;
if (para)
    c=parcluster();
    tmp=tempname();
    mkdir(tmp);
    c.JobStorageLocation=tmp;
end

FigCase = '3BExtra' ;
%FigCase = 'PNAS3B' ;
%FigCase = 'PNAS4A' ;
%FigCase = 'PNAS4B' ;
%FigCase = 'PNAS4C' ;
%FigCase = 'PNAS6' ;
%FigCase = 'PNAS5E' ;
%
%


switch FigCase
    case 'PNAS4B'
        %% Figure 4B
        
        R0 = 0.05;
        
        dAdt = linspace(0,0.03,100);
        t = linspace(0,118,100);
        dvals = [0.05, 0.027];
        
        % Total calculation 80000
        % 50 runs with 1600 each
                                
        [DADT,T,DVALS] = meshgrid(dAdt,t,dvals);
                
        R = sqrt(R0^2 + DADT.*T/pi);
        
        DADT = DADT(:);
        R = R(:);
        DVALS = DVALS(:);

        m = (400*(n-1)+1):(400*n);
        
        Ps = zeros(1,numel(m));
        
        params(10) = 0; % Set M
        paramSet = repmat(params,[numel(m),1]);
        
        paramSet(:,3) = DVALS(m);
        paramSet(:,5) = DADT(m);
        paramSet(:,6) = R(m);


        if (para)
            c.NumWorkers = 12;   parpool(c);
            parfor j = 1:numel(m)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
                fprintf(1,' %g \t %12.11f \n',m(j),Ps(j));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:numel(m)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
        end
        
        filename = ['Fig3BRunN=',num2str(n),'.dat'];
        
        fid = fopen(filename,'w');
        for j = 1:numel(m)
            fprintf(fid,' %g \t %12.11f \n',m(j),Ps(j));
        end
        fclose(fid);
    case '3BExtra'
        %% Figure 3BExtra
        
        R0 = 0.05;
        
        dAdt = linspace(0,0.03,160);
        t = linspace(0,180,160);
       
        % Total calculation 25600
        % 50 runs with 512 each
                                
        [DADT,T] = meshgrid(dAdt,t);
                
        R = sqrt(R0^2 + DADT.*T/pi);
        
        DADT = DADT(:);
        R = R(:);
        
        m = (512*(n-1)+1):(512*n);
        
        Ps = zeros(1,numel(m));
        
        params(10) = 0; % Set M
        paramSet = repmat(params,[numel(m),1]);
        
        paramSet(:,5) = DADT(m);
        paramSet(:,6) = R(m);

        if (para)
            c.NumWorkers = 12;   parpool(c);
            parfor j = 1:numel(m)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
                fprintf(1,' %g \t %12.11f \n',m(j),Ps(j));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:numel(m)                
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
        end
        
        filename = ['Fig3BExtraRunN=',num2str(n),'.dat'];
        
        fid = fopen(filename,'w');
        for j = 1:numel(m)
            fprintf(fid,' %g \t %12.11f \n',m(j),Ps(j));
        end
        fclose(fid);
    case 'PNAS4A'
        %% Figure 4A
        dAdt = 0.1;
        
        T = linspace(0,118,160);
        R = sqrt(R0^2 + dAdt.*T/pi);
        
        D_vals = [0.05 0.027];
        m = length(D_vals);
        
        Ps = zeros(1,m);
        
        params(5) = dAdt; % dAdT = 0.1;
        params(10) = 0; % Set M
        params(6) = R(n); % Set R
        
        paramSet = repmat(params,[m,1]);
        
        paramSet(:,3) = D_vals';
        
        % Set to 1 to run the loops in parallel.
        
        if (para)
            c.NumWorkers = 2;   parpool(c);
            parfor j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
        end
        fid = fopen('Fig3BData.dat','a');
        fprintf(fid,' %16.10f \t %12.11f \t %12.11f \n',params(6),Ps(1),Ps(2));
        fclose(fid);
    case 'PNAS3B'
        %% Figure 3B
        t = logspace(0,4,260);
        t = [0:0.01:0.99 t];
        dadt = logspace(-3,0,100);
        
        R = sqrt(R0^2 + dadt(n).*t/pi);
        
        params(10) = 0;
        params(5) = dadt(n); % Set DADT
        
        paramSet = repmat(params,[length(t),1]);
        paramSet(:,6) = R';
        
        Ps = zeros(1, length(t));
        
        if (para)
            
            c.NumWorkers = 12;  parpool(c);
            
            parfor j = 1:length(t)
                Ps(j) = VolumeCheck(paramSet(j,:));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:length(t)
                
                Ps(j) = VolumeCheck(paramSet(j,:));
            end
        end
        
        filename = ['Fig3DRunN=',num2str(n),'S=2.dat'];
        
        fid = fopen(filename,'w');
        for j = 1:length(t)
            fprintf(fid,'%12.11f \t %12.11f \t %12.11f \n',dadt(n), t(j), Ps(j));
        end
        fclose(fid);
    case 'PNAS4C'
        %% Figure 4C

        R0 = 0.05;
        
        r = linspace(R0,5,500);
        dvals = [0.05, 0.027];
        
        [R,DVALS] = meshgrid(r,dvals);
                
        R = R(:);
        DVALS = DVALS(:);

        m = (100*(n-1)+1):(100*n);
        
        Ps = zeros(1,numel(m));
        
        params(10) = 0;
        
        paramSet = repmat(params,[length(m),1]);
        paramSet(:,6) = R(m);
        paramSet(:,3) = DVALS(m);
        
        Ps = zeros(1, length(m));
        
        if (para)
            
            c.NumWorkers = 12;  parpool(c);
            
            parfor j = 1:length(m)
                Ps(j) = VolumeCheck(paramSet(j,:));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:length(m)            
                Ps(j) = VolumeCheck(paramSet(j,:));
            end
        end
        
        filename = ['Fig3ERunN=',num2str(n),'.dat'];
        
        fid = fopen(filename,'w');
        for j = 1:numel(m)
            fprintf(fid,' %g \t %12.11f \n',m(j),Ps(j));
        end
        fclose(fid);
    case '4Cold'
        %% Figure 4C old
        T = logspace(-8,-1,200);
        
        M_vals = [30 300];
        Ps = zeros(1,2);
        
        params(2) = 1/T(n);
        paramSet = repmat(params,[length(M_vals),1]);
        
        paramSet(:,10) = M_vals';
        
        if (para)
            c.NumWorkers = 2;   parpool(c);
            parfor j = 1:length(M_vals)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:))
            end
            delete(gcp('nocreate'))
        else
            for j = 1:length(M_vals)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:))
            end
        end
        
        % Print to file: T Ps(1) Ps(2) Ps(3) Ps(4)
        fid = fopen('Fig4CData.dat','a');
        fprintf(fid,' %16.10f \t %12.11f \t %12.11f \n',params(2),Ps(1),Ps(2));
        fclose(fid);
    case '4CAlt'
        %% Figure 4CAlt
        T = linspace(-8,-1,200);        
        M_vals = linspace(10,300,200);
        Ps = zeros(1,2);
        
        params(2) = 1/T(n);
        paramSet = repmat(params,[length(M_vals),1]);
        
        paramSet(:,10) = M_vals';
        
        if (para)
            c.NumWorkers = 2;   parpool(c);
            parfor j = 1:length(M_vals)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:))
            end
            delete(gcp('nocreate'))
        else
            for j = 1:length(M_vals)
                Ps(j) = MFPTGrowingDisk(paramSet(j,:))
            end
        end
        
        % Print to file: T Ps(1) Ps(2) Ps(3) Ps(4)
        fid = fopen('Fig4CData.dat','a');
        fprintf(fid,' %16.10f \t %12.11f \t %12.11f \n',params(2),Ps(1),Ps(2));
        fclose(fid);
    case 'PNAS6'
        %% Figure 6
        
        % Run with n = 6.
        
        kon = 100*[1.2e-2, 3.5e-3, 4.0e-6, 4.4e-6, 2.7e-5 1.3e-5];
        koff = [10.8, 4.7, 1.4, 2.6, 1.3, 2.4];
        M = [30 30 300 300 30 30];
        
        params(1) = kon(n);
        params(2) = koff(n);
        params(10) = M(n);
        
        Ps = MFPTGrowingDisk(params);
        
        % Print to file: T Ps(1) Ps(2) Ps(3) Ps(4)
        fid = fopen('Fig4DDataKON1Hun.dat','a');
        fprintf(fid,' %16.10f \t %12.11f \n',n,Ps);
        fclose(fid);
    case '4E'
        %% Figure 4E
        R = linspace(0.1,2,100);
        M_vals = [0.0 30];
        Ps = zeros(1,2);
        
        m = length(M_vals);
        
        params(6) = R(n);
        paramSet = repmat(params,[m,1]);
        paramSet(:,10) = M_vals';
        if (para)
            c.NumWorkers = 2;
            parpool(c);
            parfor j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
        end
        
        fid = fopen('Fig4EData.dat','a');
        fprintf(fid,' %16.10f \t %12.11f \t %12.11f \n',params(6),Ps(1),Ps(2));
        fclose(fid);
    case 'PNAS5E'
        %% Figure 5E
        R  = logspace(-3,0,400);
        M_vals = [0 30 30 30 30 30];
        Koff_Vals = [0 1 5 10 20 50];
        
        m = length(M_vals);
        Ps = zeros(1,m);
        
        params(6) = R(n);
        paramSet = repmat(params,[m,1]);
        paramSet(:,10) = M_vals';
        paramSet(:,2) = Koff_Vals';
        
        if (para)
            c.NumWorkers = 6;
            parpool(c);
            parfor j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
            delete(gcp('nocreate'))
        else
            for j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:));
            end
        end
        
        fid = fopen('Fig4FData.dat','a');
        fprintf(fid,' %16.10f \t %12.11f \t %12.11f \t %12.11f \t %12.11f \t %12.11f\n',R(n),Ps(1),Ps(2),Ps(3),Ps(4),Ps(5),Ps(6));
        fclose(fid);
    case 'Range'
        %% Range
        m = 2;
        kon_r = linspace(1e-5,1e-1,m);
        koff_r = linspace(0.5,50,m);
        
        [KON,KOFF] = meshgrid(kon_r,koff_r);
        
        GridKON = reshape(KON,m^2,1);
        GridKOFF = reshape(KOFF,m^2,1);
        
        k1 = GridKON(1+m*(n-1):m*n);
        k2 = GridKOFF(1+m*(n-1):m*n);
        
        paramSet = repmat(params,[m,1]);
        paramSet(:,1) = reshape(k1,m,1);
        paramSet(:,2) = reshape(k2,m,1);
        
        Ps = zeros(m,1);
        if (para)
            c.NumWorkers = 2;    parpool(c);
            parfor j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:))
            end
            delete(gcp('nocreate'))
        else
            for j = 1:m
                Ps(j) = MFPTGrowingDisk(paramSet(j,:))
            end
        end
        
        for j= 1:m
            fid = fopen('FigRangeData.dat','a');
            fprintf(fid,' %16.10f \t %12.11f \t %12.11f \n',paramSet(j,1),paramSet(j,2),Ps(j));
        end
        fclose(fid);
end



