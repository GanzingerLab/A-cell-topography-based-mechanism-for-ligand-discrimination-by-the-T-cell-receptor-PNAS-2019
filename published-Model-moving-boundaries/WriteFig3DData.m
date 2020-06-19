function WriteFig3DData
%new PNAS figure 3B
n = 36;

t  = linspace(0,118,250);
dadt = logspace(-2,1,36);

[T,DADT] = meshgrid(t,dadt);

Ps = zeros(size(T));

for j = 1:n
    filename = ['./Fig3DData/Fig3DRunN=',num2str(j),'S=2-1.dat'];
    fid = fopen(filename,'r');
    A = fscanf(fid,'%f %f %f',[3,Inf])';
    fclose(fid);

    Ps(j,:) = A(:,3);
      
end
