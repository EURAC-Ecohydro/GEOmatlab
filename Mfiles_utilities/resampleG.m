function B = resampleG(A,nx)

%function B = resampleG(A,nx)
% Copyright, 2009 Giacomo Bertoldi
%
% Given a array of columns A of Nx rows
% we obtain a matrix B averaged
% at resolution Nx/nx
% limits: Nx/nx should be integers

[Nx, Ny] = size(A);

Ni=floor(Nx/nx);

if (Ni ~= Nx/nx)
    disp 'Warning Nx/nx not integer!'
end

B=ones(Ni,Ny);

%Takes average of block and stores it

for x = nx:nx:Nx
       B(x/nx,:) = mean(A(x-nx+1:x,:));   
end

