file0 = 'tmpm01_0000.h5';
file1 = 'tmpm01_0099.h5';

ndiv = 50;
Lx = 1.0;
Lz = 0.16;
K  = 10.0e4;
F  = 1.0e1;
Nu = 0.3;
E  = (3.0*(1.0-2.0*Nu))*K;
Dx = 1.5*Lx/ndiv;


Pos0 = h5read(file0,'/Position');
P0x = Pos0(1:3:end);
P0y = Pos0(2:3:end);
P0z = Pos0(3:3:end);

Pos1 = h5read(file1,'/Position');
P1x = Pos1(1:3:end);
P1y = Pos1(2:3:end);
P1z = Pos1(3:3:end);

idx =  find(P0z<0.6*Lz & P0z>0.4*Lz);

x = P0x(idx)-Dx;
w = P0z(idx)-P1z(idx);
wt= 2.0*F*x.^2.*(3.0*(Lx-Dx)-x)/(E*Lz^4);

plot(x,w,'o',x,wt,'-');