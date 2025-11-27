Lx = 5.2;
Ly = 5.2;
R  = 0.1*0.018;
totparvol = 128.455;

df = dir("*tper_b_0*.h5");
filename = df(end).name;
P = h5read(filename,'/Position');
Px= P(1:3:end);
Py= P(2:3:end);
Pz= P(3:3:end);
V = h5read(filename,'/PVelocity');
Vx= V(1:3:end);
Vy= V(2:3:end);
Vz= V(3:3:end);

figure(1);
plot(Vx,Pz,'o');
hold on;
 
db = dir("*tper_b_bf_0*.h5");
Mu = [];
P  = [];
Phi= [];
for i=1:length(db)
    
    filename = db(i).name;
    B = h5read(filename,'/Branch');
    Bx= B(1:3:end);
    By= B(2:3:end);
    Bz= B(3:3:end);
    
    N = h5read(filename,'/Normal');
    Nx= N(1:3:end);
    Ny= N(2:3:end);
    Nz= N(3:3:end);
    
    T = h5read(filename,'/Tangential');
    Tx= T(1:3:end);
    Ty= T(2:3:end);
    Tz= T(3:3:end);
    
%     I2 = h5read(filename,'/ID2');
%     idx = find((I2~=length(Px)-1)&(I2~=length(Px)-2));
%     Bx = Bx(idx);
%     By = By(idx);
%     Bz = Bz(idx);
%     Nx = Nx(idx);
%     Ny = Ny(idx);
%     Nz = Nz(idx);
%     Tx = Tx(idx);
%     Ty = Ty(idx);
%     Tz = Tz(idx);
      
    filename = df(i).name;
    Pos = h5read(filename,'/Position');
    Px  = Pos(1:3:end);
    Py  = Pos(2:3:end);
    Pz  = Pos(3:3:end);
    
    Lz = max(Pz)-min(Pz)-2*R;      
    Vol = Lx*Ly*Lz;
    S = (1/Vol)*[sum(Bx.*(Nx+Tx)) sum(Bx.*(Ny+Ty)) sum(Bx.*(Nz+Tz));sum(By.*(Nx+Tx)) sum(By.*(Ny+Ty)) sum(By.*(Nz+Tz));sum(Bz.*(Nx+Tx)) sum(Bz.*(Ny+Ty)) sum(Bz.*(Nz+Tz))];
    
    mu = -3*sqrt(S(3,1)^2+S(1,3)^2)/trace(S)/sqrt(2);
    
    Mu = [Mu mu];
    
    P  = [P  -trace(S)/3];
    
    Phi= [Phi totparvol/Vol];
    
end

figure(2)
plot(Mu);
hold on;

figure(3);
plot(P);
hold on;

figure(4);
plot(Phi);
hold on;
