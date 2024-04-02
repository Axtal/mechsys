deltatheta =pi/5;
deltaphi=pi/5;
Theta=0:deltatheta:pi;
Phi=0:deltaphi:2*pi;
[phi,theta]=meshgrid(Phi,Theta);
r = sin(theta);
X=r.*cos(phi);
Y=r.*sin(phi);
Z=cos(theta);
figure(1);
surf(X,Y,Z);
hold on;

filename = 'ttt_b_bf_0000.h5';
Info = hdf5info(filename);                                %   read particle information from data
Fnormal = hdf5read(Info.GroupHierarchy.Datasets(5));       %   Normal force vector
N = Info.GroupHierarchy.Datasets(5).Dims/3;
Fnx=zeros(N,1);
Fny=zeros(N,1);
Fnz=zeros(N,1);
idx=1;
R = zeros(size(r,1),size(r,2));
for i=1:N
    Fnx(i)=Fnormal(idx);
    Fny(i)=Fnormal(idx+1);
    Fnz(i)=Fnormal(idx+2);
    idx = idx+3;
    % calculate angle theta1 phi1
    phi1 = atan2(Fny(i),Fnx(i))+pi;
    theta1 = acos(Fnz(i)/(Fnx(i)^2+Fny(i)^2+Fnz(i)^2)^0.5);
    %         [~, j] = min(abs(Phi-phi1));
    %         [~, k] = min(abs(Theta-theta1));
    k = floor(phi1/deltaphi)+1;
    j = floor(theta1/deltatheta)+1;
    R(j,k)=R(j,k)+1;
    % duplicate force;
    phi1 = atan2(-Fny(i),-Fnx(i))+pi;
    theta1 = acos(-Fnz(i)/(Fnx(i)^2+Fny(i)^2+Fnz(i)^2)^0.5);
    k = floor(phi1/deltaphi)+1;
    j = floor(theta1/deltatheta)+1;
    R(j,k)=R(j,k)+1;
end
%divide on area to get the frequency

R = R(1:end-1,1:end-1);

for j=1:size(R,1)
    for k=1:size(R,2)
        R(j,k) = R(j,k)/(cos(Theta(j))-cos(Theta(j+1)))/deltaphi;
    end
end

thetap = theta(1:end-1,1:end-1)+0.5*deltatheta;
phip   = phi  (1:end-1,1:end-1)+0.5*deltaphi;

spherobar(R,thetap,phip);
axis square;


% x=R.*sin(thetap).*cos(phip);
% y=R.*sin(thetap).*sin(phip);
% z=R.*cos(thetap);
% figure(2);
% surf(x,y,z);
% hold on;