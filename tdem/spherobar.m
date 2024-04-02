function y=spherobar(R,thetap,phip)
deltaphi   = phip  (1,2)-phip  (1,1);
deltatheta = thetap(2)-thetap(1);
figure;
for i=1:size(R,1)
    for j=1:size(R,2)
        plotphi = [phip(i,j)-0.5*deltaphi phip(i,j)+0.5*deltaphi phip(i,j)+0.5*deltaphi phip(i,j)-0.5*deltaphi];
        plottheta = [thetap(i,j)-0.5*deltatheta  thetap(i,j)-0.5*deltatheta  thetap(i,j)+0.5*deltatheta  thetap(i,j)+0.5*deltatheta];
        plotr = [R(i,j) R(i,j) R(i,j) R(i,j)];
        plotx = plotr.*sin(plottheta).*cos(plotphi);
        ploty = plotr.*sin(plottheta).*sin(plotphi);
        plotz = plotr.*cos(plottheta);
        fill3(plotx,ploty,plotz,'blue');
        hold on;
        
        trix  = [0 0 0 0;plotx(1) plotx(2) plotx(3) plotx(4);plotx(2) plotx(3) plotx(4) plotx(1)];
        triy  = [0 0 0 0;ploty(1) ploty(2) ploty(3) ploty(4);ploty(2) ploty(3) ploty(4) ploty(1)];
        triz  = [0 0 0 0;plotz(1) plotz(2) plotz(3) plotz(4);plotz(2) plotz(3) plotz(4) plotz(1)];
        fill3(trix,triy,triz,'blue');
    end
end