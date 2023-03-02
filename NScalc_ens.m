clc
clear all
% NS = zeros(1000,2);
index = 0;
for filenum=1000:2000
index = index +1;
cd particle_info
fname = strcat('part_data_',int2str(filenum)); 
A = load(fname,'-ascii');
cd ..

cd particle_vel_info
fname=strcat('part_vel_data_',int2str(filenum)); 
B = load(fname,'-ascii');
cd ..

C=[[1:6084]', A(1:6084,2:3), B(1:6084,2:3)];

x = C(:,2);
y = C(:,3);
vx = C(:,4);
vy = C(:,5);
Fx = scatteredInterpolant(x,y,vx);
Fy = scatteredInterpolant(x,y,vy);
dx = 0.01;dy = 0.01; L = 0.4;
m= L/dx;n=L/dy;
tx = 0:dx:L;
ty = 0:dy:L;
[xq,yq] = meshgrid(tx,ty);
vxq = Fx(xq,yq);
vyq = Fy(xq,yq);

for j = 1:n
    for i = 1:m	
        U(index,i,j)=vxq(i,j);
        V(index,i,j)=vyq(i,j);
    end
end

clearvars -except index U V

end
% averaging of velocities
dx = 0.01;dy = 0.01; L = 0.4;
m= L/dx;n=L/dy;
X=[dx/2:dx:L]';
Y=[dy/2:dy:L]';
vxq = zeros(m+1,n+1);
vyq = zeros(m+1,n+1);

for j = 1:n
    for i = 1:m	
        for ind =1:index
            vxq(i,j)=vxq(i,j)+ U(ind,i,j);
            vyq(i,j)=vyq(i,j)+V(ind,i,j);
        end
    end
end
vxq= vxq./index;vyq=vyq./index;
%dissipation calculation

for j = 2:n
	for i = 2:m			
		dU = vxq(i+1,j) - vxq(i,j);
		dV = vyq(i,j+1) - vyq(i,j);
		dUdX = dU/dx; 	
		dVdY = dV/dy; 
		%%******************************** below  dUdY
		uuP  = (vxq(i+1,j)+vxq(i,j))/2;

		uS  = (vxq(i+1,j-1)+vxq(i,j-1))/2;
		Usf = (uuP + uS)/2;		

		uuN  = (vxq(i+1,j+1)+vxq(i,j+1))/2;
		Unf = (uuP + uuN)/2;
		
		dUdY = (Unf-Usf)/dy;
        %%******************************** below dVdX
	
		vP  = (vyq(i,j+1)+vyq(i,j))/2;

		vW  = (vyq(i-1,j+1)+ vyq(i-1,j))/2;
		Vwf = (vP + vW)/2;
		

		vE  = (vyq(i+1,j+1) + vyq(i+1,j))/2;
		Vef = (vP + vE)/2;
		dVdX = (Vef-Vwf)/dx;
		NsViscous(i,j) = (2*((dUdX^2)+ (dVdY^2))+(dUdY+dVdX)^2);
	end
end
Nstot = 0;
for j = 2:n
	for i = 2:m
		Nstot = Nstot + NsViscous(i,j)*dx*dy;		
	end
end
fid=fopen('NS_ens.txt','wt');
fprintf(fid,'%8.5f \n',Nstot);
fclose(fid);
% time averaged velocity
VXY = [];
%VXY = zeros(m*n,7);
for j = 1:n
	for i = 1:m
		VXY = [VXY;[i j X(i,1) Y(j,1) vxq(i,j) vyq(i,j)] sqrt(vxq(i,j)^2 + vyq(i,j)^2)];	
	end
end
xlswrite('VXY.csv',VXY);

