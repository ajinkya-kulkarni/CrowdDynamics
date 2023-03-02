clc
clear all
% NS = zeros(1000,2);
for filenum=1:2000

cd particle_info
if(filenum<10)
    fname = strcat('part_data_000',int2str(filenum));
elseif(filenum<100)
    fname = strcat('part_data_00',int2str(filenum));
elseif(filenum<1000)
    fname = strcat('part_data_0',int2str(filenum));
else
    fname = strcat('part_data_',int2str(filenum));
end 
A = load(fname,'-ascii');
cd ..

cd particle_vel_info

if(filenum<10)
    fname = strcat('part_vel_data_000',int2str(filenum));
elseif(filenum<100)
    fname = strcat('part_vel_data_00',int2str(filenum));
elseif(filenum<1000)
    fname = strcat('part_vel_data_0',int2str(filenum));
else
    fname = strcat('part_vel_data_',int2str(filenum));
end
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
tx = 0:dx:L;
ty = 0:dy:L;
m= L/dx;n=L/dy;
[xq,yq] = meshgrid(tx,ty);
vxq = Fx(xq,yq);
vyq = Fy(xq,yq);

for j = 2:n-1
	for i = 2:m-1			
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
for j = 2:n-1
	for i = 2:m-1
		Nstot = Nstot + NsViscous(i,j)*dx*dy;		
	end
end
NS(filenum,1)=filenum;
NS(filenum,2)=Nstot;
clearvars -except filenum NS
%%  Stream function 
% psi = zeros(m,n);
% tecplot = zeros(m*n,3);
% %Calculating Stream Function
% for i=1:m    
%     for j=1:n-1        
%         psi(i,j+1) = psi(i,j)+(Vx(i,j)*dy);
%     end
%     
% end
% 
% contour(psi);
% %Formatting data for Tecplot
% k = 0;
% 
% for j=1:n  
% 
%     for i=1:m        
%         k = k+1;
%         tecplot(k,1) = (i-1)*dx;
%         tecplot(k,2) = (j-1)*dy;
%         tecplot(k,3) = psi(i,j);
%     end
% end
% 
% %Add the Header Data in .csv format & then save it in .dat format
% %Header Data:    Zone I nx J ny
% xlswrite('tecplot.csv',tecplot);
end
xlswrite('NS.csv',NS);
