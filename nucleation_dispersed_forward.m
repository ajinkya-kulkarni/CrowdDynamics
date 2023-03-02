clc;
clear all;

rng('shuffle')

xi=0.4;
yi=0.4;

nint=6120;
total=7289;
d=0.005;

m=dlmread('part_data_2000');
x=m((1:nint),2);
y=m((1:nint),3);

n=dlmread('part_vel_data_2000');
vx=n((1:nint),2);
vy=n((1:nint),3);

m1=dlmread('part_data_2000');
xwall=m1((nint+1:total),2);
ywall=m1((nint+1:total),3);

%%%%%%%%%%%%%%%%%%%%%%%% gamma %%%%%%%%%%%%%%%%%%%%%%%%

percentage=10;

k1=(percentage/100)*nint;

kk=round(k1);

gamma2 = dlmread('gamma_input.txt');

gamma = (zeros(nint,1));

gamma(randsample(nint,kk)) = gamma2;

gamma(nint+1:total) = 0;

temp_beta = horzcat((0:total-1)',gamma);

dlmwrite('gamma.txt', temp_beta,'delimiter','\t','precision',8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f2 = figure('visible', 'off');
hold on
scatter(x,y,[],gamma(1:nint),'filled')
scatter(xwall,ywall,'filled')
colormap(flipud(gray))
colorbar
axis equal
box on

print -depsc scatter_forward.eps
close(f2)


