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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

percentage=10;

k1=(percentage/100)*nint;

kk=round(k1);

gamma2 = dlmread('gamma_input.txt');

%%%%% isolation of particles in the strip defined by r to r+delta_r %%%%%

r = sqrt((x-0.4).*(x-0.4) + (y-0.4).*(y-0.4));
R = max(r);

temp = dlmread('radial_locs.txt');
radial_locs = temp*R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
for width=1:1:1000
    
Imod1 = find((r>(radial_locs-width*d)) & (r<(radial_locs+width*d)));

    if (length(Imod1) > kk)
        break
    end
    
end

Imod2 = find((r>(radial_locs-(width+0.01)*d)) & (r<(radial_locs+(width+0.01)*d)));


%%%%%%%%%%%%%%%%%%%%%%%% gamma %%%%%%%%%%%%%%%%%%%%%%%%

gamma = (zeros(nint,1));

%gamma(randsample(nint,kk)) = gamma2;

gamma(Imod2) = gamma2;

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


