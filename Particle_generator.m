% clc;
% clear all;
% 
% delete init_positions.txt;
% delete init_velocities.txt;
% 
% xi=1.0%input('x co-ord of the cavity:');
% yi=1.0%input('y co-ord of the cavity:');
% r0=0.0025%input('radius of the ball:');
% d0=2*r0;
% R=0.8;%input('Radius of the cavity:');
% 
% p=[];
% q=[];   
% x=[];
% y=[];
% xw=[];
% yw=[];
% 
% %************************************ Interior particles generation ******************************************
% 
% %%Uncomment the following lines to get a random particle distribution %
% 
% % nint=input('Enter an even number of interior particles: ');
% % 
% % a=linspace(r0,(R-4*r0),nint);
% % weight=linspace(r0,(R-4*r0),nint);
% % 
% % rho = randsample(a, nint, true, weight);
% % rho=rho';
% % %rho = 0.1*R + randperm().(R-0.1*R);
% % %rho = (R-4*r0).*rand(nint,1);
% % %rho = flipud(rho1);
% % theta=(2*pi).*rand(nint,1);
% % 
% % for i=1:1:nint
% %     p(i)=xi+rho(i)*cos(theta(i)); 
% %     q(i)=yi+rho(i)*sin(theta(i));
% % end
% % 
% % %interior_particles = [[0:nint-1]',p',q',ones(nint,1)*20000]; 
% % 
% % N_int=nint;
% 
% %************************************ Interior particles generation *******************************
% 
% %%Uncomment the following lines to get a ordered particle distribution %
% 
% for r=r0:2*r0:R-r0
% 	delTheta = 2*asin(r0/r);
% 	for theta = 0*pi:delTheta:2*pi
% 	    p = [p;xi + r*cos(theta)];
%     	q = [q;yi + r*sin(theta)];
%     end
% end
% 
% N_int_temp = size(p,1);
% 
% p=p';
% q=q';
% 
% if mod(N_int_temp,2)==1 % ensuring no odd number of particles since we have to give pairwise random velocities later.
%     p(N_int_temp)=[];
%     q(N_int_temp)=[];
% end
% 
% N_int=size(p',1)
% %************************************ moving wall formation ************************************
% 
% delR = r0 * sqrt(3);
% for layer = 0:3;
% 	r = R + layer*delR;
% 	delTheta = 2*asin(r0/r);
% 	ini_angle = 0.0 + layer*(delTheta/2);
% 	final_angle = 2*pi + layer*(delTheta/2);
% 	for theta = ini_angle : delTheta : final_angle
% 	    xw = [xw;xi + r*cos(theta)];
%     	yw = [yw;yi + r*sin(theta)];
% 	end
%  	if (layer == 0)
%  		xw1 = size(xw,1)
% 	end
% 	if (layer == 1)
% 		xw2 = size(xw,1) - xw1
% 	end
% 	if (layer == 2)
% 		xw3 = size(xw,1) - xw2 - xw1
%     end
%     if (layer == 3)
%  		xw4 = size(xw,1) -xw3  -xw2 - xw1
%     end
% end
% 
% xwm=size(xw,1)
% 
% %****************************** positions arrays generation ******************************
% 
% interior_particles = [[0:N_int-1]', p', q', ones(N_int,1)*20000];
% 
% particles_belt = [[N_int:N_int+xwm-1]', xw, yw, ones(xwm,1)*(1e+25)];
% 
% Total_particles = vertcat(interior_particles,particles_belt);
% 
% dlmwrite('init_positions.txt', Total_particles, 'delimiter','\t','precision',7)
% 
% N_total = size(Total_particles,1)
% 
% 
% %************************* Velocity array generation ***********************%
% thetacap = (2).*rand(N_int,1)-1;
% epsilon = 10e-8;
% 
% list = randperm(N_int,N_int)';
% p1=(1).*rand(1,1);
% p2=(1).*rand(1,1);
% 
% for i=1:2:N_int-1
%      Vx(list(i))=p1*epsilon*(thetacap(list(i))); 
%      Vx(list(i+1))=-Vx(list(i)); % pairwise equal and opposite random velocities
%      Vy(list(i))=p2*epsilon*(thetacap(list(i)));
%      Vy(list(i+1))=-Vy(list(i));
% end
% 
% for j=N_int+1:1:N_total
%      Vx(j)=0; 
%      Vy(j)=0;
% end
% 
% vel_particles = [[0:N_total-1]',Vx',Vy'];
% 
% dlmwrite('init_velocities.txt',vel_particles, 'delimiter', '\t', 'precision', 7)
% 
% %******************************* plots *****************************%
%  
% hold on
% scatter(p,q)
% scatter(x,y)
% scatter(xw,yw)
% %quiver(xw,yw,Vx,Vy)
% hold off
% axis equal

%%%%%%%%%%%% Part pertaining to rectangular domain %%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

rng 'shuffle'

delete init_positions.txt;
delete init_velocities.txt;
delete beta.txt;

xi=0.4;%input('x co-ord of the cavity:');
yi=0.4;%input('y co-ord of the cavity:');
r0=0.0025;%input('radius of the particle:');
d0=2*r0;
halfwidth=0.2;%input('half width of the cavity:');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% external particle generation %%%%%%%%

% 1st bottom wall %
x1=[]; %% d= x co-ord of the wall under the belt particle
y1=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

x1=xi-halfwidth-d0:d0:xi+halfwidth+d0;
temp=size(x1);
y1=ones(temp)*(yi-halfwidth-d0);

%%

% 1st right wall %
x2=[]; %% d= x co-ord of the wall under the belt particle
y2=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

y2=yi-halfwidth-d0:d0:yi+halfwidth+d0;
temp=size(y2);
x2=ones(temp)*(xi+halfwidth+d0);

%%

% 1st top wall %
x3=[]; %% d= x co-ord of the wall under the belt particle
y3=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

x3=xi+halfwidth+d0:-d0:xi-halfwidth-d0;
temp=size(x3);
y3=ones(temp)*(yi+halfwidth+d0);

%%

% 1st left wall %
x4=[]; %% d= x co-ord of the wall under the belt particle
y4=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

y4=yi+halfwidth+d0:-d0:yi-halfwidth-d0;
temp=size(y4);
x4=ones(temp)*(xi-halfwidth-d0);

%%


% 2nd bottom wall %
x5=[]; %% d= x co-ord of the wall under the belt particle
y5=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

x5=xi-halfwidth-2*d0:d0:xi+halfwidth+2*d0;
temp=size(x5);
y5=ones(temp)*(yi-halfwidth-2*d0);

%%

% 2nd right wall %
x6=[]; %% d= x co-ord of the wall under the belt particle
y6=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

y6=yi-halfwidth-2*d0:d0:yi+halfwidth+2*d0;
temp=size(y6);
x6=ones(temp)*(xi+halfwidth+2*d0);

%%

% 2nd top wall %
x7=[]; %% d= x co-ord of the wall under the belt particle
y7=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

x7=xi+halfwidth+2*d0:-d0:xi-halfwidth-2*d0;
temp=size(x7);
y7=ones(temp)*(yi+halfwidth+2*d0);

%%

% 2nd left wall %
x8=[]; %% d= x co-ord of the wall under the belt particle
y8=[]; %% e= y co-ord of the wall under the belt particle
temp=[];

y8=yi+halfwidth+2*d0:-d0:yi-halfwidth-2*d0;
temp=size(y8);
x8=ones(temp)*(xi-halfwidth-2*d0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xwall=horzcat(x1,x2,x3,x4,x5,x6,x7,x8);
ywall=horzcat(y1,y2,y3,y4,y5,y6,y7,y8);

xp=[];
yp=[];

%%%%%%% internal particle generation %%% for an ordered grid
k=1; f = 1.03;

for y=yi-halfwidth:f*d0:yi+halfwidth
    for x=xi-halfwidth:f*d0:xi+halfwidth
     
        xp(k)=x;    
        yp(k)=y;
        k=k+1;
    
    end
end

nint_temp = length(xp);

if mod(nint_temp,2)==1 % ensuring no odd number of particles since we have to give pairwise random velocities later.
    xp(nint_temp)=[];
    yp(nint_temp)=[];
end

%%%%%%% internal particle generation %%% for random grid %%%%%

% n=input('Even number of interior particles in the domain:');
% 
% xp = xi-halfwidth + (2*halfwidth).*rand(n,1); % for random locations of the interior particles
% yp = yi-halfwidth + (2*halfwidth).*rand(n,1);
% 
% xp=xp';
% yp=yp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nint=length(xp)
nwall=size(xwall,2)
n_total = nwall+nint

xp=xp';
yp=yp';



%%%% beta %%%%

ratio = 0.1;
beta1 = 1;
beta2 = beta1*ratio;

beta = beta1*(ones(n_total,1));

kk=round(nint/2);

beta(randsample(nint,kk)) = beta2;

%%%%%%%%%%%%%%%%%%%%

interior_particles = [[0:nint-1]',xp,yp,ones(nint,1)*20000];

wall_particles = [[nint:nint+nwall-1]',xwall',ywall',ones(nwall,1)*(1e+25)]; 

Total_particles = vertcat(interior_particles,wall_particles);

dlmwrite('init_positions.txt', Total_particles, 'delimiter','\t','precision',8)

temp_beta = horzcat((1:n_total)',beta);
dlmwrite('beta.txt', temp_beta,'delimiter','\t','precision',8)

thetacap = (2).*rand(nint,1)-1;
epsilon = 10e-8;

list=randi(nint,nint,1);
p1=(1).*rand(1,1);
p2=(1).*rand(1,1);

for i=1:2:nint/2-1
     Vx(list(i))=p1*epsilon*(thetacap(list(i))); 
     Vx(list(i+1))=-Vx(list(i)); % pairwise equal and opposite random velocities
     Vy(list(i))=p2*epsilon*(thetacap(list(i)));
     Vy(list(i+1))=-Vy(list(i));
end

for j=nint/2+1:1:n_total
     Vx(j)=0; 
     Vy(j)=0;
end

vel_particles = [[0:n_total-1]',Vx',Vy'];

dlmwrite('init_velocities.txt',vel_particles, 'delimiter', '\t', 'precision', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('visible', 'off');

hold on
scatter(xwall,ywall,'filled')
scatter(xp,yp,[],beta(1:nint),'filled')
hold off
axis equal
box on
colormap(gray(256));
colorbar;

print -depsc scatterplot.eps
close(f1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
