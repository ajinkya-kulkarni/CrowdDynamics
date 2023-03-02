
clc;
clear all;

%%% Radial and Azimultal velocity calculations %%%

nint=6120;

j=1;

for filenum=0:1:400
filenum
cd particle_info
if (filenum<10)
    fname = strcat('part_data_000',int2str(filenum));
elseif (filenum<100)
    fname = strcat('part_data_00',int2str(filenum));
elseif (filenum<1000)
    fname = strcat('part_data_0',int2str(filenum));
end
    
A = load(fname,'-ascii');
cd ..

cd particle_vel_info
if (filenum<10)
    fname = strcat('part_vel_data_000',int2str(filenum));
elseif (filenum<100)
    fname = strcat('part_vel_data_00',int2str(filenum));
elseif (filenum<1000)
    fname = strcat('part_vel_data_0',int2str(filenum));
end

B = load(fname,'-ascii');
cd ..

C=[[1:nint]', A(1:nint,2:3), B(1:nint,2:3)];

x = C(:,2);
y = C(:,3);
vx = C(:,4);
vy = C(:,5);

for i=1:1:nint
    theta(i)=atan2((y(i)-0.4),(x(i)-0.4));
    vradial_single(i) = vy(i)*sin(theta(i)) + vx(i)*cos(theta(i));
    vtheta_single(i) = vy(i)*cos(theta(i)) - vx(i)*sin(theta(i));
    r(i) = sqrt((x(i)-0.4)^2 + (y(i) - 0.4)^2);
end

vradial=(mean(vradial_single));
vtheta=(mean(vtheta_single));

speed(j)=sqrt(vradial*vradial + vtheta*vtheta);

j=j+1;

end


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %delete Binder_Cumulants.txt

dlmwrite('speed_final.txt',speed','precision', 16);

% dlmwrite('Binder_Cumulants.txt',matrix,'delimiter', '\t','precision', 16);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
 
%clc
%clear all

%f1 = figure('visible', 'off');
%a=dlmread('speed_final.txt');

%speed = a(:,1);

%plot((1:1:length(speed)),speed,'-o')
%grid on
%box on
%hold off
%xlabel('Time')
%ylabel('$V_{avg}$','interpreter','latex')
%title('Binder Cumulants vs C_{v}')

%print -depsc speed.eps
%close(f1)
% 
%exit;s
