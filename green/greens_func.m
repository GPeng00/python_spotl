function greens_func(deps,lame)
for mm=1:length(deps)
    dep=deps(mm)
%!\rm green.bou.*
% DISTS 
ds1=.025:.01:.995;
ds2=1.05:.1:9.95;
ds3=10.25:.5:89.75;
ds4=90.5:1:179.5;
ds=vertcat(ds1', ds2', ds3', ds4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SOLUTIONS FROM MALVERN 1969 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% everything in r, theta and d
a=6371136; %earth radius
g=9.81;
nu=lame.nu; %
G=lame.G; % 
P=1*g; % load in Newtons
E=lame.E;
depth=dep*1000;
for z=depth % in meters
x=ds(:,1) %fars(:,1)
for j=1:length(x)
    delta=x(j)*pi/180;
    r=111325*x(j); % a*delta in radians
    R=sqrt(r^2+z^2);
    %%%%%%%%% STRESS TENSOR
    T(1,1)=P/(2*pi*R^2)*((-3*r^2*z)/(R^3)+(R*(1-2*nu))/(R+z));  % Trr
    T(2,2)=P*(1-2*nu)/(2*pi*R^2)*(z/R-R/(R+z));  % Ttt
    T(3,3)=-3*P*z^3/(2*pi*R^5);  % Tzz
    T(1,3)=-3*P*r*z^2/(2*pi*R^5); % Trz
    %%%%%%%%% Stress to Strain %%%%%%%%%%%
    e(1,1)=1/E*(T(1,1)-nu*(T(2,2)+T(3,3))); % e_rr
    e(2,2)=1/E*(T(2,2)-nu*(T(1,1)+T(3,3))); % e_tt
    e(3,3)=1/E*(T(3,3)-nu*(T(1,1)+T(2,2))); % e_zz
    e(1,3)=T(1,3)/(2*G); % e_rz
    bous(j,1)=e(1,1)*10^12*(a*delta)^2; % e_rr 
    bous(j,2)=e(2,2)*10^12*(a*delta)^2; % e_tt  
    bous(j,3)=e(3,3)*10^12*(a*delta)^2; % e_zz
    bous(j,4)=e(1,3)*10^12*(a*delta)^2; % e_rz
end
bous=horzcat(ds(:,1), bous);

%%%%%%%%%%%%%%% SPOTL Normalization %%%%%%%%%%%%%%%%%%%%%%%%

K=10^12;
for i=1:length(bous)
    Gtp(i,1)=bous(i,2)/K*((2*sin(bous(i,1)*pi/180*.5))/(bous(i,1)*pi/180))^2; % e_rr 
    Gtp(i,2)=bous(i,3)/K*((2*sin(bous(i,1)*pi/180*.5))/(bous(i,1)*pi/180))^2; % e_tt 
    Gtp(i,3)=bous(i,4)/K*((2*sin(bous(i,1)*pi/180*.5))/(bous(i,1)*pi/180))^2; % e_zz 
    Gtp(i,4)=bous(i,5)/K*((2*sin(bous(i,1)*pi/180*.5))/(bous(i,1)*pi/180))^2; % e_rz
end

Delta=ds;
del=vertcat(.01*ones(98,1),.1*ones(90,1),.5*ones(160,1), 1*ones(90,1));

bGF=[];
for i = 1:size(ds,1)
    fac=log((sin(Delta(i,1)*pi/180/2)*cos(del(i,1)*pi/180/4)+cos(Delta(i,1)*pi/180/2)*sin(del(i,1)*pi/180/4))/(sin(Delta(i,1)*pi/180/2)*cos(del(i,1)*pi/180/4)-cos(Delta(i,1)*pi/180/2)*sin(del(i,1)*pi/180/4))); 
    bGF(i,1)=Gtp(i,1)*fac; % e_rr --> e_tt
    bGF(i,2)=Gtp(i,2)*fac; % e_tt --> e_ll
    bGF(i,3)=Gtp(i,3)*fac; % e_zz --> e_zz
    bGF(i,4)=Gtp(i,4)*fac; % e_rz --> e_tz
end

bGF=horzcat(Delta, bGF); % (2:e_tt) (3:e_ll) (4:e_zz) (5:e_tz) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% FORMAT OUT %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load green % this is the file that you use the potentials, displacements, etc. from
%new=horzcat(green(:,1:4), vertcat(bGF(1:188,2:3),green(189:length(ds),5:6)) ,bGF(:,5), bGF(:,7), green(:,7));
new=horzcat(green(:,1:4),...
    vertcat(bGF(1:188,2:3), zeros(length(ds)-188,2)),...
    green(:,7),...
    vertcat(bGF(1:188,4:5), zeros(length(ds)-188,2)));

file=['green.bou.' num2str(dep)]
fid = fopen(file,'w');
fprintf(fid,'%s\n',['HALF SPACE GREENS FUNCTION ' num2str(dep) ' KM DEPTH']);
fprintf(fid,'%i %i %i %i %5.4f %5.4f %5.4f %c\n',9,1,4,length(ds1),ds1(1), ds1(end), ds1(2)-ds1(1),'F');
for i=1:length(ds1)
    for j=1:size(new,2)-1
        fprintf(fid,' %13.6e',new(i,j));
    end
    fprintf(fid,' %13.6e\n',new(i,j+1));
end
fprintf(fid,'%i %i %i %i %5.4f %5.4f %5.4f %c\n',9,2,4,length(ds2),ds2(1), ds2(end), ds2(2)-ds2(1),'F');
for i=1+length(ds1):length(ds1)+length(ds2)
    for j=1:size(new,2)-1
        fprintf(fid,' %13.6e',new(i,j));
    end
    fprintf(fid,' %13.6e\n',new(i,j+1));
end
fprintf(fid,'%i %i %i %i %5.4f %5.4f %5.4f %c\n',9,3,4,length(ds3),ds3(1), ds3(end), ds3(2)-ds3(1),'C');
for i=1+length(ds1)+length(ds2):length(ds1)+length(ds2)+length(ds3)
    for j=1:size(new,2)-1
        fprintf(fid,' %13.6e',new(i,j));
    end
    fprintf(fid,' %13.6e\n',new(i,j+1));
end
fprintf(fid,'%i %i %i %i %5.4f %5.4f %5.4f %c\n',9,4,4,length(ds4),ds4(1), ds4(end), ds4(2)-ds4(1),'C');
for i=1+length(ds1)+length(ds2)+length(ds3):length(ds1)+length(ds2)+length(ds3)+length(ds4)
    for j=1:size(new,2)-1
        fprintf(fid,' %13.6e',new(i,j));
    end
    fprintf(fid,' %13.6e\n',new(i,j+1));
end
fclose(fid);

% command = sprintf('%s %s %d %d', 'makegreen.scr', 'in', abs(z/1000), size(new,2));
% unix(command)

figure
subplot(2,1,1)
plot(ds, green(:,5),'g.')
hold on
plot(ds, new(:,5),'b.')
xlim([0 10])
title('e_r_r')
subplot(2,1,2)
plot(ds, green(:,6),'g.')
hold on
plot(ds, new(:,6),'b.')
xlim([0 10])
title('e_t_t')
end
end
