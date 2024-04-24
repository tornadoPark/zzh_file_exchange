clc
clear all
close all

global rho_l
global lambda_l

xMin = 0;
xMax = 0.04;

rMin = 0.001;%15e-3;
rMax = 0.004;%19e-3;
rMaxtotal = 0.006;%19e-3;
ncellx = 60;
ncellr = 20;
ncellrtotal = 30;

dt = 0.001;

% uLid = [1,0];
pTop = 100;
TGradBot = [0,1.4e4];
TTop = 700;
h=100;
T_inf=700;

Ttr = 686;
T_vapor=Ttr; 
M=0.023;
Ru=8.314;  
Q_evap = 500;% 总传热量 （单位 W）
hfg_vapor = exp(-57.566+0.18157*T_vapor-2.2885e-4*power(T_vapor,2)+1.5614e-7*power(T_vapor,3)-5.5058e-11*power(T_vapor,4)+7.8615e-15*power(T_vapor,5))*0.001;
l_cond = 0.5*xMax;
uLid = [0,-2e-3];
Q_evap_sonic=0.1;
delta_T_vapor=0.00001;
n_points_x = ncellx;
mflag=zeros(1,n_points_x);
mac=0.2;
alpha_evapcond = (2*mac/(2-mac)).*sqrt(M/2/3.1415/Ru);  %%%%  常系数
Q_in = norm(TGradBot*16.3*pi*2*rMaxtotal*0.25*xMax);


h0 = figure;
h1 = figure;
h2 = figure;
% h3 = figure;
%% initialization
generateMesh
generateMeshtotal

% determine the B.C. type of each face:
bType = -1*ones(ncellr,ncellx,4);
bType(:,1,1) = 0;	% 0 means left wall
bType(:,ncellx,3) = 0;	% 0 means right wall
bType(1,:,2) = 0;	% 0 means down wall
%bType(1,:,4) = 4;	 4 means up boundary
bType(ncellr,1:ncellx*0.3,4) = 41;	%   evap-outlet
bType(ncellr,ncellx*0.3+1:ncellx*0.5,4) = 42;   %  adiabatic
bType(ncellr,ncellx*0.5 + 1:ncellx ,4) = 43;	%  cond-inlet

disp('B.C. has been determined*********')


% determine the B.C. type of each face:
bTypetotal = -1*ones(ncellrtotal,ncellx,4);
bTypetotal(:,1,1) = 0;	% 0 means left wall
bTypetotal(:,ncellx,3) = 0;	% 0 means right wall
bTypetotal(1,:,2) = 0;	% 0 means down wall
%bType(1,:,4) = 4;	 4 means up boundary
bTypetotal(ncellrtotal,1:ncellx*0.3,4) = 41;	%   evap-outlet
bTypetotal(ncellrtotal,ncellx*0.3+1:ncellx*0.5,4) = 42;   %  adiabatic
bTypetotal(ncellrtotal,ncellx*0.5 + 1:ncellx ,4) = 43;	%  cond-inlet

disp('B.C.total has been determined*********')

u = zeros(ncellr,ncellx,2);
p = zeros(ncellr,ncellx);
T = TTop+30+zeros(ncellrtotal,ncellx);

disp('initial fields have been determined*********')


nu = 2.5e-7; % kinematic viscosity = dynamic viscosity / rho
kporous = 1e-7;
epsilon = 0.67; %porosity 
% lambda = 1*ones(ncellrtotal,ncellx);
% rhocp = 1*ones(ncellrtotal,ncellx);
Rhocp_sh=8000*1260;
alpha = getVolumeValueOfAlpha(T,epsilon,ncellr,ncellrtotal,ncellx); % heat diffusivity = k/rho/cp, alpha of water approx 1e-7
% alpha = lambda./rhocp;
alpha_f = getFaceValueOfAlpha(T,alpha,ioffset,joffset,bTypetotal,ncellr,ncellrtotal,ncellx,epsilon);

disp('properties are assigned *******')





disp('initialization is done********')
relaxP = 0.8; % Jasak suggest 0.2 on thesis page 149
relaxU = 0.5;
relaxT = 0.95;

%% calculation
% variables initialization
% F = zeros(ncellr,ncellx,4);
% VbyA = zeros(ncellr,ncellx);
% VbyA_f = zeros(ncellr,ncellx,4);
% HbyA_f = zeros(ncellr,ncellx,4);

ux= u(:,:,1);
uxOld = ux;
res_time = [];

for time = 1:100000
    uLast = u;
    TLast = T;

    subIterMax = 3000;
    for iter = 1:subIterMax

        gradP = grad(ncellx,ncellr,dx,dr,p,pTop);
        solveVelocity
        solvePressure

        gradPOld = gradP;
        gradP = grad(ncellx,ncellr,dx,dr,p,pTop);

        correctVelocity
        utotal=[zeros(ncellrtotal-ncellr,ncellx,2);u];  
        alpha = getVolumeValueOfAlpha(T,epsilon,ncellr,ncellrtotal,ncellx); % heat diffusivity = k/rho/cp, alpha of water approx 1e-7
        alpha_f = getFaceValueOfAlpha(T,alpha,ioffset,joffset,bTypetotal,ncellr,ncellrtotal,ncellx,epsilon);
        solveTemperature

        disp('here')

        
        q_lv = -lambda_l(1,:).*(T(2,:)-T(1,:))./dr;%蒸发应为负值
        T_surf = T(1,:) - (-T(1,:)+T(2,:))/2; % 外推半格温度（即气液界面温度）
        transient_cond_vapor_stage2
        TTop = T_vapor;% 指定温度方程边界条件为第一类边界条件
        % 指定温度方程边界条件为第二类边界条件
        % 指定速度方程边界条件在geruLid函数内
        uLidy = Q_evap/(2*rMin*pi*l_cond)/rho_l(1,ncellx)/hfg_vapor;
        uLid = [0,uLidy];

        p_vapor;% 指定速度方程压力边界值
        
        %% residual error and picture
        maxU = max(max((u(:,:,1))));
        if iter == 2
            res0=abs(maxU - maxUOld);
        end
        if iter > 2 && iter <= subIterMax
            res = abs(maxU - maxUOld)/maxUOld;
            %             disp(res)
            if res <= 1e-4
                disp(strcat('time=', string(time*dt),'s','     converge (iter=', string(iter),')'))
%                 solveTemperature
%                 alpha = getVolumeValueOfAlpha(T);
%                 alpha_f = getFaceValueOfAlpha(alpha,ioffset,joffset,bType);
                break
            end
            if iter == subIterMax && res >1e-4
                disp(strcat('time=', string(time*dt),'s','     Unconverge (res=', string(res),')'))
%                 solveTemperature
%                 alpha = getVolumeValueOfAlpha(T); 
%                 alpha_f = getFaceValueOfAlpha(alpha,ioffset,joffset,bType);
            end
        end
        maxUOld = maxU;
    end

    if mod(time,5)==0

        figure(h0);
        cla;
        paint(utotal(:,:,1),dx,dr)
        title(strcat(string(time*dt),' s'))
        colorbar
        drawnow

        figure(h1);
        cla;
        paint(utotal(:,:,2),dx,dr)
% surf(u(:,:,2),dx,dr)
        title(strcat(string(time*dt),' s'))
        colorbar
        drawnow
        
        figure(h2);
        cla;
        paint(T,dx,dr)
% surf(u(:,:,2),dx,dr)
        title(strcat(string(time*dt),' s'))
        colorbar
        drawnow
        
        max(max(u(:,:,2)))

        res_time(end+1) = norm(ux - uxOld,"fro");
        if res_time(end) <= 1e-5
            disp(strcat('Steady'))
            break
        end
        uxOld = ux;
    end

    if mod(time,500) == 0
        name = [num2str(time),'.mat'];
        save(name,"p","u","T");
    end
    
end
