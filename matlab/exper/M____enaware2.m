
% Run default energy-aware ergodic search experiments in "Energy-Aware 
% Ergodic Search: Continuous Exploration for Multi-agent Systems with 
% Battery Constraints"


answer=questdlg('Would you like to clean the environment?',...
    'Clean',...
	'Yes','No','No');
switch answer
case 'Yes'
    clear;
end
clear("answer");

CASADI_PATH='~/casadi';
addpath(CASADI_PATH);


%% definitions

D=2; % dimension, e.g., 2D, 3D, etc.
L=2; % period
epsilon=.5; % tollerance interval
N=1000; % horizon
k=9; % number of freq. in Fourier transform
[K(1,:,:) K(2,:,:)]=ndgrid(0:1:k,0:1:k); % set of indices
min_b=.1; % minimum battery constraint


%% experiment specific data

RECHARGE_GAIN=21;

uav_x0=[.5;.5]; % charging station at start
                % the charging station is just one;
uav_xf=uav_x0;

Z_=100; % used in main loop
uav_x0_b3=Z_; % batteries SoC, initial

N_STEP=2;
LOOP_UNTIL_OPT=1;
ALPHA_FIX=1; % keep traditional ergodic search version
results=[];
count=1;
uav_no=1;
count=1;

while 1
    x0=uav_x0;
    xf=uav_xf;
    
    results=[results;{count,uav_no,x0,xf,nan,nan,nan}];

    load('gauss_1693923525.mat'); % loads alpha, Mu, Sigma from 
                                  % configuration

    Z_=uav_x0_b3;
   
    if count>=15
        break; 
    end

    run('../enaware');
    
    Z_=x0_b(3);

    results(end,5)={Z_};
    results(end,6)={debug.alpha};
    results(end,7)={N};

    if exist('ALPHA_CONST_','var')
        ALPHA_CONST=ALPHA_CONST_;
    end
    
    [~,ALPHA_CONST_]=sort(debug.alpha);
    uav_x0_b3(uav_no)=Z_;

    save(strcat('data2/part_no_',string(count),'.mat')); % saving partial 
                                                         % results
    savefig(figure(count),strcat('data2/part_no_',string(count),'.fig'));

    clear('alpha','Mu','Sigma');
    
    uav_x0=debug.x(:,end);
    uav_x0_b3=uav_x0_b3+RECHARGE_GAIN;

    count=count+1;
end

save(strcat('data2/results.mat'),'results'); % saving results summary table

