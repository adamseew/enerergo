
% Run default energy-aware ergodic search experiments in "Energy-aware 
% ergodic search: Continuous long-term exploration for multiagent 
% systems"


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
epsilon=0.05; % tollerance interval
k=9; % number of freq. in Fourier transform
[K(1,:,:) K(2,:,:)]=ndgrid(0:1:k,0:1:k); % set of indices
min_b=.1; % minimum battery constraint


%% experiment specific data

RECHARGE_GAIN=7;

uav_x0=[[.1;.3] ...
        [.9;.7] ...
        [.1;.7] ...
        [.9;.3]]; % charging station at start
                  % uavs fly two-by-two, always landing at each others
                  % bases;
uav_xf=[uav_x0(:,2) ...
        uav_x0(:,1) ...
        uav_x0(:,4) ...
        uav_x0(:,3)];

Z_=100; % used in main loop
uav_x0_b3=Z_*ones(4,1); % batteries SoC, initial

N_STEP=2;
LOOP_UNTIL_OPT=1;

results=[];
count=1;
uav_no=1;

while 1
    x0=uav_x0(:,uav_no);
    xf=uav_xf(:,uav_no);
    uav_x0(:,uav_no)=xf;
    uav_xf(:,uav_no)=x0;
    
    results=[results;{count,uav_no,x0,xf,nan,nan,nan}];

    switch uav_no
        case {1,3}
            load('gauss_1689953587.mat'); % loads alpha, Mu, Sigma from 
                                          % configuration
            
        case {4,2}
            load('gauss_1689949481.mat');
    end

    if count>4 % changes initial alpha guesses from previous iterations
        alpha=results{end-4,6};
    end

    Z_=uav_x0_b3(uav_no);
    N=500; % horizon

    run('../enaware');
    
    Z_=x0_b(3);
    
    if Z_<min_b % finish once any of the UAVs' batteries is depleted
        break; 
    end

    results(end,5)={Z_};
    results(end,6)={debug.alpha};
    results(end,7)={N};

    if exist('ALPHA_CONST_','var')
        ALPHA_CONST=ALPHA_CONST_;
    end
    
    [~,ALPHA_CONST_]=sort(debug.alpha);
    uav_x0_b3(uav_no)=Z_;

    save(strcat('data/part_no_',string(count),'.mat')); % saving partial 
                                                        % results
    savefig(figure(count),strcat('data/part_no_',string(count),'.fig'));

    clear('alpha','Mu','Sigma');

    uav_no=uav_no+1; % iterate the 4 UAVs
    if uav_no>4 
        uav_x0_b3=uav_x0_b3+RECHARGE_GAIN*ones(4,1);
        uav_no=1;
    end

    count=count+1;
end

save(strcat('data/results.mat'),'results'); % saving results summary table

