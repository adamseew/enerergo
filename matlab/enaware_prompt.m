
% Energy-aware ergodic search 1/2
% i.e., an ergodic controller that uses the
% controller from econtrol.m to derive a control action that respects the
% battery constraints while maximazing the coverage


answer=questdlg('Would you like to clean the environment?',...
    'Clean',...
	'Yes','No','No');
switch answer
case 'Yes'
    clear;
end
clear("answer");


%% definitions

D=2; % dimension, e.g., 2D, 3D, etc.
L=2; % period
epsilon=0.05; % tollerance interval
approx_f=1; % approximation factor, i.e., 1 no approximation. Use powers
            % of 10

if D~=2
    error("myComponent:notImplemented",strcat("Error. \nlinear transf",...
          "ormation matrices Am, indeces K are not yet implemented for",...
          " dimensions other than 2 (dimension is %d)"), D);
end

disp("starting gauss app")
gauss_app=gauss("on"); % starting Guassian mixture model designer app
try
    waitfor(gauss_app,"closed","yes"); % wait for Guassian mixture model 
                                       % designer app to terminate
catch
    delete(gauss_app.ErgodiccontrollerdesignerUIFigure)
    error("Data from gauss.mlapp are required")
end
Mu=gauss_app.Mu;
Sigma=gauss_app.Sigma;
alpha=gauss_app.alpha;
k=gauss_app.k; % number of freq. in Fourier transform
[K(1,:,:) K(2,:,:)]=ndgrid(0:1:k,0:1:k); % set of indices
N=gauss_app.N; % horizon
x0=gauss_app.x0; % initial guesses
xf=gauss_app.xf; % desired final point, e.g., recharge station

delete(gauss_app.ErgodiccontrollerdesignerUIFigure)
if length(alpha)>1
    disp(strcat("there are ",string(length(alpha)),...
        " gaussians in the model"))
else
    disp("there is one gaussian in the model")
end
clear("gauss_app")

LOOP_UNTIL_OPT=1; % set to 1 so that the controller is refined until an
                  % optimal configuration is found. Set to 0 otherwise
N_STEP=3;

% call the main enaware script

min_b=.05; % minimum battery constraint

enaware

