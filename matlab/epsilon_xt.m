function epsilon_xt_val = epsilon_xt(xt,args,debug)
%EPSILON_XT_VAL goal function of the ergodic control, i.e., function to
%minimize
%   Function which defines a metric in the spectral domain, decomposing the
%   trajectory and the spatial distribution phi_x
% Inputs:
%   xt   : trajectory, as a row vector, e.g., in 2D with N points, you
%          will pass 2 rows times N columns (works with CasADi variables)
%   args : additional arguments. You must pass at least K (indices), and D 
%          (dimension). There are additional parameters that might be 
%          needed for the dependencies (c_k_xt, phi_k, and f_k_x)
%   debug: input empty. Returns detailed debug trajectories for all the 
%          components, i.e., c_k_xt, phi_k and Lambda_k for each value in K
% Outputs:
%   epsilon_xt_val: 1x1 metric (in R1/C1, depending on f_k_x) in spectral
%                   domain

    if size(args.K,1)~=args.D
        error("myComponent:invalidArgument",strcat("Error. \nEach set ",...
              "of indices in K has to have the same size as the dimens",...
              "ion (D=%d)"), args.D);
    end

    epsilon_xt_val=0;

    % weights assigning more importance to low frequency components
    Lambda_k=zeros(length(args.K));
    
    % re-setting debug variables
    debug.c_k_xt=nan(length(args.K),1);
    debug.phi_k=nan(length(args.K),1);
    debug.Lambda_k=nan(length(args.K),1);

    % iterating each index in K
    for k=1:length(args.K)
        
        % fourier series coefficients (FSC) along the trajectory x
        c_k_xt_val=c_k_xt(args.K(:,:,k),xt,args);
        phi_k_val=phi_k(args.K(:,:,k),args); % FSC of the spatial 
                                             % distribution (in the model)

        debug.c_k_xt=[debug.c_k_xt;c_k_xt_val];
        debug.phi_k=[debug.phi_k;phi_k_val];
        
        % building the Lambda_k diagonal matrix
        for j=1:length(args.K)
            Lambda_k(j,j)=(sum(args.K(:,j,k).^2)+1).^(-(args.D+1)/2);
            
            debug.Lambda_k=[debug.Lambda_k;Lambda_k(j,j)];
        end

        %                         '*        *
        epsilon_xt_val=epsilon_xt_val+...
            (c_k_xt_val-phi_k_val)'*Lambda_k*(c_k_xt_val-phi_k_val);
    end

    epsilon_xt_val=epsilon_xt_val/2;
end

