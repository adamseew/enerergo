function c_k_xt_val = c_k_xt(k,xt,args)
%C_K_XT fourier series coefficients along the trajectory x
%   Function which returns the fourier series coefficients along the 
%   trajectory x of the duration N in args (horizon)
% Inputs:
%   k   : current index from the set of indices K
%   xt  : trajectory as a row vector, e.g., in 2D with N points, you will
%         pass 2 rows times N columns (works with CasADi variables)
%   args: additional arguments. You must pass at least N (horizon)
% Outputs:
%   c_k_xt_val: length(k)x1 fourier series coefficients along the
%               trajectory x

    if args.N<1
        error("myComponent:invalidArgument",strcat("Error. \nInvalid v",...
              "alue for the horizon N"));
    end

    c_k_xt_val=0;
    
    % discretized version
    for s=1:floor(args.N)
        c_k_xt_val=c_k_xt_val+f_k_x(k,xt(:,s),args);
    end

    c_k_xt_val=c_k_xt_val/floor(args.N);
end

