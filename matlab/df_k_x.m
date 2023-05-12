function df_k_x_val = df_k_x(k,x,args) 
%DF_K_X gradient of the f_k_x function to define the control action
%   Fuction which returns the coefficients of the gradient with respect to
%   x of the f_k_x function to define the control action
% Inputs:
%   k   : current index from the set of indices K
%   x   : point as a row vector (works with CasADi variables)
%   args: additional arguments. You must pass at least D (dimension), and 
%         L (period)
% Outputs:
%   df_k_x_val: 2xlength(k) gradient coefficients

    if args.D~=2
        error("myComponent:notImplemented",strcat("Error. \nlinear df_",...
              "k_x not yet implemented for dimensions other than 2 (di",...
              "mension is %d)"), args.D);
    end

    % warning('Only the real part of df_k_x is currently implemented');

    % ?        .*                            .*
    df_k_x_val=-2*pi*args.L^(-args.D-1)*...
        [...
         k(1,:).*sin(2*pi*x(1)*k(1,:)/args.L).*cos(2*pi*x(2)*k(2,:)/args.L);...
         k(2,:).*cos(2*pi*x(1)*k(1,:)/args.L).*sin(2*pi*x(2)*k(2,:)/args.L) ...
        ];
end

