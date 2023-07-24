function [dx, Voc] = thevenin(x, I, args)
%THEVENIN implements the battery Thevenin ECM model from 
%https://www.doi.org/10.1109/TCST.2016.2542115
% Inputs:
%   x   : battery state. A 3x1 vector with R v1, v2 and battery SoC
%   I   : load in amperes
%         pass 2 rows times N columns (works with CasADi variables)
%   args: additional arguments. You must pass at least V (battery 
%         voltage), Rs (R in the ECM), R1 and C1 (first RC element 
%         in the ECM), R2 and C2, Q (battery capacity)
% Outputs:
%   dx:   derived battery state
%   Voc:  battery Voc

    if size(x,1)~=3
        error("myComponent:invalidArgument",strcat("Error. \nThe state",...
              " must be a row vector of three values. v1, v2, and z (s",...
              "tate of charge)"));
    end
    
    Voc=args.V+x(1)+x(2)+I*args.Rs;
    dx=[-1/(args.R1*args.C1) 0 0; 0 -1/(args.R2*args.C2) 0; 0 0 0]*x+...
        [1/args.C1;1/args.C2;-1/args.Q]*I;
end
