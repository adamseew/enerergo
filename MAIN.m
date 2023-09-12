function MAIN()
% Runs all the experiments in this folder
% Use econtrol for spatial distribution ergofic controller,
%     timeopti for  time-optimal ergodic controller with IPOPT,
%     eneraware     energy-aware ergodic controller with IPOPT,
%     calinon       ergodic controller from Calinon's work (src:
%                   https://doi.org/10.1007/978-3-030-23876-6_3)

% Depends on CasADi framework (src: https://web.casadi.org/get/)
%     Specify CasADi folder in CASADI_PATH

% Copyright (c) Adam Seewald, IA & GRAB Labs at Yale University
% Department of Mechanical Engineering and Materials Science 
% Distributed under CC BY-NC-SA licence
% Details: http://creativecommons.org/licenses/by-nc-sa/4.0/

    [indx,tf]=listdlg('PromptString',{'Modality',''},...
        'SelectionMode','single','ListSize',[300 160],'ListString',...
        {'pure ergodic control',...
         'time optimal ergodic control',...
         'energy-aware ergodic control',...
         'Calinon ergodic control'});

    CASADI_PATH='~/casadi';
    addpath(CASADI_PATH);
    cd matlab;

    switch indx
        case 1
            econtrol
        case 2
            timeopti
        case 3
            enaware_prompt
        case 4
            calinon
    end
end
