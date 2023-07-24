
% Generates plots im "Energy-aware ergodic search: Continuous long-term 
% exploration for multiagent systems"

% ! run after generating all the parts in data (run MAIN_eneaware.m first)


answer=questdlg('Would you like to clean the environment?',...
    'Clean',...
	'Yes','No','No');
switch answer
case 'Yes'
    clear;
end
clear("answer");


parts={'data/part_no_1.mat', 'data/part_no_2.mat',...
       'data/part_no_3.mat', 'data/part_no_4.mat',...
       'data/part_no_5.mat', 'data/part_no_6.mat',...
       'data/part_no_7.mat', 'data/part_no_8.mat',...
       'data/part_no_9.mat', 'data/part_no_10.mat',...
       'data/part_no_11.mat','data/part_no_12.mat',...
       'data/part_no_13.mat','data/part_no_14.mat',...
       'data/part_no_15.mat','data/part_no_16.mat',...
       'data/part_no_17.mat','data/part_no_18.mat',...
       'data/part_no_19.mat','data/part_no_20.mat',...
       'data/part_no_21.mat','data/part_no_22.mat'...
       'data/part_no_23.mat','data/part_no_24.mat',...
       'data/part_no_25.mat','data/part_no_26.mat',...
       'data/part_no_27.mat','data/part_no_28.mat'};
colors={'b','r'};
[X_,Y_]=meshgrid(0:.005:1);
plot_no=0;

WIDTH=122;
HEIGHT=WIDTH;
FIXED_SPACE=220;
SCREEN_HEIGHT=get(0,'screensize')*.75;
SCREEN_HEIGHT=SCREEN_HEIGHT(4)-FIXED_SPACE;

disp('started plotting')

for part_no=1:2:length(parts)

    fig=figure;
    set(fig,'color','w');
    ax1=gca;
    ax1.YGrid='on';
    ax1.Layer='top';
    ax1.GridLineStyle=':';
    ax1.GridAlpha=.5;
    set(ax1,'XTick',get(ax1,'YTick'));
    xticks([0 .2 .4 .6 .8 1]);
    yticks([0 .2 .4 .6 .8 1]);
    xticklabels({'0','0.6','1.2','1.8','2.4','3'});
    yticklabels({'0','0.6','1.2','1.8','2.4','3'}); % stretching to 3 m
    grid on;
    hold on; 

    xlim([0 1]);
    ylim([0 1]);
    axis square; % setting axes

    if part_no>=3
        for prv_part_no=1:part_no-1 % plotting previous trajectories
            load(parts{prv_part_no},'debug');

            plot(debug.x(1,:),debug.x(2,:),'color',[.35 .35 .35 .8]);
            % if prv_part_no>part_no-3 % if uncommented, plots previous
                                       % base stations in light gray
                % scttr=scatter(debug.x(1,1),debug.x(2,1),.85e3,...
                %     [.35 .35 .35],'filled'); % plotting previous base st.
                % clear('alpha');
                % alpha(scttr,.45);
            % end
        end
    end

    load(parts{part_no},'Sigma','Mu','debug');
    p=  debug.alpha(1)*mvnpdf([X_(:) Y_(:)],Mu(:,1)',Sigma(:,:,1))+...
        debug.alpha(2)*mvnpdf([X_(:) Y_(:)],Mu(:,2)',Sigma(:,:,2));
    load(parts{part_no+1},'Sigma','Mu','debug');
    p=p+debug.alpha(1)*mvnpdf([X_(:) Y_(:)],Mu(:,1)',Sigma(:,:,1))+...
        debug.alpha(2)*mvnpdf([X_(:) Y_(:)],Mu(:,2)',Sigma(:,:,2));
    p=reshape(p,size(X_));
    pclr=pcolor(X_,Y_,p);
    shading interp;
    set(pclr,'facealpha',0.65);
    colormap(flipud(bone(15))); % adds the prob. distributions (gauss. mix.
                                % in background)
    
    for part_no_=1:length(colors)
       
        load(parts{part_no+part_no_-1},'debug');
        
        % scttr=scatter(debug.x(1,1),debug.x(2,1),.85e3,...
        %     'filled',strcat(colors{part_no_},'o')); % plotting base st.
        % clear('alpha');
        % alpha(scttr,.15);
        
        plot(debug.x(1,:),debug.x(2,:),colors{part_no_});
        plot(debug.x(1,1),debug.x(2,1),strcat(colors{part_no_},'s'),...
            'MarkerFaceColor',colors{part_no_});

        plot(debug.x(1,end),debug.x(2,end),strcat(colors{part_no_},'o'),...
            'MarkerFaceColor',colors{part_no_},'MarkerSize',3);
    end

    for part_no_=1:length(colors)

        load(parts{part_no+part_no_-1},'Mu');

        for j=1:size(Mu,2) 
            plot(Mu(1,j),Mu(2,j),strcat(colors{part_no_},'s')); % plotting 
                                                                % Mus
        end
    end

    pause(.5);
    set(fig,'units','points','position',...
        [FIXED_SPACE+mod(plot_no,6)*(WIDTH+35),...
         SCREEN_HEIGHT-floor(plot_no/6)*(HEIGHT+110),...
         WIDTH,...
         HEIGHT]);
    plot_no=plot_no+1;
    saveas(fig,strcat('data/','fig_',...
        string(2*plot_no-1),'+',string(2*plot_no)),'svg');
end

disp('done plotting. Press a key to close...')
pause;

close all;

