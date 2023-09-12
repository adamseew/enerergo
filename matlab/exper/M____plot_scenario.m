
% Generates the "scenario plot" in "Energy-Aware 
% Ergodic Search: Continuous Exploration for Multi-agent Systems with 
% Battery Constraints" (Fig. 3)

% ! run after generating all the parts in data (run MAIN_eneaware.m first)


answer=questdlg('Would you like to clean the environment?',...
    'Clean',...
	'Yes','No','No');
switch answer
case 'Yes'
    clear;
end
clear("answer");


parts={'data/part_no_1.mat', 'data/part_no_2.mat'};
colors={'c','c'};
[X_,Y_]=meshgrid(0:.005:1);
plot_no=0;

WIDTH=122;
HEIGHT=WIDTH;
FIXED_SPACE=220;
SCREEN_HEIGHT=get(0,'screensize')*.75;
SCREEN_HEIGHT=SCREEN_HEIGHT(4)-FIXED_SPACE;

uav_x0=[[.1;.3] ...
        [.9;.7] ...
        [.1;.7] ...
        [.9;.3]]; % charging stations


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

    load(parts{part_no},'Sigma','Mu','debug');
    p=  .25*mvnpdf([X_(:) Y_(:)],Mu(:,1)',Sigma(:,:,1))+...
        .25*mvnpdf([X_(:) Y_(:)],Mu(:,2)',Sigma(:,:,2));
    load(parts{part_no+1},'Sigma','Mu','debug');
    p=p+.25*mvnpdf([X_(:) Y_(:)],Mu(:,1)',Sigma(:,:,1))+...
        .25*mvnpdf([X_(:) Y_(:)],Mu(:,2)',Sigma(:,:,2));
    p=reshape(p,size(X_));
    pclr=pcolor(X_,Y_,p);
    shading interp;
    set(pclr,'facealpha',0.65);
    colormap(flipud(bone(15))); % adds the prob. distributions (gauss. mix.
                                % in background)

    for part_no_=1:length(colors)

        load(parts{part_no+part_no_-1},'Mu');

        for j=1:size(Mu,2) 
            plot(Mu(1,j),Mu(2,j),strcat(colors{part_no_},'s')); % plotting 
                                                                % Mus
        end
    end

    %plot(uav_x0(1,1),uav_x0(2,1),strcat(colors{1},'o'),...
    %    'MarkerFaceColor',colors{1},'MarkerSize',3)
    %plot(uav_x0(1,2),uav_x0(2,2),strcat(colors{1},'o'),...
    %    'MarkerFaceColor',colors{1},'MarkerSize',3)
    %plot(uav_x0(1,3),uav_x0(2,3),strcat(colors{1},'o'),...
    %    'MarkerFaceColor',colors{1},'MarkerSize',3)
    %plot(uav_x0(1,4),uav_x0(2,4),strcat(colors{1},'o'),...
    %    'MarkerFaceColor',colors{1},'MarkerSize',3)

    pause(.5);
    set(fig,'units','points','position',...
        [FIXED_SPACE+mod(plot_no,6)*(WIDTH+35),...
         SCREEN_HEIGHT-floor(plot_no/6)*(HEIGHT+110),...
         WIDTH,...
         HEIGHT]);
    plot_no=plot_no+1;
    set(fig,'Renderer','painters');
    saveas(fig,strcat('data/','fig_scenario'),'svg');
end

disp('done plotting. Press a key to close...')
pause;

close all;



