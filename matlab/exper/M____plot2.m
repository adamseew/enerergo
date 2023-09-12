
% Generates the "intermittend exploration plot" in "Energy-Aware 
% Ergodic Search: Continuous Exploration for Multi-agent Systems with 
% Battery Constraints" (Fig. 4)

% ! run after generating all the parts in data (run M____eneaware2.m first)


answer=questdlg('Would you like to clean the environment?',...
    'Clean',...
	'Yes','No','No');
switch answer
case 'Yes'
    clear;
end
clear("answer");

parts_={'data/part_no_1.mat', 'data/part_no_2.mat'};

parts={'data2/part_no_1.mat', 'data2/part_no_2.mat',...
       'data2/part_no_3.mat', 'data2/part_no_4.mat',...
       'data2/part_no_5.mat', 'data2/part_no_6.mat',...
       'data2/part_no_7.mat', 'data2/part_no_8.mat',...
       'data2/part_no_9.mat', 'data2/part_no_10.mat',...
       'data2/part_no_11.mat','data2/part_no_12.mat',...
       'data2/part_no_13.mat','data2/part_no_14.mat'};
colors={'b','r'};
[X_,Y_]=meshgrid(0:.005:1);
plot_no=0;

WIDTH=122;
HEIGHT=WIDTH;
FIXED_SPACE=220;
SCREEN_HEIGHT=get(0,'screensize')*.75;
SCREEN_HEIGHT=SCREEN_HEIGHT(4)-FIXED_SPACE;

disp('started plotting')

for part_no=1:length(parts)

    if mod(part_no-1,4)~=0
        continue;
    end

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

    if part_no>=2
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
        debug.alpha(2)*mvnpdf([X_(:) Y_(:)],Mu(:,2)',Sigma(:,:,2))+...
        debug.alpha(3)*mvnpdf([X_(:) Y_(:)],Mu(:,3)',Sigma(:,:,3))+...
        debug.alpha(4)*mvnpdf([X_(:) Y_(:)],Mu(:,4)',Sigma(:,:,4));
    p=reshape(p,size(X_));
    pclr=pcolor(X_,Y_,p);
    shading interp;
    set(pclr,'facealpha',0.65);
    colormap(flipud(bone(15))); % adds the prob. distributions (gauss. mix.
                                % in background)
    load(parts{part_no},'debug');
        
    plot(debug.x(1,:),debug.x(2,:),'b');
    plot(debug.x(1,1),debug.x(2,1),'bs',...
        'MarkerFaceColor','b');

    plot(debug.x(1,end),debug.x(2,end),'bo',...
        'MarkerFaceColor','b','MarkerSize',3);

    for part_no_=1:length(colors)

        load(parts_{part_no_},'Mu');

        for j=1:size(Mu,2) 
            plot(Mu(1,j),Mu(2,j),strcat('c','s')); % plotting 
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
    set(fig,'Renderer','painters');
    saveas(fig,strcat('data2/','fig_',...
        string(plot_no)),'svg');
end

disp('done plotting. Press a key to close...')
pause;

close all;

