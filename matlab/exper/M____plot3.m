

fig=figure;
set(fig,'color','w');
ax1=gca;
ax1.YGrid='on';
ax1.Layer='top';
ax1.GridLineStyle=':';
ax1.GridAlpha=.5;
set(ax1,'xdir','reverse');
grid on;
hold on;

curves=[];
maxsize_=inf;
load('data/part_no_1.mat','debug');
debug0=debug;
for k=1:2:28
    load(strcat('data/part_no_',string(k),'.mat'),'debug');
    debug1=debug;
    load(strcat('data/part_no_',string(k+1),'.mat'),'debug');
    debug2=debug;
    
    if min([length(debug1.epsilon) length(debug2.epsilon)])<maxsize_
        maxsize_=min([length(debug1.epsilon) length(debug2.epsilon)]);
    end
end

for k=1:2:28
    load(strcat('data/part_no_',string(k),'.mat'),'debug');
    debug1=debug;
    load(strcat('data/part_no_',string(k+1),'.mat'),'debug');
    debug2=debug;

    maxsize=min([length(debug1.epsilon) length(debug2.epsilon)]);
    plot(debug0.z(1:maxsize_)/120,.5*debug1.epsilon(:,1:maxsize_)+.5*debug2.epsilon(:,1:maxsize_),...
        'color','r')
    curves=[curves;...
        .5*debug1.epsilon(:,1:maxsize_)+.5*debug2.epsilon(:,1:maxsize_)];
end
max_curves=max(curves,[],1);
min_curves=min(curves,[],1);

ylim([0 .15]);

curve_avg = avgcurve(maxsize_, {[debug0.z(1:maxsize_)/120 curves(1,:)'],...
    [debug0.z(1:maxsize_)/120 curves(2,:)'],...
    [debug0.z(1:maxsize_)/120 curves(3,:)'],...
    [debug0.z(1:maxsize_)/120 curves(4,:)'],...
    [debug0.z(1:maxsize_)/120 curves(5,:)'],...
    [debug0.z(1:maxsize_)/120 curves(6,:)'],...
    [debug0.z(1:maxsize_)/120 curves(7,:)'],...
    [debug0.z(1:maxsize_)/120 curves(8,:)'],...
    [debug0.z(1:maxsize_)/120 curves(9,:)'],...
    [debug0.z(1:maxsize_)/120 curves(10,:)'],...
    [debug0.z(1:maxsize_)/120 curves(11,:)'],...
    [debug0.z(1:maxsize_)/120 curves(12,:)'],...
    [debug0.z(1:maxsize_)/120 curves(13,:)'],...
    [debug0.z(1:maxsize_)/120 curves(14,:)']});



WIDTH=132;
HEIGHT=WIDTH*2/3;
FIXED_SPACE=220;
SCREEN_HEIGHT=get(0,'screensize')*.75;
SCREEN_HEIGHT=SCREEN_HEIGHT(4)-FIXED_SPACE;

pause(.5);
set(fig,'units','points','position',...
    [FIXED_SPACE+mod(1,6)*(WIDTH+35),...
     SCREEN_HEIGHT-floor(1/6)*(HEIGHT+110),...
     WIDTH,...
     HEIGHT]);
set(fig,'Renderer','painters');
saveas(fig,strcat('data/','fig_erg-soc_1'),'svg');


fig2=figure;
set(fig2,'color','w');

ax1=gca;
ax1.YGrid='on';
ax1.Layer='top';
ax1.GridLineStyle=':';
ax1.GridAlpha=.5;
set(ax1,'xdir','reverse');
grid on;
hold on;
ylim([0 .15]);

x2=[debug0.z(1:maxsize_)/120; flip(debug0.z(1:maxsize_)/120,1)];
y2=[min_curves'; flip(max_curves',1)];
hf=fill(x2,y2,'b');
set(hf,'facealpha',.4)

plot(debug0.z(1:maxsize_)/120,curve_avg(:,2),'color','b','MarkerSize',3)

WIDTH=132;
HEIGHT=WIDTH/3;
FIXED_SPACE=220;
SCREEN_HEIGHT=get(0,'screensize')*.75;
SCREEN_HEIGHT=SCREEN_HEIGHT(4)-FIXED_SPACE;

pause(.5);
set(fig2,'units','points','position',...
    [FIXED_SPACE+mod(1,6)*(WIDTH+35),...
     SCREEN_HEIGHT-floor(1/6)*(HEIGHT+110),...
     WIDTH,...
     HEIGHT]);
set(fig2,'Renderer','painters');
saveas(fig2,strcat('data/','fig_erg-soc_2'),'svg');
