clc
clear all

load sim_mobilebot.mat   % load Simulation output
data = ans; clear ans
tout = data(1,:);
yout = data(2:13,:);



% FPS = 30;                    % Frames per second
% tout = 0:1/FPS:toutraw(end);
% 
% for n=1:length(tout)
%     yout(:,n) = interp1(toutraw',youtraw',tout(n))';
% end

load dimensions_4_3   % Load dimensions from .mat file create in symbolics script

[F1, V1, C1] = rectverts(RP(4,:),[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0]); %green
[F2, V2, C2] = rectverts(RP(5,:),[1 0 0],[1 0 0],[1 0 0],[1 0 0],[1 0 0],[1 0 0]); %red
[F3, V3, C3] = rectverts(RP(6,:),[1 0 1],[1 0 1],[1 0 1],[1 0 1],[1 0 1],[1 0 1]); %purple

% rotate from nx3 to 3xn and add row of ones to make homogeneous points(4xn)
HV1 = [V1';ones(1,length(V1))];
HV2 = [V2';ones(1,length(V2))];
HV3 = [V3';ones(1,length(V3))];

% OPTIONAL:  read in a .stl file exported from solidworks
% [F, V, C] = STLverts('Link2_Coarse.stl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

s = get(0,'ScreenSize');
h = 0.7*s(4);     % 70% of screen height
dy = (s(4)-h)/2;    % vertical offset
dx = (s(3)-h)/2;    % horizontal offset
figure('Position',[dx dy h h]) %centered
grid on; hold on; axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')
light                           % add a default light
daspect([1 1 1])                % Setting the aspect ratio
view(3)                         % Isometric view
axis([-3 3 -3 3 0 3]);
xlabel ('X (m)','FontSize',16,'Color','r')
ylabel ('Y (m)','FontSize',16,'Color','r')
zlabel ('Z (m)','FontSize',16,'Color','r')

% add a pause and continue button:
h = uicontrol('Position',[20 20 200 40],'String','Continue',...
    'Callback','uiresume(gcbf)');
h = uicontrol('Position',[20 60 200 40],'String','Pause',...
    'Callback','uiwait(gcbf)');
h = uicontrol('Position',[220 20 200 40],'String','Quit',...
    'Callback','close(gcbf)');

dstr = datestr(now,'yyyy_mm_dd_HH_MM_SS'); % create a new filename for each video, to avoid overwriting
namestr = ['mobilebot_video_',dstr]
writerObj = VideoWriter(namestr,'MPEG-4');
set(gca,'nextplot','replacechildren');  % this and the following command required for video recording
set(gcf,'Renderer','zbuffer');
framerate = 30;
writerObj.FrameRate = framerate; 
open(writerObj);

for n=1:length(yout)
           
    % rotate the points using homogeneous transformations
    PV1 = gs4func(yout(1:6,n))*HV1;
    PV2 = gs5func(yout(1:6,n))*HV2;
    PV3 = gs6func(yout(1:6,n))*HV3;
    
    % remove the homogeneous 1 and transpose
    PV1 = PV1(1:3,:)';
    PV2 = PV2(1:3,:)';
    PV3 = PV3(1:3,:)';
    
    % clear old plot before re-patching
    cla
    
    % create patch objects
    p1 = patch('faces', F1, 'vertices' ,PV1);
    p2 = patch('faces', F2, 'vertices' ,PV2);
    p3 = patch('faces', F3, 'vertices' ,PV3);
    
    % Set the face color flat
    set([p1 p2 p3], 'facec', 'flat');
    
    % Set the color (from file)
    set(p1, 'FaceColor', 'r');
    set(p2, 'FaceColor', 'b');
    set(p3, 'FaceColor', 'g');
    
    % Use for transparency (1 equals solid)
    set([p1 p2 p3], 'facealpha',0.5)
    
    % Set the edge color (use 'none' for STL files)
    set([p1 p2 p3], 'EdgeColor','k');    % use none for STL files
    
    % Plot time on screen
    text(-L1,-L1-L2,L1,[num2str(tout(n),'%10.1f'),' s'],'Fontsize',16);
    
    drawnow   % force all elements to draw before capturing frames
    
    writeVideo(writerObj,getframe)
    
end

close(writerObj);

close all;
load sim_mobilebot.mat   % load Simulation output
data = ans; clear ans
tout = data(1,:);
yout = data(2:13,:);

load PIDData.mat
PIDData = ans;
load GoalData.mat
GoalData = ans;
load FloorCollisionData.mat
FloorCollisionData = ans;

figure(2);
hold on;
grid on;

title('\fontsize{24}Position vs Desired Goal');
xlabel('\fontsize{18} Time t (sec)');
ylabel('\fontsize{18} Position \theta_1 - \theta_6 (m or rad)')
plot(tout, data(2,:), 'r-', 'LineWidth', 2);
plot(tout, data(3,:), 'g-', 'LineWidth', 2);
plot(tout, data(4,:), 'b-', 'LineWidth', 2);
plot(tout, data(5,:), 'c-', 'LineWidth', 2);
plot(tout, data(6,:), 'm-', 'LineWidth', 2);
plot(tout, data(7,:), 'y-', 'LineWidth', 2);
plot(tout, GoalData(2,:), 'r--', 'LineWidth', 2);
plot(tout, GoalData(3,:), 'g--', 'LineWidth', 2);
plot(tout, GoalData(4,:), 'b--', 'LineWidth', 2);
plot(tout, GoalData(5,:), 'c--', 'LineWidth', 2);
plot(tout, GoalData(6,:), 'm--', 'LineWidth', 2);
plot(tout, GoalData(7,:), 'y--', 'LineWidth', 2);
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', '\theta_6');
set(gca, 'fontsize',16);

figure(3);
hold on;
title('\fontsize{24}Floor Collision Force (\theta_3)');
xlabel('\fontsize{18} Time t (sec)');
ylabel('\fontsize{18} Force (N)');
plot(FloorCollisionData(1,:), FloorCollisionData(4,:), 'b-', 'LineWidth', 2);
set(gca, 'fontsize',16);

figure(4);
hold on;
title('\fontsize{24}PID Force & Torque');
xlabel('\fontsize{18} Time t (sec)');
ylabel('\fontsize{18} Force/Torque (N, N-m)');
plot(PIDData(1,:), PIDData(2,:), 'r-', 'LineWidth', 2);
plot(PIDData(1,:), PIDData(3,:), 'g-', 'LineWidth', 2);
plot(PIDData(1,:), PIDData(4,:), 'b-', 'LineWidth', 2);
plot(PIDData(1,:), PIDData(5,:), 'c-', 'LineWidth', 2);
plot(PIDData(1,:), PIDData(6,:), 'm-', 'LineWidth', 2);
plot(PIDData(1,:), PIDData(7,:), 'y-', 'LineWidth', 2);
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', '\theta_6');

set(gca, 'fontsize',16);
%plot(tout, data(2,:), 'k-');






