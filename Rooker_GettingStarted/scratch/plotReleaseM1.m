% 
% 
% function plotRelease(releaseNum, version) 
% % calls for release number as int
% %   version as raw/L@
% 
% if nargin < 2
%     version = 'path';
% end
% 
% % Load the data
% 
% 
% if releaseNum == 1
% 
%     if strcmp(version, 'raw')
%         Data = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat");
%     elseif strcmp(version, 'L0')
%         Data = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_AquadoppHR_L0.mat");
%     else
%         disp('Invalid Version')
%     end
%     
% elseif releaseNum == 2 || releaseNum == 3
% 
%     if strcmp(version, 'L0')
%         %Data = load(%Release 2 path);
%     elseif strcmp(version, 'L1')
%         %Data = load(%Release 2 path);
%     else
%         disp('Invalid Version')   
%     end
%     
% elseif releaseNum == 0
%     Data = version;
%     version = 'out';
% else
%     Data = load(releaseNum);
% 
% end
% 
% 
% % Limit data to during dye release
TRange = readtable("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\info\dye_mixing_cals_and_releases\dye_release_times.csv");
% 
% if ~exist('Data.time', 'var')
%     time = Data.date;
% else 
%     time = Data.time;
% end
% if releaseNum > 0
%dye = find((Time>=datenum(TRange.StartTime_UTC_(1))) & (time<=datenum(TRange.EndTime_UTC_(1))));
% else
%    dye = ':';
% end
% 
% 
% % time to plot the b1:b3 data as speed vs time (remeber time is in decimal
% % days)
% % columns are bins (depth), rows are time at 600 second intervals, and
% % values are velocity
releaseNum = 'M1';

version = 'out';
Data = load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\2024_PROCESSED_DATA\M1\L0\ADCP\ADCP_M1_002.mat");
dye = find((Data.Time>=datenum(TRange.StartTime_UTC_(1))) & (Data.Time<=datenum(TRange.EndTime_UTC_(1))));
% else
 for bindex = 1:size(Data.Velocity_East,1)
     % Define B, A, C, as matrix of seconds x beam:
     %%B = [Data.b1(dye,bindex), Data.b2(dye,bindex), Data.b3(dye,bindex)];
     %V = [Data.Velocity_East(dye,bindex), Data.Velocity_North(dye,bindex), Data.Velocity_Up(dye,bindex)];
     ENU = [Data.Velocity_East(bindex,:); Data.Velocity_North(bindex, :); Data.Velocity_Up(bindex, :)];
     %A = [Data.(dye,bindex), Data.a2(dye,bindex), Data.a3(dye,bindex)];
     %C = [Data.c1(dye,bindex), Data.c2(dye,bindex), Data.c3(dye,bindex)];
     ENU = ENU';
     for beam = 1:3
     % Generate a figure
         tsf = figure('name',[ 'Beam ' num2str(beam) ' Bin ' num2str(bindex)],'NumberTitle', 'off');
         histfig = figure('name',[ 'Beam ' num2str(beam) ' Bin ' num2str(bindex)],'NumberTitle', 'off');
     % Create Axes using subplot
         figure(tsf);
         ax1 = subplot(4,1,1);
         ax2 = subplot(4,1,2);
         ax3 = subplot(4,1,3);
         ax4 = subplot(4,1,4);
     % Plot data
         plot(ax1,ENU(:, beam),'b.');
         ylim(ax1,[-0.5 0.5]);
         ylabel(ax1, 'ENU Velocity(m/s)');
         
%          plot(ax2,B(:,beam),'m.');
%          ylim(ax2,[-0.5 0.5]);
%          ylabel(ax2, 'Velocity (m/s)');
%          
%          plot(ax3,A(:,beam),'r.');
%          ylim(ax3,[0 200]);
%          ylabel(ax3, 'Amplitude');
%          
%          plot(ax4,C(:,beam),'g.');
%          ylabel(ax4, 'Correlation (%)');
%          xlabel(ax4, 'Time (s)');
     % Export Figure
         figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\timeseries\"+ version + "_ts_r" + releaseNum + "_beam" + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)

     % plot histogram
         figure(histfig);
         histogram(ENU(:,beam),[-0.75:0.01:0.75]);
         ylabel('Counts')
         xlabel('Velocity (m/s)')
     % Export Histogram
         figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\histogram\" + version + "_hist_r" + releaseNum + "_beam" + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)
     end
 end
