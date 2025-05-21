% Load the L0 data
load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_AquadoppHR_L0.mat")
%load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Aquadopp_L0.mat")
% time to plot the b1:b3 data as speed vs time (remeber time is in decimal
% days)
%columns are bins (depth), rows are time at 600 second intervals, and
%values are velocity

% Limit data to during dye release


dye = find((date>=datenum('03-Jul-2024 18:38:00')) & (date>=datenum('03-Jul-2024 19:52:00')));
 for bindex = 1:size(b1,2)
     % Define B, A, C, as matrix of seconds x beam:
     B = [b1(dye,bindex), b2(dye,bindex), b3(dye,bindex)];
     A = [a1(dye,bindex), a2(dye,bindex), a3(dye,bindex)];
     C = [c1(dye,bindex), c2(dye,bindex), c3(dye,bindex)];
     for beam = 1:3
     % Generate a figure
         tsf = figure('name',[ 'Beam ' num2str(beam)],'Visible', 'off');
         histfig = figure('name',[ 'Beam ' num2str(beam)]);
     % Create Axes using subplot
         figure(tsf);
         ax1 = subplot(3,1,1);
         ax2 = subplot(3,1,2);
         ax3 = subplot(3,1,3);
     % Plot data
         plot(ax1,B(:,beam),'b.');
         ylim(ax1,[-0.5 0.5]);
         ylabel(ax1, 'Velocity (m/s)');
         plot(ax2,A(:,beam),'r.');
         ylabel(ax2, 'Amplitude');
         plot(ax3,C(:,beam),'g.');
         ylabel(ax3, 'Correlation (%)');
         xlabel(ax3, 'Time (s)');
     % Export Figure
         figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release1\timeseries\ts_r1_beam" + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)

     % plot histogram
         figure(histfig);
         histogram(B(:,beam),[-0.75:0.01:0.75]);
         ylabel('Counts')
         xlabel('Velocity (m/s)')
     % Export Histogram
         figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release1\histogram\hist_r1_beam" + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)
     end
 end
