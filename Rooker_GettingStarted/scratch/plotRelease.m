% Load the data

% load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\Release1\raw\KELP1_AquadoppHR_raw.mat")
% load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_AquadoppHR_L0.mat")
% load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release2\L0\KELP2_Aquadopp_L0.mat")
% load("C:\Users\jkr6136\OnecdDrive - UNC-Wilmington\Kelp_data\data\Release1\raw\KELP1_AquadoppHR_raw.mat")

function plotRelease(inputPath, outputFolder, startDate, endDate);
    if ~exist('outputFolder', 'dir')
        mkdir(outputFolder);
    end
load(inputPath)
% this function plots beam, velocity, aplitude, and correlation over time
% and saves the figure in the given output file

if exist('time')
    dye = find((time>=datenum(startDate)) & (time<=datenum(endDate)));
else
    dye = find((date>=datenum(startDate)) & (date<=datenum(endDate)));
end 
outputFolder
 for bindex = 1:size(b1,2)
     % Define B, A, C, as matrix of seconds x beam:
     B = [b1(dye,bindex), b2(dye,bindex), b3(dye,bindex)];
     V = [v1(dye,bindex), v2(dye,bindex), v3(dye,bindex)];
     A = [a1(dye,bindex), a2(dye,bindex), a3(dye,bindex)];
     C = [c1(dye,bindex), c2(dye,bindex), c3(dye,bindex)];
     for beam = 1:3
     % Generate a figure
         tsf = figure('name',[ 'Beam ' num2str(beam)],'Visible', 'off');
         histfig = figure('name',[ 'Beam ' num2str(beam)]);
     % Create Axes using subplot
         figure(tsf);
         ax1 = subplot(4,1,1);
         ax2 = subplot(4,1,2);
         ax3 = subplot(4,1,3);
         ax4 = subplot(4,1,4);
     % Plot data
         plot(ax1,B(:,beam),'b.');
         ylim(ax1,[-0.5 0.5]);
         ylabel(ax1, 'Velocity (m/s)');
         
         plot(ax2,B(:,beam),'m.');
         ylim(ax2,[-0.5 0.5]);
         ylabel(ax2, 'Velocity (m/s)');
         
         plot(ax3,A(:,beam),'r.');
         ylabel(ax3, 'Amplitude');
         
         plot(ax4,C(:,beam),'g.');
         ylabel(ax4, 'Correlation (%)');
         xlabel(ax4, 'Time (s)');
     % Export Figure
         figname = [string(outputFolder) + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)

     % plot histogram
         figure(histfig);
         histogram(V(:,beam),[-0.75:0.01:0.75]);
         ylabel('Counts')
         xlabel('Velocity (m/s)')
     % Export Histogram
         figname = [string(outputFolder) + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)
     end
 end
end 