ENUlabel = ['E' 'N' 'U'];

 for bindex = 1
     % Define B, A, C, as matrix of seconds x beam:
     dims = [Data.heading(dye,bindex), Data.pitch(dye,bindex), Data.roll(dye,bindex)];
     %V = [Data.v1(dye,bindex), Data.v2(dye,bindex), Data.v3(dye,bindex)];
     %ENU = [Data.east(dye,bindex), Data.north(dye,bindex), Data.up(dye,bindex)];
     %A = [Data.a1(dye,bindex), Data.a2(dye,bindex), Data.a3(dye,bindex)];
     %C = [Data.c1(dye,bindex), Data.c2(dye,bindex), Data.c3(dye,bindex)];
     for beam = 1
     % Generate a figure
         tsf = figure('name',[ 'Beam ' num2str(beam) ' Bin ' num2str(bindex)],'NumberTitle', 'off');
         
         %histfig = figure('name',[ 'Beam ' num2str(beam) ' Bin ' num2str(bindex)],'NumberTitle', 'off');
     % Create Axes using subplot
         figure(tsf);
         ax1 = subplot(3,1,1);
         ax2 = subplot(3,1,2);
         ax3 = subplot(3,1,3);
         %ax4 = subplot(4,1,4);     
     % Plot data
%          plot(ax1,ENU(:,beam),'b.');
%          ylim(ax1,[-0.5 0.5]);
%          ylabel(ax1, ['Velocity, '  ENUlabel(beam)  '(m/s)']);
%          
         plot(ax1,dims(:,1),'g');
         ylim(ax1,[-360 360]);
         ylabel(ax1, 'nuat heading');
         
         plot(ax2,dims(:,2),'g');
         ylim(ax2,[-180 180]);
         ylabel(ax2, 'pitch');
         
         plot(ax3,dims(:,3),'g');
         ylim(ax3, [-180 180])
         ylabel(ax3, 'roll');
         xlabel(ax3, 'Time (samples)');
         
        hold on
         
         fetch_M1(1)
         dims = [Data.heading(dye,bindex), Data.pitch(dye,bindex), Data.roll(dye,bindex)];
         figure(tsf);
         plot(ax1,dims(:,1),'b');
         
         
         plot(ax2,dims(:,2),'b');
         
         
         plot(ax3,dims(:,3),'b');
         
         hold off
     % Export Figure
%          figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\timeseries\"+ version + "_ts_r" + releaseNum + "_beam" + beam + "_bin" + bindex + ".pdf"];
%          exportgraphics(gcf, figname);
%          close(gcf)

     % plot histogram
%          figure(histfig);
%          histogram(dims(:,beam),[-0.75:0.01:0.75]);
%          ylabel('Counts')
%          xlabel('Velocity (m/s)')
     % Export Histogram
%          figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\histogram\" + version + "_hist_r" + releaseNum + "_beam" + beam + "_bin" + bindex + ".pdf"];
%          exportgraphics(gcf, figname);
%          close(gcf)
     end
 end

