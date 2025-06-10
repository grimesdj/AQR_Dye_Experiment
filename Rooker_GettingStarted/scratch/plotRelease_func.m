ENUlabel = ['E' 'N' 'U'];

 for bindex = 1:size(Data.b1,2)
     % Define B, A, C, as matrix of seconds x beam:
     B = [Data.b1(dye,bindex), Data.b2(dye,bindex), Data.b3(dye,bindex)];
     %V = [Data.v1(dye,bindex), Data.v2(dye,bindex), Data.v3(dye,bindex)];
     ENU = [Data.east(dye,bindex), Data.north(dye,bindex), Data.up(dye,bindex)];
     A = [Data.a1(dye,bindex), Data.a2(dye,bindex), Data.a3(dye,bindex)];
     C = [Data.c1(dye,bindex), Data.c2(dye,bindex), Data.c3(dye,bindex)];
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
         plot(ax1,ENU(:,beam),'b.');
         ylim(ax1,[-0.5 0.5]);
         ylabel(ax1, ['Velocity, '  ENUlabel(beam)  '(m/s)']);
         
         plot(ax2,B(:,beam),'m.');
         ylim(ax2,[-0.5 0.5]);
         ylabel(ax2, 'Beam Velocity (m/s)');
         
         plot(ax3,A(:,beam),'r.');
         ylim(ax3,[30 200]);
         ylabel(ax3, 'Signal Strength (dB/dB)');
         
         plot(ax4,C(:,beam),'g.');
         ylim(ax4, [70 100])
         ylabel(ax4, 'Pulse Correlation (%)');
         xlabel(ax4, 'Time (samples)');
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

