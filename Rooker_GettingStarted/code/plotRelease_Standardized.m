ENUlabel = ['E' 'N' 'U'];
dt = datetime(datestr(Data.Time(dye)));
 for bindex = binStart:binEnd
     % Define B, A, C, as matrix of seconds x beam:
     B = [Data.Velocity_Beam1(dye,bindex), Data.Velocity_Beam2(dye,bindex), Data.Velocity_Beam3(dye,bindex)];
     %V = [Data.v1(dye,bindex), Data.v2(dye,bindex), Data.v3(dye,bindex)];
     ENU = [Data.Velocity_East(dye,bindex), Data.Velocity_North(dye,bindex), Data.Velocity_Up(dye,bindex)];
     A = Data.Amplitude_Minimum(dye,bindex);
     C = Data.Correlation_Minimum(dye,bindex);
     for beam = 1:3
     % Generate a figure
         tsf = figure('name',[ 'Beam ' num2str(beam) ' Bin ' num2str(bindex)],'NumberTitle', 'off');
         histfig = figure('name',[ 'Beam ' num2str(beam) ' Bin ' num2str(bindex)],'NumberTitle', 'off');
     % Create figure
         figure(tsf);     
     % Plot data
         subplot(4,1,1);
         plot(dt,ENU(:,beam),'b.');
         ylim([-0.5 0.5]);
         ylabel(['Velocity, '  ENUlabel(beam)  '(m/s)']);
         
         subplot(4,1,2);
         plot(dt,B(:,beam),'m.');
         ylim([-0.5 0.5]);
         ylabel('Beam Velocity (m/s)');
         
         subplot(4,1,3);
         plot(dt,A,'r.');
         ylim([30 200]);
         ylabel('Signal Strength (dB/dB)');
         
         subplot(4,1,4);
         plot(dt, C,'g.');
         ylim([70 100])
         ylabel('Pulse Correlation (%)');
         xlabel('Time (samples)');
     % Export Figure
         figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\timeseries\"+ version + "_ts_r" + releaseNum + "_beam" + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)

     % plot histogram
         figure(histfig);
         histogram(ENU(:,beam),[-0.75:0.01:0.75]);
         ylabel('Counts')
         xlabel(['Velocity, '  ENUlabel(beam)  '(m/s)'])
     % Export Histogram
         figname = ["C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\figures\Release" + releaseNum + "\histogram\" + version + "_hist_r" + releaseNum + "_beam" + beam + "_bin" + bindex + ".pdf"];
         exportgraphics(gcf, figname);
         close(gcf)
     end
 end

