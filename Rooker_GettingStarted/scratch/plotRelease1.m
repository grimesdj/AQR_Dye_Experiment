% Load the L0 data
load("C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\data\Release1\L0\KELP1_AquadoppHR_L0.mat")
% time to plot the b1:b3 data as speed vs time (remeber time is in decimal
% days)
%columns are bins (depth), rows are time at 600 second intervals, and
%values are velocity

B = [b1(:,1), b2(:,1), b3(:,1)];
A = [a1(:,1), a2(:,1), a3(:,1)];
C = [c1(:,1), c2(:,1), c3(:,1)];
for x = 1:3
% Generate a figure
    figure('name',[ 'Beam ' num2str(x)]);
% Create Axes using subplot
    ax1 = subplot(3,1,1);
    ax2 = subplot(3,1,2);
    ax3 = subplot(3,1,3);
% Plot data
    plot(ax1,B(:,x),'b.')
    ylim(ax1,[-0.5 0.5]);
    ylabel(ax1, 'Velocity (m/s)')
    plot(ax2,A(:,x),'r.')
    ylabel(ax2, 'Amplitude')
    plot(ax3,C(:,x),'g.')
    ylabel(ax3, 'Correlation (%)')
    xlabel(ax3, 'Time (s)')
end