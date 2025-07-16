% working on some Dye/Temp analysis

%load("DyeReleaseLanderData.mat")
%S = load("S/L1/mooring_S.mat");

S.Temperature = S.Temperature';
meanTemp = mean(R1.Temperature(:, 3));
Tprime = abs(S.Temperature(:,2)-meanTemp);


figure;
ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
plot(ax1, datetime(S.Time, 'ConvertFrom', 'datenum'), Tprime, 'y')
xline(ax1, datetime('03-Jul-2024 18:38:00'), 'g')
xline(ax1, datetime('03-Jul-2024 19:52:00'), 'r')
plot(ax2, datetime(S.Time, 'ConvertFrom', 'datenum'), S.Dye, 'b')
xline(ax2, datetime('03-Jul-2024 18:38:00'), 'g')
xline(ax2, datetime('03-Jul-2024 19:52:00'), 'r')
linkaxes([ax1 ax2], 'x')