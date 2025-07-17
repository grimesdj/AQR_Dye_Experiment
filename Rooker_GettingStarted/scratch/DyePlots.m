% working on some Dye/Temp analysis

load("../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat")
S = load("../../../../Kelp_data/data/2024_PROCESSED_DATA/S/L1/mooring_S.mat");

S.Temperature = S.Temperature';
meanTemp = mean(R1.Temperature(:, 3));
Tprime = abs(S.Temperature(:,2)-meanTemp);

% %
% R1.Start = datetime('03-Jul-2024 18:38:00');
% R1.End   = datetime('03-Jul-2024 19:52:00');
% %
% R2.Start = datetime('08-Jul-2024 17:41:00');
% R2.End   = datetime('08-Jul-2024 20:08:00');
% %
% R3.Start = datetime('11-Jul-2024 17:28:00');
% R3.End   = datetime('11-Jul-2024 19:55:00');
% %

save("../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat")


figure;
ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
plot(ax1, datetime(S.Time, 'ConvertFrom', 'datenum'), Tprime, 'y')
xline(ax1, datetime('03-Jul-2024 18:38:00'), 'g')
xline(ax1, datetime('03-Jul-2024 19:52:00'), 'r')
plot(ax2, datetime(S.Time, 'ConvertFrom', 'datenum'), SW.Dye, 'b')
xline(ax2, datetime('03-Jul-2024 18:38:00'), 'g')
xline(ax2, datetime('03-Jul-2024 19:52:00'), 'r')
linkaxes([ax1 ax2], 'x')