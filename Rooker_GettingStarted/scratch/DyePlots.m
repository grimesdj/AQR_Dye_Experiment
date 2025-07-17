% working on some Dye/Temp analysis
clear all
close all

load("../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat")

% Release 1 Temperature
T_bar = mean(R1.Temperature(:,3));

% Get files by mooring
files = dir('../../../../Kelp_data/data/2024_PROCESSED_DATA/');
sites = {'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'};

% Filter out . and .. and non-directories
files = files([files.isdir]);
files = files(~ismember({files.name}, {'.', '..'}));

for i = 1:length(files)
    if ismember(files(i).name, sites)
        matFile = dir(fullfile(files(i).folder, files(i).name, 'L1', '*.mat'));
        if ~isempty(matFile)
            Moor = load(fullfile(matFile(1).folder, matFile(1).name));            
            idx = find(Moor.Temperature_mab == 1.0);
            Tprime = abs(T_bar-Moor.Temperature(idx, :));
            figure(i);
            ax1 = subplot(2, 1, 1);
            ax2 = subplot(2, 1, 2);
            
            plot(ax1, datetime(Moor.Time, 'ConvertFrom', 'datenum'), Tprime, 'y')
            ylabel(ax1, "T' (^{\circ}C)")
            
                xline(ax1, R1.Start, 'g')
                xline(ax1, R1.End, 'r')
            
                xline(ax1, R2.Start, 'g')
                xline(ax1, R2.End, 'r')
            
                xline(ax1, R3.Start, 'g')
                xline(ax1, R3.End, 'r')
            
            plot(ax2, datetime(Moor.Time, 'ConvertFrom', 'datenum'), Moor.Dye(1,:), 'b')
            ylabel(ax2, 'Dye Conc (ppb)')

                xline(ax2, R1.Start, 'g')
                xline(ax2, R1.End, 'r')
            
                xline(ax2, R2.Start, 'g')
                xline(ax2, R2.End, 'r')
            
                xline(ax2, R3.Start, 'g')
                xline(ax2, R3.End, 'r')
            
            linkaxes([ax1 ax2], 'x')
            title(ax1, [Moor.Name ' mooring'])
            xlim([R1.Start '06-Jul-2024 00:00:00'])
            exportgraphics(gcf, ['../../../../Kelp_data/Summer2025/Rooker/figures/Dye_Concentrations/TempDye_' Moor.Name '.pdf'] )
            
        end
    end
end
        




%%%% Command line Vomit

% Clean up data
valid = ~isnan(CTD1.dye) & ~isnan(CTD1.temp);
T = CTD1.temp(valid);
D = CTD1.dye(valid);

% Normalize dye values to make them sum to 1 (i.e., turn into weights)
w = D / nansum(D);

% Compute KDE of T, weighted by dye
[x, f] = ksdensity(T, 'Weights', w);  % <-- this is key

% Plot
plot(x, f, 'LineWidth', 2)
xlabel('Temperature (°C)')
ylabel('Weighted PDF of Dye')
title('PDF of Temperature Weighted by Dye Concentration')
grid on

% Clean up data
valid = ~isnan(CTD1.dye) & ~isnan(CTD1.temp);
T = CTD1.temp(valid);
D = CTD1.dye(valid);
% Normalize dye values to make them sum to 1 (i.e., turn into weights)
w = D / nansum(D);
% Compute KDE of T, weighted by dye
[x, f] = ksdensity(T, 'Weights', w);  % <-- this is key
% Plot
plot(f, x, 'LineWidth', 2)
xlabel('Temperature (°C)')
ylabel('Weighted PDF of Dye')
title('Probability of Dye')
datetime(CTD1.time(1), 'ConvertFrom', 'datenum')

ans = 

  datetime

   03-Jul-2024 20:46:47

plot(f-mean(R1.Temperature(:,3)), x, 'g')
load("../../../../Kelp_data/data/2024_PROCESSED_DATA/DyeReleaseLanderData.mat")
hold on, plot(f-mean(R1.Temperature(:,3)), x, 'g')
help ksdensity
 ksdensity - Kernel smoothing function estimate for univariate and bivariate data
    This MATLAB function returns a probability density estimate, f, for the
    sample data in the vector or two-column matrix x.

    Syntax
      [f,xi] = ksdensity(x)
      [f,xi] = ksdensity(x,pts)
      [f,xi] = ksdensity(___,Name,Value)
      [f,xi,bw] = ksdensity(___)

      ksdensity(___)
      ksdensity(ax,___)

    Input Arguments
      x - Sample data
        column vector | two-column matrix
      pts - Points at which to evaluate f
        vector | two-column matrix
      ax - Axes handle
        handle

    Name-Value Arguments
      Bandwidth - Bandwidth of kernel smoothing window
        "normal-approx" (default) | "plug-in" | scalar value |
        two-element vector
      BoundaryCorrection - Boundary correction method
        'log' (default) | 'reflection'
      Censoring - Logical vector
        vector of 0s (default) | vector of 0s and 1s
      Function - Function to estimate
        'pdf' (default) | 'cdf' | 'icdf' | 'survivor' | 'cumhazard'
      Kernel - Type of kernel smoother
        'normal' (default) | 'box' | 'triangle' | 'epanechnikov' |
        function handle | character vector | string scalar
      NumPoints - Number of equally spaced points
        100 (default) | scalar value
      Support - Support for density
        "unbounded" (default) | "positive" | "nonnegative" | "negative" |
        two-element vector | two-by-two matrix
      PlotFcn - Function used to create kernel density plot
        'surf' (default) | 'contour' | 'plot3' | 'surfc'
      Weights - Weights for sample data
        vector

    Output Arguments
      f - Estimated function values
        vector
      xi - Evaluation points
        100 equally spaced points | 900 equally spaced points | vector |
        two-column matrix
      bw - Bandwidth of smoothing window
        scalar value

    Examples
      Estimate Density
      Estimate Density with Boundary Correction
      Estimate Cumulative Distribution Function at Specified Values
      Plot Estimated Cumulative Distribution Function for Given Number of Points
      Estimate Survivor and Cumulative Hazard for Censored Failure Data
      Estimate Inverse Cumulative Distribution Function for Specified Probability Values
      Return Bandwidth of Smoothing Window
      Plot Kernel Density Estimate of Bivariate Data

    See also histogram, mvksdensity, kde

    Introduced in Statistics and Machine Learning Toolbox before R2006a
    Documentation for ksdensity
    Other uses of ksdensity

sum(x)

ans =

    6.2754

hold on, plot(f-mean(R1.Temperature(:,3)), x/sum(x), 'g')
sum(x/sum(x))

ans =

     1

R1

R1 = 

  struct with fields:

                Time: [8881×1 double]
         Temperature: [8881×5 double]
            Pressure: [8881×5 double]
      PressureOffset: [10.0773 10.1051 10.1338 10.0974 10.2765]
            Salinity: [8881×1 double]
        Conductivity: [8881×1 double]
        Depth_deploy: [8.1571 8.4919 8.7136 8.9453 9.2306]
                 mab: [1.5560 1.2320 1 0.7680 0.5000]
             Density: [8881×5 double]
                  N2: [8881×4 double]
    Temperature_bins: [12 12.2500 12.5000 12.7500 13 13.2500 13.5000 13.7500 14 14.2500 14.5000 14.7500 15 … ] (1×33 double)
     Temperature_pdf: [0 0 0 0 0.0016 0.0393 0.1141 0.0507 0.0707 0.1105 0.1884 0.2128 0.2120 0 0 0 0 0 0 … ] (1×33 double)
               Start: 03-Jul-2024 18:38:00
                 End: 03-Jul-2024 19:52:00

