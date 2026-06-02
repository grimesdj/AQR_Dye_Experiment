function [y, fsd, idx] = hamming_filter(x, wpass, fs, fig, ds)
    % 
    % hamming_filter
    %     applies a low-pass filter to data using a hamming filter
    % 
    % USAGE: [y, fsd] = hamming_filter(t, x, wpass, fs, fig, ds)
    %        [y, fsd] = hamming_filter(x, wpass, fs);
    % 
    %     INPUTS:
    %             x     = vector data to be low-passed
    %             wpass = frequency cutoff (Hz)
    %             fs    = original sampling rate (Hz) (defaults to 1 Hz)
    %             fig   = (optional) figure flag: 1 for figures, 0 for none (defaults to 0)
    %             ds    = (optional) downsample flag: 1 for downsampling, 0 for none (defaults to 1)
    %             
    % 
    %     OUTPUTS:
    %             y     = filtered data vector
    %             fsd   = new sample rate (Hz)
    %             D     = Downsample factor
    %             
    % 
    %             written by Jason Rooker, May 2026


                
% Defaults


if nargin < 5 || isempty(ds)
    ds = 1;

    if nargin < 4 || isempty(fig)
        fig = 0;

        if nargin < 3 || isempty(fs)
            fs = 1;
        end
    end
end

% window length
N = round((1/wpass) * fs);

w = hamming(N);
w = w/sum(w);

mask = ~isnan(x);

x0 = x;
x0(~mask) = 0;

% normalization factor
norm = conv(double(mask),w,'same');
y = conv(x0,w,'same') ./ norm;

% NaN-dense filter to protect statistical integrety for burst data
y(norm < 0.5) = NaN;

% downsample
if ds
    D = max(1, round(fs / (4 * wpass)));
else
    D = 1;
end

idx = 1:D:length(x);
y = y(idx);
fsd = fs/D;

% (optional) figs to make sure it worked
if fig

    xfig = 1:length(x);
    figure
    ax1 = subplot(2, 1, 1);
    plot(xfig, x, '.', 'LineWidth', 2)
    ylabel({'Original', 'Data'})
    set(gca, 'FontSize', 18)
    grid minor
    
    ax2 = subplot(2, 1, 2);
    plot(xfig(1:D:end), y, '.', 'LineWidth', 2)
    ylabel({'Original', 'filtered'})
    set(gca, 'FontSize', 18)
    ylim(ax1.YLim)
    linkaxes([ax1 ax2], 'x')
    grid minor

end

end
