%% Grid Velocity

clear all
close all

moorings = {'M1', 'M2', 'M3'};
for mooring_ID = 1:length(moorings)
    %% Load
    mooring = moorings{mooring_ID};
    
    % ADCP
    fprintf('loading %s data...\n', mooring)
    fpath = fullfile('..', '..', '..', '..', 'Kelp_data', 'data', '2024_PROCESSED_DATA', mooring, 'L0', 'ADCP');
    fname = "ADCP_" + mooring + "_L0_10min.mat";
    M = load(fullfile(fpath,fname));
    %M.Config = load('../../../../Kelp_data/data/2024_PROCESSED_DATA/M1/L0/ADCP/ADCP_M1_config.mat');
    
    
    % apply qcFlag
    fprintf('applying qcFlag...\n')
    M.qcFlag(isnan(M.qcFlag)) = 0;
    qcFlag = round(M.qcFlag); % omg the qcFlag got smoothed
    fields = fieldnames(M);
    for i = 1:length(fields)
        field = fields{i};
        dum = M.(field);
        if size(dum) == size(qcFlag) & ~strcmp(field, 'qcFlag')
            dum(~qcFlag) = NaN;
            M.(field) = dum;
        end
    end

    %% Grid data
    card = {'East', 'North'};
        for d = 1:2
        direction = card{d};
        field = "Velocity_" + direction;
        
    
        
        % original reference
        hbar = nanmean(M.Pressure); % mean depth
        bin_dep = hbar - M.bin_mab; % each bin at mean depth
        bin_dep = [hbar; bin_dep]; % add bottom 0 bin
        o = M.(field)(1, :);%zeros(1, length(M.Time));
        U = [o;M.(field)]; % add 0 velocity at seafloor
        
        % reference grid
        [t_grid, z_grid] = meshgrid(M.Time, bin_dep);
        
        % clear nans
        mask = ~isnan(U);
        U = U(mask);
        t_grid = t_grid(mask);
        z_grid = z_grid(mask);
        
        % interpolate to a surface
        F = scatteredInterpolant(t_grid(:), z_grid(:), U(:), 'linear', 'nearest');
        
        % query grid
        z_max = max(bin_dep);
        z_min = hbar - min(M.Pressure);
        step = 0.5;
        dz = z_min+2:step:z_max;
        [Tq, Zq] = meshgrid(M.Time, dz);
        
        % Eval at query grid
        Vq.(field) = F(Tq, Zq);
        end

        %% Rotate to maj and minor axes

        
        fprintf('Rotating %s to Principal Axes\n', mooring)
        x = Vq.Velocity_East;
        y = Vq.Velocity_North;
        [U, V, R] = pca_rotation(x, y);
        M.U = U;
        M.V = V;
        

        card = {'U', 'V'};
        for d = 1:2
        direction = card{d};
        field = direction;
        % mean prfiles
        mean_profile(d).x = nanmean(M.(field), 2);
        std_profile(d).x = nanstd(M.(field),[], 2);
        
        figure
        subplot(1, 2, 1);
        plot(mean_profile(d).x, dz,'-s', 'LineWidth', 2)
        xlabel('$u$ [m/s]', 'Interpreter','latex')
        ylabel('$h$ [m]', 'Interpreter','latex')
        xline(0, 'k--', 'LineWidth', 1.5)
        set(gca, 'FontSize', 18)
        title('$\mu$ Profile', 'Interpreter','latex')
        grid minor
        axis ij
        
        subplot(1, 2, 2);
        plot(std_profile(d).x, dz, '-s', 'LineWidth', 2)
        xlabel('$u$ [m/s]', 'Interpreter','latex')
        ylabel('$h$ [m]', 'Interpreter','latex')
        set(gca, 'FontSize', 18)
        title('$\sigma$ Profile', 'Interpreter','latex')
        grid minor
        axis ij
        
        sgtitle(sprintf('%s %s Velocity', mooring, direction))
        
        end
        U_grid = M.U;
        V_grid = M.V;
        mean_profile_U = mean_profile(1).x;
        mean_profile_V = mean_profile(2).x;
        

        spath = fullfile(fpath, '..', '..', 'L1', 'ADCP');
        save(fullfile(spath, mooring + "_10min_gridded_PCA.mat"), 'Tq', 'Zq', 'U_grid', 'V_grid', 'dz', 'mean_profile_U', 'mean_profile_V', 'std_profile')


end