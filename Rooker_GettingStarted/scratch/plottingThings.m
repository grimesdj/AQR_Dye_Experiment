

% Make Stats Table

cols = {'Avg_P', 'Std_P', 'Avg_U', 'Avg_theta'};
versions = {'L0', 'M1', 'Vector'};

T = array2table(zeros(2, 4), 'VariableNames', cols, 'RowNames', versions(1:2));

% Avg Tables for L0 and M1
for j = 1:3
    ind = sprintf('AvgR%d', j);
    clear T
    T = array2table(zeros(2, 4), 'VariableNames', cols, 'RowNames', versions(1:2));
    for i = 1:length(versions)-1
        Stat = aquaplot(j, versions{i}, 'off');
        newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
        T{versions(i),:} = newRow;
        T.Properties.RowNames{i} = [ind ' ' versions{i}];
    end
    Table.(ind) = T;
end




% Tables for Stats at Vector Point
for j = 1:3
    ind = sprintf('PtR%d', j);
    clear T
    T = array2table(zeros(3, 4), 'VariableNames', cols, 'RowNames', versions);
    %L0
    Stat = aquaplot(j, versions{1}, 'off', 34, 35);
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
    T{versions(1),:} = newRow;
    T.Properties.RowNames{1} = [ind ' ' versions{1}];
    %M1
    Stat = aquaplot(j, versions{2}, 'off', 11, 12);
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
    T{versions(2),:} = newRow;    
    T.Properties.RowNames{2} = [ind ' ' versions{2}];
    %Vector
    Stat = aquaplot(j, versions{3}, 'off');
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
    T{versions(3),:} = newRow;    
    T.Properties.RowNames{3} = [ind ' ' versions{3}];
    Table.(ind) = T;
end


Results = [Table.AvgR1; Table.AvgR2; Table.AvgR3; Table.PtR1; Table.PtR2; Table.PtR3]; 
writetable(Results, 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\tables\Stats.xlsx', 'WriteRowNames', true);



