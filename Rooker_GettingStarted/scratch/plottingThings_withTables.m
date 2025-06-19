

% Make Stats Table

cols = {'Avg_P', 'Std_P', 'Avg_U', 'Avg_theta'};
versions = ["L0", "M1", "Vector"];

T = array2table(zeros(3, 4), 'VariableNames', cols, 'RowNames', versions);

% Avg Tables for L0 and M1
for j = 1:3
    ind = sprintf('AvgR%d', j);
    for i = 1:length(versions)-1
        Stat = aquaplot(j, versions{i}, 'off');
        newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
        T{versions(i),:} = newRow;
        C = versions(i)+ind;
        T.Properties.RowNames{versions(i)} = C;
    end
    Table.(ind) = T;
end


% Tables for Stats at Vector Point
for j = 1:3
    ind = sprintf('PtR%d', j);
    
    %L0
    Stat = aquaplot(j, versions{1}, 'off', 34, 35);
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
    T{versions(1),:} = newRow;
    T.Properties.RowNames{versions(1)} = {versions(1) (ind)};
    %M1
    Stat = aquaplot(j, versions{2}, 'off', 11, 12);
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
    T{versions(2),:} = newRow;    
    T.Properties.RowNames{versions(2)} = {versions(2) (ind)};
    %Vector
    Stat = aquaplot(j, versions{3}, 'off');
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.theta];
    T{versions(3),:} = newRow;    
    T.Properties.RowNames{versions(3)} = {versions(3) (ind)};
    Table.(ind) = T;
end
Results = join(Table.AvgR3, join(Table.AvgR2, join(Table.AvgR2, join(Table.PtR3, join(Table.PtR2, Table.PtR1)))));
writetable(Results, 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\tables\Stats.csv')
