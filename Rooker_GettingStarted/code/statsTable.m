

% Make Stats Table

cols = {'Avg_P', 'Std_P', 'Avg_U','Std_U','Avg_E','Std_E','Avg_N','Std_N','Avg_theta'};
versions = {'L0', 'M1', 'Vector'};



% Avg Tables for L0 and M1
for j = 1:3
    ind = sprintf('AvgR%d', j);
    clear T
   T = array2table(zeros(2, length(cols)), 'VariableNames', cols, 'RowNames', versions(1:2));
    
        Stat = aquaplot(j, versions{1}, 'off');
        newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.Mag(2), Stat.East(1), Stat.East(2), Stat.North(1), Stat.North(2), Stat.theta];
        T{versions(1),:} = newRow;
        T.Properties.RowNames{1} = [ind ' ' versions{1}];
    figure, quiver(0, 0, Stat.East(1), Stat.North(1), 0, 'b')
    hold on
        Stat = aquaplot(j, versions{2}, 'off', 1, 12);
        newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.Mag(2), Stat.East(1), Stat.East(2), Stat.North(1), Stat.North(2), Stat.theta];
        T{versions(2),:} = newRow;
        T.Properties.RowNames{2} = [ind ' ' versions{2}];
        quiver(0, 0, Stat.East(1), Stat.North(1), 0, 'r')
    hold off
    
    Table.(ind) = T;
end




% Tables for Stats at Vector Point
for j = 1:3
    ind = sprintf('PtR%d', j);
    clear T
    T = array2table(zeros(3, length(cols)), 'VariableNames', cols, 'RowNames', versions);
    %L0
    Stat = aquaplot(j, versions{1}, 'off', 34, 35);
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.Mag(2), Stat.East(1), Stat.East(2), Stat.North(1), Stat.North(2), Stat.theta];
    T{versions(1),:} = newRow;
    T.Properties.RowNames{1} = [ind ' ' versions{1}];
    figure, quiver(0, 0, Stat.East(1), Stat.North(1), 0, 'b')
    hold on
    %M1
    Stat = aquaplot(j, versions{2}, 'off', 5, 6);
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.Mag(2), Stat.East(1), Stat.East(2), Stat.North(1), Stat.North(2), Stat.theta];
    T{versions(2),:} = newRow;    
    T.Properties.RowNames{2} = [ind ' ' versions{2}];
    quiver(0, 0, Stat.East(1), Stat.North(1), 0, 'r')
    %Vector
    Stat = aquaplot(j, versions{3}, 'off');
    newRow = [Stat.Pres(1), Stat.Pres(2), Stat.Mag(1), Stat.Mag(2), Stat.East(1), Stat.East(2), Stat.North(1), Stat.North(2), Stat.theta];
    T{versions(3),:} = newRow;    
    T.Properties.RowNames{3} = [ind ' ' versions{3}];
    Table.(ind) = T;
    quiver(0, 0, Stat.East(1), Stat.North(1), 0, 'g')
    hold off
end




Results = [Table.AvgR1; Table.AvgR2; Table.AvgR3; Table.PtR1; Table.PtR2; Table.PtR3]; 
writetable(Results, 'C:\Users\jkr6136\OneDrive - UNC-Wilmington\Kelp_data\Summer2025\Rooker\tables\Stats.xlsx', 'WriteRowNames', true);



