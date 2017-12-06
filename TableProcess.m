function table = TableProcess(timeline, J, JC)
% TABLEPROCESS merges processes. Jumps in processes don't occur at the same
% times (EXP)
% IN : 
% timeline      vector 1xS          time linspace
% J             vector 1xM          new simulation jumps
% JC            matrix NxM          new simulation jump chain
% OUT : 
% tablenew      matrix (1+N)xS      merged simulations
% -------------------------------------------------------------------
% --- formal checks
if(numel(JC)==length(JC)) % we want JC to not be a vector, N = 1 : useless
    error('Input more than one particle, please.') 
end
if(timeline(1)<J(1)) % check that time is legit
    error('Cannot sample before process')
end

S = length(timeline);
N = size(JC, 1); 
                                                                             % --- define the new table
table = [J , timeline; ...
         JC, -1*ones(N, S)]; 

table = sortrows(table')';                                                  % sort along the times

                                                                            % --- look for repetitions in time 
reptime = table(1, :);                                                      % repeating time
noreptime = unique(table(1,:));                                             % non repeating time 
if(length(noreptime)<length(table(1,:)))                                    % get rid of repeating times
    compte = histc(reptime, noreptime);                                     % count time apparitions
    indx = compte~=1;                                                       % find repeating times
    quel = noreptime(indx);                                                 % find which times are repeating
    for i = quel
        table(:, reptime==i) = [max(table(:, reptime==i),[],2), NaN*ones(N+1, 1)]; % merge and NaN the extra column
    end
    table = table(:,all(~isnan(table)));                                    % delete the NaN columns
end

for i = 2: size(table, 1)                                                   % --- complete the missing (-1) times 
    for j = 2: size(table, 2)
        if(table(i, j)==-1)
            table(i, j) = table(i, j-1); 
        end
    end
end
                                                                            % --- get rid of time values that are not in timeline
indx = ismember(table(1,:), timeline);                                      % find the columns with time = timeline
table = table(:, indx);                                                     % keep only the associated columns
table = table(2:size(table,1), :);                                          % delete time from table
end