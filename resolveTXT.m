function [metabolites, transport, reactions] = resolveTXT()
    path = 'inputs/InputFile.txt';
    lines = {};
    metabolites = {};
    transport = {};
    reactions = {};
    fp = fopen(path, 'r');

    line = fgetl(fp);

    while ischar(line)

        if ~startsWith(line, '#') && strlength(line) > 0
            lines{end + 1} = line;
        end

        line = fgetl(fp);
    end

    [~, mpos] = ismember('METABOLITES', lines);
    [~, tpos] = ismember('TRANSPORT', lines);
    [~, rpos] = ismember('REACTIONS', lines);

    metabolites = lines(mpos + 1:tpos - 1);

    cur = 0;

    for i = tpos + 1:rpos - 1
        [in, id] = ismember(lines{i}, metabolites);

        if in
            transport(end + 1, :) = {lines(i), containers.Map};
            cur = id;
        else
            pair = strip(split(lines(i), ':'));

            if str2double(pair{2}) ~= 0
                transport{end, 2}(pair{1}) = str2double(pair{2});
            else
                transport{end, 2}(pair{1}) = pair{2};
            end

        end

    end

    % format for each reaction :
    name_vec = {'r_vec', 'p_vec', 'nr_vec', 'np_vec', 'kp', 'km', 'Kr_vec', 'Kp_vec', 'DGro', 'ap', 'am', 'aeq'};
    input_size = 10; % no km and am
    name_map = containers.Map(name_vec, num2cell(1:length(name_vec)));
    reactions = cell((length(lines) - rpos) / input_size, length(name_vec));

    for i = 0:(length(lines) - rpos) / input_size - 1

        for j = rpos + 1 + i * input_size:rpos + input_size + i * input_size
            pair = strip(split(lines(j), ':'));

            if endsWith(pair{1}, 'vec')
                
                pair{2} = strip(split(pair{2},','))';
                if ~isnan(sum(str2double(pair{2})))
                    pair{2} = str2double(pair{2});
                end
                

            elseif ~isnan(str2double(pair{2}))
                pair{2} = str2double(pair{2});
            end
            reactions{i+1, name_map(pair{1})} = pair{2};


            % if str2double(pair{2}) ~= []
            %     reactions{i+1, name_map(pair{1})} = str2double(pair{2});
            % else
            %     reactions{i+1, name_map(pair{1})} = pair{2};
            % end
        end

        % calculate the km and am
        reactions{i + 1, name_map('km')} = reactions{i + 1, name_map('kp')} / Keq(reactions{i + 1, name_map('DGro')}, 300);
        reactions{i + 1, name_map('am')} = reactions{i + 1, name_map('ap')} / reactions{i + 1, name_map('aeq')};
    end
    

end

function output = Keq(DGro, T)
    R = 8.314; %J/(mol*K)
    output = exp(-DGro /( R * T));
end
