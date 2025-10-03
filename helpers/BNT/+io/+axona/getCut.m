function clust = getCut(cutfile)
    fid = fopen(cutfile, 'r');
    if fid == -1
        error('Failed to open file %s', cutfile);
    end
    clust = zeros(1000000, 1);
    counter = 0;
    while ~feof(fid)
        string = strtrim(fgets(fid));
        if ~isempty(string)
            if string(1) == 'E'
                break;
            end
        end
    end
    while ~feof(fid)
      string = strtrim(fgets(fid));
      if ~isempty(string)
         content = sscanf(string, '%u')';
         N = length(content);

         clust(counter+1:counter+N) = content;
         counter = counter + N;
      end
    end
    fclose(fid);
    if counter > 0
        clust = clust(1:counter);
    else
        clust = [];
    end
end
