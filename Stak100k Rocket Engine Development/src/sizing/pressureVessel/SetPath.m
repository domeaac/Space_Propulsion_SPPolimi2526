function SetPath()

thisFile   = mfilename('fullpath');
rootFolder = fileparts(thisFile);

p = genpath(rootFolder);

exclude = {'.git','private','@','+','slprj'};

for k = 1:numel(exclude)
    p = regexprep(p, ...
        ['[^;]*' filesep exclude{k} '[^;]*;'], '');
end

addpath(p);

end
