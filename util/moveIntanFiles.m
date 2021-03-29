D = uigetdir('Choose directory with files you want to move');
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % List of subfolders of D.
new_path = 'E:\ICMS18\Mar26_ICMS18_ELEC\allSingle';
mkdir(new_path)
% new_path = 'E:\ICMS15\Mar26_ICMS15_ELEC\pairedTest';
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*')); % Specify the file extension.
    C = {T(~[T.isdir]).name}; % Files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        movefile(F, new_path);
    end
end