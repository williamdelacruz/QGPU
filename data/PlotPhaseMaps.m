close all
clc


% * * * * * * * * * * * * * * * * * *
%  Phase, correlation, mask,
%  unwrapped, residues, branch-cuts
% * * * * * * * * * * * * * * * * * *


filename = 'peaks512x512';
ysize = 512;
xsize = 512;
type = 'float';
 

% - - - - - - - - - - - - - -
%        Phase map
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.phase'));
if fid ~= -1
    data = fread(fid, type);
    fclose(fid);

    data = reshape(data, ysize, xsize);
    data1 = data;
    data = data/max(max(data));
    figure(1), 
     subplot(2,2,1), 
    imagesc(data), colormap(gray)
    title('Phase map')
    %set(gca,'FontSize',15)
else
    fprintf('File %s does not exists.\n', strcat(filename,'.phase'))
end


% - - - - - - - - - - - - - -
%     Correlation map
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.qual'));
if fid ~= -1
    data = fread(fid,'uint64');
    fclose(fid);

    data =reshape(data, ysize, xsize);
    data = data/max(max(data));
    figure(1), subplot(2,2,2), imagesc(data), colormap(gray)
    title('Quality map')
else
    fprintf('File %s does not exists.\n', strcat(filename,'.corr'))
end



% - - - - - - - - - - - - - -
%        Mask map
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.mask'));
if fid ~= -1
    data = fread(fid);
    fclose(fid);

    data = reshape(data, ysize, xsize);
    data = data/max(max(data));
    figure(1), subplot(2,2,3), imagesc(data), colormap(gray)
    title('Mask map')
else
    fprintf('File %s does not exists.\n', strcat(filename,'.mask'))
end


% - - - - - - - - - - - - - -
%    Unwrapped phase map
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.out'));
if fid ~= -1
    data = fread(fid,'float');
    fclose(fid);

    data = data - min(min(data));
    data = reshape(data, ysize, xsize);
    data = data/max(max(data));
    figure(1), subplot(2,2,4), imagesc(data), colormap(gray)
    title('Unwrapped phase map')
else
    fprintf('File %s does not exists.\n', filename)
end

