function showMat(varargin)
% usage: showMat(mat1,[mat2 [,mat3 [,mat4]]]), where [] means optional
% 1nd parameter expects img to be shown in image slices
% 2st parameter expects mat to be shown in iso-lines
% 3rd parameter expects voxel sizes (3*1 or 1*3 matrix), unit cm, default
% (0.1; 0.1; 0.1) cm
% if 3nd or 4nd parameter is set to 0, showMat(...) returns immediately
% instead of holding the GUI
%
% for instance, dose, img are two 3D matrix repenting dose and CT images
% showMat(img) only displays the CT images
% showMat(0, dose) only displays the dose in iso-lines
% showMat(img, dose) displays both of them.
% showMat(img, 0, [0.1, 0.1, 0.15]) displays the CT with given voxel sizes
% showMat(0, dose, 0) only displays the dose in a detached process.

if numel(varargin) == 0
    fprintf('Need at least one parameter in function showMat\n');
    return;
end

if numel(varargin)>=1 && ~isequal(varargin{1},0) hasImg= 1;
else hasImg = 0;
end

if numel(varargin) >=2 && ~isequal(varargin{2},0) hasMat = 1;
else hasMat = 0;
end

vs=[0.1;0.1;0.1];
keepGUI = 1;
if numel(varargin) >=3 
    if numel(varargin{3})==3
        vs=[varargin{3}(1);varargin{3}(2);varargin{3}(3)];
    end
    if isequal(varargin{3},0)
        keepGUI = 0;
    end
end

if numel(varargin) >=4 && isequal(varargin{4},0)
    keepGUI = 0;
end


if hasMat ==0 && hasImg ==0
    fprintf('cannot have mat and img empty at the same time\n');
    return;
else
    if hasImg
        sz = size(varargin{1});
    end
    if hasMat
        sz = size(varargin{2});
    end
end
% write the data
if ispc
    fpath = [getenv('userprofile'),'\temp.mats'];
else
    fpath = [getenv('HOME'),'/temp.mats'];
end
fid=fopen(fpath,'wb');
if fid < 0
    warning('Cannot open the file temp.mats to write.\n');
    return;
end
fwrite(fid,uint32(sz(1)),'uint32');
fwrite(fid,vs(1),'double');
fwrite(fid,uint32(sz(2)),'uint32');
fwrite(fid,vs(2),'double');
fwrite(fid,uint32(sz(3)),'uint32');
fwrite(fid,vs(3),'double');
fwrite(fid,uint32(hasMat),'uint32');
if hasMat
    fwrite(fid,varargin{2},'double');
end
fwrite(fid,uint32(hasImg),'uint32');
if hasImg
    fwrite(fid,varargin{1},'double');
end
fclose(fid);
fpath = ['"', fpath,'"']; % in case there's space character in the path
if keepGUI
    system(['3DMatViewer ', fpath]);
else
    system(['3DMatViewer ', fpath, ' &']);
end
