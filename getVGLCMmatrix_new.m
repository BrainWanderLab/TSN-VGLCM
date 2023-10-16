function [GLCM] = getVGLCMmatrix_new(ROI,levels,R,d)
levelTemp = max(levels)+1;ROI(isnan(ROI)) = levelTemp;
levels = [levels,levelTemp];
nL=max(levels);
%% offset: distance and orientation
direct = [0 d; -d d; -d 0; -d -d];%Angle: 0,45,90 135
%% X-Y plane
dim=size(ROI);
matrix_X=squeeze(ROI(dim(1)-R,:,:));
[GLCM(:,:,:,1),~] = graycomatrix(matrix_X,'NumLevels',nL,'Offset',direct,'G',[],'Symmetric',true);
matrix_Y=squeeze(ROI(:,dim(2)-R,:));
[GLCM(:,:,:,2),~] = graycomatrix(matrix_Y,'NumLevels',nL,'Offset',direct,'G',[],'Symmetric',true);
matrix_Z=squeeze(ROI(:,:,dim(2)-R));
[GLCM(:,:,:,3),~] = graycomatrix(matrix_Z,'NumLevels',nL,'Offset',direct,'G',[],'Symmetric',true);
GLCM=GLCM(1:end-1,1:end-1,:,:);
m=size(GLCM);
GLCM=reshape(GLCM,[m(1),m(2),m(3)*m(4)]);
