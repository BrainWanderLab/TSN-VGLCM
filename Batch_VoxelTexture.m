function Batch_VoxelTexture(datafile,GM_maskfile,WM_maskfile,resultdir)
type = 'GLCM';
Ng = 8;
R = 2;
d = 1;
del = 1;
% Orig
filter = 0;
VoxelTexture_Batch_individual(datafile,GM_maskfile,WM_maskfile,resultdir,type,filter,Ng,R,d,del);
% Wavelet
filter = 1;
VoxelTexture_Batch_individual(datafile,GM_maskfile,WM_maskfile,resultdir,type,filter,Ng,R,d,del);

end
%% main code
function VoxelTexture_Batch_individual(datafile,GM_maskfile,WM_maskfile,resultdir,type,filter,Ng,R,d,del)
%--------------------------------------------------------------------------
%% Parameter
%--------------------------------------------------------------------------
[~,filename,ext]=fileparts(datafile);
% Level
levels=[1:1:Ng];
%% prepare volume
V = spm_vol(datafile);
img = spm_read_vols(V);
img(isnan(img)) = 0;
img(isinf(img)) = 0;
% GM-mask
Vmask = spm_vol(GM_maskfile);
GM_imgmask = spm_read_vols(Vmask);
GM_imgmask(isnan(GM_imgmask)) = 0;
GM_imgmask(isinf(GM_imgmask)) = 0;
% WM-mask
Vmask = spm_vol(WM_maskfile);
WM_imgmask = spm_read_vols(Vmask);
WM_imgmask(isnan(WM_imgmask)) = 0;
WM_imgmask(isinf(WM_imgmask)) = 0;
% imgmask
imgmask = GM_imgmask+WM_imgmask;
imgmask(imgmask<0.5) = 0;
imgmask(imgmask>0) = 1;
maskfilename = [resultdir,filesep,'individual_mask.mat'];
if exist(maskfilename)~=2
    save(maskfilename,'imgmask','-v7.3');
end
imgmask = repmat(imgmask,[1,1,1,size(img,4)]);
img = img.*imgmask;
%% filter
if filter
    disp('...Doing Wavelet');
    imgstruct = waveletfilter(img);
    for i = 1:size(imgstruct,2)
        img(:,:,:,i) = imgstruct{1,i};  
    end
else
    disp('...Doing Orig');
    img(:,:,:,1) = img;
end
clear imgstruct i
% img
img = abs(img);
img(~imgmask) = NaN;
%% calculation type
% struct_texture = cell(length(find(imgmask==1)),5);
dim = size(imgmask);
for q = 1:size(img,4)
    if filter
         switch q
            case 1
                outdir = [resultdir,filesep,'text_wave_LLL'];
            case 2
                outdir = [resultdir,filesep,'text_wave_HLL'];
            case 3
                outdir = [resultdir,filesep,'text_wave_LHL'];
            case 4
                outdir = [resultdir,filesep,'text_wave_HHL'];
            case 5
                outdir = [resultdir,filesep,'text_wave_LLH'];
            case 6
                outdir = [resultdir,filesep,'text_wave_HLH'];
            case 7
                outdir = [resultdir,filesep,'text_wave_LHH'];
            case 8
                outdir = [resultdir,filesep,'text_wave_HHH'];
         end
    else
            outdir = [resultdir,filesep,'text_Orig'];
    end
    %
    if ~exist(outdir)
        mkdir(outdir);
    end
    % Start    
    tmpdata = squeeze(img(:,:,:,q));  
    num = 0;
    knum = 0;
    ll = 0;
    coords=[];
    features=[];
   % t1=clock;
   % tic
    for i = 1+R:dim(1)-R
        for j = 1+R:dim(2)-R
            for k = 1+R:dim(3)-R
                if imgmask(i,j,k) == 1
                    num = num+1;
                    coords(num,:)=[i,j,k];
                    I_range = [i-R:i+R];
                    J_range = [j-R:j+R];
                    K_range = [k-R:k+R];
                    ROI= tmpdata(I_range,J_range,K_range);
                    ROIvec=ROI(:);
                    ROIvec(isnan(ROIvec))=[];
                    ROImin=min(ROIvec);
                    ROImax=max(ROIvec);
                    ROI=(Ng-1)*(ROI-ROImin)./(ROImax-ROImin);
                    ROI=round(ROI)+1;
                    %% GLCM
                    switch upper(type)
                        case 'GLCM'
                            matrix = getVGLCMmatrix_new(ROI,levels,R,d);
                            glcm = GLCM_Features4(matrix);
                            if num==1
                                features_name=fieldnames(glcm);
                            end
                            tmp_feature=struct2array(glcm);                        
                            knum = knum+1;
                            features(knum,:) = TC_nanmean(tmp_feature, 1);
                            %features(knum,:)=mean(tmp_feature);
                            if knum ==10000
                                ll = ll+1;
                                save([outdir,filesep,'tmp_features_',num2str(ll),'.mat'],'features','-v7.3');
                                features = [];
                                knum = 0;
                            end
%                             if mod(num,1000)==0
%                                 t2=clock;
%                                 du=t2(5)*60+t2(6)-(t1(5)*60+t1(6));
%                                 fprintf('time used %0.0f sec/1000 voxels\n',du);
%                                 t1=t2;
%                             end
                    end
                end
            end
        end
    end
    ll = ll+1;
    save([outdir,filesep,'tmp_features_',num2str(ll),'.mat'],'features','-v7.3');
   % toc
    clear glcm i j k I_range J_range K_range knum matrix num ROI ROImax ROImin ROIvec tmp_feature
    coordvox=coord2scalar(coords,dim);
    clear coords
    features = [];
    for i = 1:ll
        Y = load([outdir,filesep,'tmp_features_',num2str(i),'.mat']);
        features = [features;Y.features];
        clear Y
    end
    tmpimg=zeros(dim);
    for f=1:length(features_name)
        outfile=[outdir,filesep,filename,'_',features_name{f},ext];
        data=reshape(tmpimg,[],1);
        data(coordvox)=features(:,f);
        data=reshape(data,dim(1),dim(2),dim(3));
        V.fname=outfile;
        V.dt=[64,0];
        spm_write_vol(V,data);
        clear data
    end
    clear tmpimg
    if del
        files = dir([outdir,filesep,'tmp_features_*.mat']);
        for f = 1:length(files)
            delete([outdir,filesep,files(f).name]);
        end
    end
    fprintf('.');
end
fprintf('\n');
end

function [volumestruct] = waveletfilter(volume)
lev = 1;% Decomposition Level
WT = wavedec3(volume, lev, 'sym4');
volumestruct = cell(1,8);
volumestruct{1} = waverec3(WT,'LLL');
volumestruct{2} = waverec3(WT,'HLL');
volumestruct{3} = waverec3(WT,'LHL');
volumestruct{4} = waverec3(WT,'HHL');
volumestruct{5} = waverec3(WT,'LLH');
volumestruct{6} = waverec3(WT,'HLH');
volumestruct{7} = waverec3(WT,'LHH');
volumestruct{8} = waverec3(WT,'HHH');
%fprintf('Done.\n');
end

function TC = TC_nanmean(x, dim)
if(nargin == 1)
    dim = 1;
    if(size(x, 1) == 1)
        dim = find(size(x) > 1, 1);
    end
end

% find nan entries
nanind = isnan(x);

% set nan entries to be zeros and sum
x(nanind) = 0;
xsum = sum(x, dim);

% count nan-nan entries
count = size(x, dim) - sum(nanind, dim);

% nanmean
TC = xsum ./ count;


end

