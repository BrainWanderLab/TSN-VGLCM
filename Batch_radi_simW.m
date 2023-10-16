function Batch_radi_simW(outdir,vglcm_datadir,mT1datadir,subj,MNI_datafile)
disp('...MNI_2_Individual_space');
[~,roiname,~] = fileparts(MNI_datafile);
mT1_datafile = [mT1datadir,filesep,subj,filesep,'miT1.nii'];
indiviudal_defmap_file = [mT1datadir,filesep,subj,filesep,'iy_T1.nii'];
subj_vglcm_datadir = [vglcm_datadir,filesep,subj];
MNI_2_Individual_space(subj_vglcm_datadir,MNI_datafile,mT1_datafile,indiviudal_defmap_file);
disp('...texture_normalise');
[data,idmask] = texture_normalise(vglcm_datadir,subj);
disp('...texture_similarity');
subj_ROImaskfile = [subj_vglcm_datadir,filesep,'individual_',roiname,'.nii'];
texture_ROI = texture_extract(data,idmask,subj_ROImaskfile);
CorrMat = corr(texture_ROI);
CorrMat(isnan(CorrMat)) = 0;
CorrMat(isinf(CorrMat)) = 0;
CorrMat(CorrMat>= 1) = 1 - 1E-16; %. Supress the voxels with extremely high correlation values .
% Fisher's r-z transform
zCorrMat = (0.5 * log((1 + CorrMat)./(1 - CorrMat)));
zCorrMat(1:size(zCorrMat,1)+1:end) = 0;
save([outdir,filesep,'Sim_glcm_',subj,'.txt'],'CorrMat','-ascii', '-tabs');
save([outdir,filesep,'zSim_glcm_',subj,'.txt'],'zCorrMat','-ascii', '-tabs');
end
%%
function [texture_ROI] = texture_extract(data,idmask,ROImaskfile)
% data: vox*num_sum_features
%% mask
Vmask = spm_vol(ROImaskfile);
ROImask = spm_read_vols(Vmask);
ROImask(isnan(ROImask)) = 0;
ROImask(isinf(ROImask)) = 0;
tempmask = reshape(zeros(Vmask.dim),[],1);
tempmask(idmask,1) = 1;
tempmask = reshape(tempmask,Vmask.dim(1),Vmask.dim(2),Vmask.dim(3));
ROImask = round(ROImask).*tempmask;
label = unique(reshape(ROImask,[],1));
label(label==0) = [];
texture_ROI = zeros(size(data,2),length(label));
for r = 1:length(label)
    tmp = find(ROImask==label(r));
    row = [];
    for ff = 1:length(tmp)
        row(ff,1) = find(idmask == tmp(ff));
    end
    texture_ROI(:,r) = mean(data(row,:))';
end
end
%%
function [data_all,idmask] = texture_normalise(datadir,subj)
% data_all: voxel*num_sum_features
Y = load([datadir,filesep,subj,filesep,'individual_mask.mat']);
imgmask = Y.imgmask;clear Y
idmask = find(imgmask==1);
data_all = zeros(length(find(imgmask==1)),22,9);
%% data
files = dir([datadir,filesep,subj,filesep,'text_*']);
for f = 1:length(files)
    sub_tmpdir = [datadir,filesep,subj,filesep,files(f).name];
    data = [];
    img_files = dir([sub_tmpdir,filesep,'miT1_*.nii.gz']);
    for i = 1:length(img_files)
       % [~,fname,~] = fileparts([sub_tmpdir,filesep,img_files(i).name]);
        V = spm_vol([sub_tmpdir,filesep,img_files(i).name]);
        img = spm_read_vols(V);
        img(isnan(img)) = 0;
        img(isinf(img)) = 0;
        imgvec = img(find(imgmask==1));
        imgvec = (imgvec-mean(imgvec))./std(imgvec); % Z-trans
        data = [data imgvec]; 
        %
       % temp = zeros(V.dim);
      %  temp(find(imgmask==1)) = imgvec;
      %  V.fname = [sub_tmpdir,filesep,'zscore_',fname];
     %   spm_write_vol(V,temp);   
    end
    data_all(:,:,f) = data;
end
data_all = reshape(data_all,size(data_all,1),size(data_all,2)*size(data_all,3))';
% Z-trans
% data_all = (data_all-repmat(mean(data_all),size(data_all,1),1))./repmat(std(data_all),size(data_all,1),1);
data_all = data_all';
end