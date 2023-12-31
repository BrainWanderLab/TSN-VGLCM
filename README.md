# TSN-VGLCM

1. Theory
 
  The structural covariance network (SCN) captures shared morphological covariance patterns between brain regions according to morphological measures derived from structural magnetic resonance imaging (sMRI) [1,2]. Texture similarity network (TSN) is the individual-level construction approach of the SCN based on the sMRI-derived voxel-wise 3D gray‐level co‐occurrence matrix (VGLCM) texture feature maps [3], which can provide detailed spatial information and dramatically diminishes the influence of insufficient areas definition on texture feature calculation. Compared with the tranditional SCN approach, the individualized TSN is independent of the sample size or priori hypotheses, which characterizes individual structural heterogeneity among the population using adequate VGLCM features.

In terms of TSN, the VGLCM was performed on the skull-stripped brain map to calculate each subject's texture maps in native space [4-6]. For a certain voxel, a GLCM matrix was reconstructed from the 5×5 voxel square region of interest (ROI) centering on this voxel on three orthogonal planes, respectively. All voxels within the ROI in each plane were scaled to 8 gray levels to reconcile VGLCM precision with computational complexity. Each cell in the GLCM matrix represented the number of co-occurrence times of any gray level pair between two neighboring voxels in four direction angles (0°, 45°, 90°, and 135°). Then, the GLCMs of the three orthogonal planes were averaged, and 20 texture metrics were estimated based on the average GLCM for each voxel. 

In addition, a multilevel 3D wavelet decomposition was performed on the original skull-stripped brain map, resulting in 8 wavelet-transformed brain maps from low-pass to high-pass components for each subject. And the above-mentioned VGLCM texture maps were repeatedly calculated for each wavelet-transformed map. Thus, 180 texture feature maps were ultimately generated from the sMRI data for each subject.

2. Usage
   
  This package processing as follows:

  Step 1: All sMRI data (T1-weighted brain image) were preprocessed using the CAT12 toolbox (http://dbm.neuro.uni-jena.de/cat/) implemented in SPM12 (http://www.fil.ion.ucl.ac.uk/spm/). And then, the native gray matter (e.g., p1*.nii) and white matter (e.g., p2*.nii) concentration maps were obtained for each subject using the CAT12 toolbox. In addtion, the inverted nonlinear deformation-field parameter (e.g., iy_*.nii) tranforming from the standard space to the individual native spce was remained for the following the individual-level network construction.

  Step 2: The standard whole brain atlas (e.g., BNA_maxprob_thr25_1.5mm.nii, 246 cerebral parcellations from Human Brainnetome Atlas) was pre-defined as the nodes of TSN for the following the individual-level network construction.

  Step 3: Voxel-wise GLCM texture maps were carried out using a self-developed program (e.g., batch_vox_glcm.sh) based on Matlab 2016b. A total of 180 texture feature maps was ultimately generated from the 9 typies of sMRI data (1 original and 8 wavelet-transformed brain maps) for each subject. This script contains 6 inputs as follows: (1) T1 rawdata map (corrected for bias-field); (2) the native gray matter map; (3) the native white matter map; (4) the result directory for saving texture feature maps; (5) Subject ID; (6) the full path of the package.
  
  E.g.,
  batch_vox_glcm.sh ./Subj001/mT1.nii ./Subj001/p1T1.nii ./Subj001/p2T1.nii ./result/Subj001 Subj001 ./TSN-VGLCM

  Step 4: TSN construction were also carried out using a self-developed program (e.g., batch_vox_glcm_matrix.sh) based on Matlab 2016b in Linux system. This script constructed the subject-level brain TSN using the pre-defined standard whole brain atlas from Step 2. The standard brain atlas was first warped into each subject's native space using the nonlinear deformation warp map (e.g., iy_*.nii) generated at Step 1. Then the feature vector of each parcellation of each subject was extracted from the 180 VGLCM maps. After that, a Pearson correlation was used to calculate the covariance coefficient of the feature vectors between each pair of areas. Then Fisher's r-to-z transformation algorithm was used to convert the covariance coefficient to approximately normally distributed, resulting in a symmetrical covariance matrix (termed TSN). This script contains 6 inputs as follows: (1) the directory for saving texture feature maps from Step 3; (2) the directory for saving T1 rawdata map and the nonlinear deformation warp map; (3) the full-path filename of the standard brain atlas; (4) the result directory for saving TSN; (5) Subject ID; (6) the full path of the package.
 
  E.g.,
  batch_vox_glcm_matrix.sh ./result/Subj001 ./Subj001 ./TSN-VGLCM/BNA_maxprob_thr25_1.5mm.nii Subj001 ./TSN-VGLCM
  
Reference:

[1] Pichet Binette A, Gonneaud J, Vogel JW, La Joie R, Rosa-Neto P, Collins DL, et al. (2020): Morphometric network differences in ageing versus Alzheimer's disease dementia. Brain: a journal of neurology. 143:635-649.

[2] Alexander-Bloch A, Giedd JN, Bullmore E (2013): Imaging structural covariance between human brain regions. Nat Rev Neurosci. 14:322-336.

[3] Ding H, Zhang Y, Xie YY, Du XT, Ji Y, Lin LY, Chang ZY, Zhang B, Liang M, Yu CS, Qin W: Individualized texture similarity network in schizophrenia (Under Review) 

[4] Maani R, Yang YH, Emery D, Kalra S (2016): Cerebral Degeneration in Amyotrophic Lateral Sclerosis Revealed by 3-Dimensional Texture Analysis. Frontiers in neuroscience. 10:120.

[5] Ta D, Khan M, Ishaque A, Seres P, Eurich D, Yang YH, et al. (2020): Reliability of 3D texture analysis: A multicenter MRI study of the brain. J Magn Reson Imaging. 51:1200-1209.

[6] Ishaque A, Mah D, Seres P, Luk C, Johnston W, Chenji S, et al. (2019): Corticospinal tract degeneration in ALS unmasked in T1-weighted images using texture analysis. Hum Brain Mapp. 40:1174-1183.
