function MNI_2_Individual_space(outdir,MNI_datafile,mT1_datafile,indiviudal_defmap_file)
cfg_batch{1}.spm.util.defs.comp{1}.def = {indiviudal_defmap_file};
cfg_batch{1}.spm.util.defs.comp{2}.id.space = {mT1_datafile};
cfg_batch{1}.spm.util.defs.out{1}.pull.fnames = {MNI_datafile};
cfg_batch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {outdir};
cfg_batch{1}.spm.util.defs.out{1}.pull.interp = 0;
cfg_batch{1}.spm.util.defs.out{1}.pull.mask = 1;
cfg_batch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
cfg_batch{1}.spm.util.defs.out{1}.pull.prefix = 'individual_';
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run',cfg_batch);
end