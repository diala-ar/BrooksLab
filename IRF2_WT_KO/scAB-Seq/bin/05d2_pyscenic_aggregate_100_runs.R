args = commandArgs(trailingOnly=TRUE)
cond = args[1]
occurrence_cutoff = as.numeric(args[2]) # percentage of TF occurrence in 100 runs.

library(parallel)

#############################################
# find high confidence regulons
#############################################

# save modules as list of TFs with their targets
extract_modules = function(reg) {
  TFs = reg$TF
  targets_str = reg$TargetGenes
  modules = mclapply(targets_str, function(x) {
    x = gsub("[(|[(|)|)]|'", '', x)
    x = strsplit(x, split=', ')[[1]]
    x = x[seq(1, length(x), 2)]
  }, mc.cores=10)
  names(modules) = TFs
  modules
}

# merge all targets of same TF in one element of a list
merge_modules = function(modules) {
  TFs = unique(names(modules))
  mod = mclapply(TFs, function(TF) unique(unlist(modules[names(modules)==TF])), mc.cores=10)
  names(mod) = TFs
  mod
}

# identify TFs occurring more than 80% of the time 
high_confidence_TFs = function(all_runs, cond) {
  TFs = lapply(1:length(all_runs), function(i) names(all_runs[[i]]))
  TF_freq = table(unlist(TFs))
  TF_freq_df = data.frame(TF_freq)
  names(TF_freq_df)[1] = 'TF'
  write.csv(TF_freq_df, file=paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/TF_freq_', cond, '.csv'), row.names=F)
  names(TF_freq)[TF_freq>=(length(all_runs)*occurrence_cutoff)]
}

# identify pairs of TF_targets occurring more than 80% of the time
high_confidence_modules = function(all_runs) {
  mod_df = NULL
  for (i in 1:length(all_runs)) {
    temp_mods = all_runs[[i]]
    temp_TFs  = names(temp_mods)
    mod_df[[i]] = do.call(rbind, lapply(1:length(temp_TFs), function(j) cbind(temp_TFs[j], temp_mods[[j]])))
    mod_df[[i]] = apply(mod_df[[i]], 1, function(x) paste(x, collapse='|'))
  }
  tf_target = unlist(mod_df)
  tf_target_freq = table(tf_target)
  tf_targets = names(tf_target_freq)[tf_target_freq>=(length(all_runs)*occurrence_cutoff)]
  tf_targets_df = matrix(unlist(strsplit(tf_targets, '|', fixed=T)), ncol=2, byrow=T)
  TFs = unique(tf_targets_df[, 1])
  hg_conf_mod = lapply(TFs, function(tf) tf_targets_df[tf_targets_df[,1]==tf, 2])
  names(hg_conf_mod) = TFs
  hg_conf_mod
}


redund_mod = merged_mod = NULL
for (i in 1:100) { #10:32
  reg = read.csv(paste0('results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/reg_', cond, '_run_', i,'.csv'), header=T, skip=1)
  names(reg)[1:2] = c('TF', 'MotifID')
  reg = reg[-1, ]
  
  temp_redund_mod = extract_modules(reg)
  merged_mod[[i]] = merge_modules(temp_redund_mod)
  if ((i %% 5)==0) print(paste(cond, i))
}

hg_conf_TFs = high_confidence_TFs(merged_mod, cond) # TFs occurring more than 80% of the times
mod_with_hg_conf_TFs = lapply(1:length(merged_mod), function(i) {
  temp = merged_mod[[i]][hg_conf_TFs]
  temp = temp[!sapply(temp, is.null)]
})
hg_conf_mods = high_confidence_modules(mod_with_hg_conf_TFs)


hg_conf_mods_output = lapply(1:length(hg_conf_mods), function(i) paste0(names(hg_conf_mods)[i], '\t', 'na\t', paste(hg_conf_mods[[i]], collapse = '\t'),'\t'))
sink(paste0("results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/pyscenic_100_runs/high_confidence_reg_", cond,".gmt"))
writeLines(unlist(lapply(hg_conf_mods_output, paste, collapse=" ")))
sink()


