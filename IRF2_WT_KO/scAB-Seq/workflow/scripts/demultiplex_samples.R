library(data.table)

ncounts_file = snakemake@input[['ncounts_file']]
tags_file    = snakemake@input[['tags_file']]
metadata_file  = snakemake@params[['samples_metadata']]

raw_counts = as.data.frame(fread(ncounts_file, skip=8))
tags       = as.data.frame(fread(tags_file))
metadata   = read.csv(metadata_file)

# change antidody derived tags (adt) long names with shorter ones
adt_index = grep('pAbO$', colnames(raw_counts))
adt = colnames(raw_counts)[adt_index]
adt = gsub(':.*$', '-adt', adt) # replace all characters after the ':' by '-adt'
adt = gsub('\\|.*$', '-adt', adt) # replace all characters after the '|' by '-adt'
colnames(raw_counts)[adt_index] = adt

tagged_raw_counts = merge(raw_counts, tags)
tagged_raw_counts = tagged_raw_counts[!tagged_raw_counts$Sample_Name %in% c('Multiplet', 'Undetermined'), ]
table(tagged_raw_counts$Sample_Name)

# rename samples according to content of metadata_file
for (i in 1:nrow(metadata)) {
  tagged_raw_counts[tagged_raw_counts$Sample_Name==metadata$sample[i], 'Cell_Index'] = 
    paste0(metadata$condition[i],'-', tagged_raw_counts[tagged_raw_counts$Sample_Name==metadata$sample[i], 'Cell_Index']) 
}
head(tagged_raw_counts[, 1:5])
tagged_raw_counts = tagged_raw_counts[, -c(ncol(tagged_raw_counts), ncol(tagged_raw_counts)-1)]  # remove Sample_Tag, Sample_Name columns

# personalize
fwrite(tagged_raw_counts, file=snakemake@output[['tagged_raw_counts']], row.names=F)


