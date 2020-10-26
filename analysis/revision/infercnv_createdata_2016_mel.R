library(BoutrosLab.plotting.general);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/raw/melanoma_2016');
mycount <- read.table('GSE72056_melanoma_single_cell_revised_v2.txt', header = TRUE);
annot <- data.frame(t(mycount[1:3, ]));
colnames(annot) <- c('tumor', 'malignant', 'non');
annot <- annot[-1, ];
annot$malignant <- gsub(' ', '', annot$malignant)
annot$type <- ifelse(annot$malignant=='1', 'ref', 'others');
saveRDS(annot, file = generate.filename('melanoma_2016', 'annotation', 'rds'))
mycount <- mycount[-(1:3), ];
mycount <- mycount[!rownames(mycount)%in%c('8030', '23205'), ];
rownames(mycount) <- as.vector(mycount$Cell);
mycount <- mycount[, -1];
write.table(mycount, generate.filename('infercnv_count', 'melanoma_2016', 'txt'), quote = FALSE, col.names = NA, sep = '\t')
write.table(annot[, 'type', drop = FALSE], generate.filename('infercnv_type', 'melanoma_2016', 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');

library(infercnv);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/raw/melanoma_2016');
mycounti <- read.table('2020-03-10_infercnv_count_melanoma_2016.txt');
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= mycounti,
    annotations_file= '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/raw/melanoma_2016/2020-03-10_infercnv_type_melanoma_2016.txt',
    delim="\t",
    gene_order_file="/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt",
    ref_group_names=c('ref')
    );

#out_dir=paste0(Sys.Date(), "_infercnv_raw_ref_", name);
out_dir=paste0(Sys.Date(), "_infercnv_null", 'melanoma_2016');
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=out_dir,
    cluster_by_groups=T,
    mask_nonDE_genes = T,
    include.spike = T
    );
save(infercnv_obj, file = paste0(out_dir, '/', 'infercnv_obj_output.rda'));
####
mycount <- mycounti
mycnv <- read.table('2020-03-10_infercnv_nullmelanoma_2016/infercnv.observations.txt', row.names = 1, as.is = TRUE, header=TRUE);
mycnv.ref <- read.table('2020-03-10_infercnv_nullmelanoma_2016/infercnv.references.txt', row.names = 1, as.is = TRUE, header=TRUE);
for(i in unique(annot$tumor)){
	mycounti <- mycount[, rownames(annot[annot$tumor==i, ])];
	itype <- annot[annot$tumor==i, 'type', drop = FALSE];
	write.table(mycounti, generate.filename('infercnv_count', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
	write.table(itype, generate.filename('infercnv_type', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};
####
annot <- readRDS('2020-03-10_melanoma_2016_annotation.rds');



ss=`ls 2020-03-10_infercnv_type_*txt|grep -v melanoma_2016`
for s in $ss
do
        name=`echo $s|cut -f 1 -d .|sed 's/.*_//g'`
        sbatch /cluster/projects/hansengroup/sujunc/scRNA/script/revision/infercnv_2016_per.sh $name
done

#!/bin/bash
#SBATCH --mem=30G
#SBATCH -J run_cnv_raw
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o /cluster/projects/hansengroup/sujunc/scRNA/log/%x-%j.out
module load R/3.5.0
cd /cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/raw/v0.8.2/
Rscript /cluster/projects/hansengroup/sujunc/scRNA/script/revision/infercnv_2016_per.R $1
