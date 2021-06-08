# VIKAS PEDDU 5/27/2021
# BME237 
library("sleuth")
library("biomaRt")
library('tximport')
library('ggplot2')

setwd('/Users/vikas/Documents/UCSC/classes/BME237/final_project')
sample_id <- dir(file.path("/Users/vikas/Documents/UCSC/classes/BME237/final_project/kallisto_results_real/"))

s2c <- read.table(file.path("sleuth_metadata.txt"), header = TRUE, stringsAsFactors=FALSE)
#s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c$path<-NA
for(i in 1:length(sample_id)){ 
  s2c$path[i]<-file.path(getwd(),'kallisto_results_real',sample_id[i])
}

# Human gene annotation from ensembl to gene symbol
human_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
human_t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"),
                            #filters = grep('refseq',  listFilters(human_mart)$name, ignore.case=T, value=T),
                            mart = human_mart)
human_t2g <- dplyr::rename(human_t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
human_t2g$ext_gene[which(human_t2g$ext_gene == '')] <- human_t2g$ens_gene[which(human_t2g$ext_gene == '')]

# Chimp gene annotation from ensembl to gene symbol
chimp_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "ptroglodytes_gene_ensembl",
                         host = 'ensembl.org')
chimp_t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                           "external_gene_name"), mart = chimp_mart)
chimp_t2g <- dplyr::rename(chimp_t2g, target_id = ensembl_transcript_id,
                           ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
chimp_t2g$ext_gene[which(chimp_t2g$ext_gene == '')] <- chimp_t2g$ens_gene[which(chimp_t2g$ext_gene == '')]

# Orangutan gene annotation from ensembl to gene symbol
orang_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                               dataset = "pabelii_gene_ensembl",
                               host = 'ensembl.org')
orang_t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                           "external_gene_name"), mart = orang_mart)
orang_t2g <- dplyr::rename(orang_t2g, target_id = ensembl_transcript_id,
                           ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
orang_t2g$ext_gene[which(orang_t2g$ext_gene == '')] <- orang_t2g$ens_gene[which(orang_t2g$ext_gene == '')]



runsleuth<-function(df, species, compgroup,t2g){ 
  specdf<-df[(grepl(species,df$sample)) & (df$condition ==  compgroup[1] | df$condition == compgroup[2]),]
          
  #print(df[which(grepl(as.character(species)))])
  #print(specdf)
  outobj <- sleuth_prep(specdf, extra_bootstrap_summary = TRUE,dropInfReps=TRUE, num_cores = 8, target_mapping = t2g, gene_mode = TRUE, aggregation_column = 'ext_gene')
  outobj <- sleuth_fit(outobj, ~condition, 'full')
  outobj <- sleuth_fit(outobj, ~1, 'reduced')
  outobj <- sleuth_lrt(outobj, 'reduced', 'full')

  return(outobj)
  }
cat('running human')
so_cm_human<-runsleuth(s2c,'hips',c('cyto','mono'), human_t2g)
so_cp_human<-runsleuth(s2c,'hips',c('cyto','polh'), human_t2g)

cat('running chimp')
so_cm_chimp<-runsleuth(s2c,'cips',c('cyto','mono'), human_t2g)
so_cp_chimp<-runsleuth(s2c,'cips',c('cyto','polh'), human_t2g)

cat('running orang')
so_cm_orang<-runsleuth(s2c,'oips',c('cyto','mono'), human_t2g)
so_cp_orang<-runsleuth(s2c,'oips',c('cyto','polh'), human_t2g)



# q is the FDR 
bind_and_fiter<-function( outdf,name,comp,cond){ 
  wt<-sleuth_wt(comp, 
            which_beta = cond)
  indf<-sleuth_results(wt, test = cond, show_all = FALSE)
  indf$sample<-as.character(name)
  tempfilter<-dplyr::filter(indf, qval <= 0.01)
  tempout<-rbind(outdf,tempfilter)
  return(tempout)
}
wt<-sleuth_wt(so_cm_human, which_beta = 'conditionmono')
sleuth_table <- sleuth_results(wt, test = 'conditionmono')
sleuth_table$sample<-'cm_human'
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.01)

sleuth_significant<-bind_and_fiter(sleuth_significant,'cp_human',so_cp_human,'conditionpolh')
sleuth_significant<-bind_and_fiter(sleuth_significant,'cm_chimp',so_cm_chimp,'conditionmono')
sleuth_significant<-bind_and_fiter(sleuth_significant,'cp_chimp',so_cp_chimp,'conditionpolh')
sleuth_significant<-bind_and_fiter(sleuth_significant,'cm_orang',so_cm_orang,'conditionmono')
sleuth_significant<-bind_and_fiter(sleuth_significant,'cp_orang',so_cp_orang,'conditionpolh')
sleuth_significant<-sleuth_significant[complete.cases(sleuth_significant$ext_gene),]

save.image(file = "sleuth_analysis.rdata")
sleuth_significant<-sleuth_significant[which(sleuth_significant$sample[sleuth_significant$b > 2 | sleuth_significant$b < -2,])

write.csv(sleuth_significant, 'sleuth_p05_significant.csv')

ggplot(sleuth_significant, aes(x = b, y = -log10(qval), color = sample)) + 
  geom_point() + 
  theme_minimal() + 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#tximport aggregation
test = getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", mart= chimp_mart)

x<-tximport('/Users/vikas/Documents/UCSC/classes/BME237/final_project/kallisto_results_real/hips_polh_jd957_1.fastq.gz/abundance.h5',
         type = 'kallisto',
         tx2gene = human_t2g,
         ignoreTxVersion = TRUE)
