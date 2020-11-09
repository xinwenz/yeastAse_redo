# This is file xinwen@YAYaNUCUbintu:~/cloud/yeastAse_redo/projectGuide.txt
Goal: Assemble genomes
	Code: HPC ?
	Input: ?
	Output:	/home/xinwen/mnt/nnp/0_3_0_pilon_busco/rm11_B.fasta   
			/home/xinwen/mnt/nnp/0_3_0_pilon_busco/yps128_5.fasta

Goal: Find all snps between rm11_B and yps128_5
	Code: [HPC] run mummer delta, show-snps, code?????
	Input:	/home/xinwen/mnt/nnp/0_3_0_pilon_busco/rm11_B.fasta   
			/home/xinwen/mnt/nnp/0_3_0_pilon_busco/yps128_5.fasta
	Output:	/home/xinwen/mnt/nnp/1_1_0_cali_bowtie/0_5Bsnp/ryB5.posmap
			/home/xinwen/mnt/nnp/1_1_0_cali_bowtie/0_5Bsnp/yr5B.posmap
		
		
Goal: Keep intersect SNP pair
		# in ry compare and yr compare, keep intersect SNP pair
		# labe snps in one chrom into a group
		# find yr-map, condense indels into one position
	Code:/home/xinwen/cloud/yeastAse_redo/Scripts/Z1_snp_nofilter.R
	Input:	/home/xinwen/mnt/nnp/1_1_0_cali_bowtie/0_5Bsnp/ryB5.posmap
			/home/xinwen/mnt/nnp/1_1_0_cali_bowtie/0_5Bsnp/yr5B.posmap
	Output:	/home/xinwen/cloud/yeastAse_redo/Z1_out/posmap52.RData
			/home/xinwen/cloud/yeastAse_redo/Z1_out/rm11_B_snpls_nofilter
			/home/xinwen/cloud/yeastAse_redo/Z1_out/yps128_5_snpls_nofilter

Goal: Hisat2 DNA mapping, collct snps counts in snp list
	Code: HPC ? 
	Input: 	/home/xinwen/mnt/nnp/1_0_1_calibrate/deMx_trim/D_*_[AH]_r1.fq.gz
			/home/xinwen/mnt/nnp/1_0_1_calibrate/deMx_trim/D_*_[AH]_r2.fq.gz

			/home/xinwen/cloud/yeastAse_redo/Z1_out/rm11_B_snpls_nofilter
			/home/xinwen/cloud/yeastAse_redo/Z1_out/yps128_5_snpls_nofilter
	Output:
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_rm11_B/D_*_[AH]_nosp.bam
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_yps128_5/D_*_[AH]_nosp.bam

			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_rm11_B/D_*_[AH]_nosp.pileup.BPNAME
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_yps128_5/D_*_[AH]_nosp.pileup.BPNAME

			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_rm11_B/D_*_[AH]_nosp.pileup.BPNAME.snpCount
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_yps128_5/D_*_[AH]_nosp.pileup.BPNAME.snpCount
			 

Goal: Plot of sum of counts in each allele
		# Also Get cal.Rata, a summary of snp read counts
	Code:	/home/xinwen/cloud/yeastAse_redo/Scripts/Z1_hst_cal_sum.R
	Input: 	/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_rm11_B/D_*_[AH]_nosp.pileup.BPNAME.snpCount
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_yps128_5/D_*_[AH]_nosp.pileup.BPNAME.snpCount
	Output:	/home/xinwen/cloud/yeastAse_redo/Z1_out/allele_sum_nofilter.pdf
			/home/xinwen/cloud/yeastAse_redo/Z1_out/cal.RData

Goal: Filter out bad snps, get new snp list, also snp with group number
	Code:	/home/xinwen/cloud/yeastAse_redo/Scripts/Z1_filter_bad.R
	Input:	/home/xinwen/cloud/yeastAse_redo/Z1_out/cal.RData
	Output:	/home/xinwen/cloud/yeastAse_redo/Z1_out/cal_evd_p_allinfo.RData

			/home/xinwen/cloud/yeastAse_redo/Z1_out/evidence_p_snps_counts_sum.pdf
			/home/xinwen/cloud/yeastAse_redo/Z1_out/evidence_snps_counts_sum.pdf

			/home/xinwen/cloud/yeastAse_redo/Z1_out/rm11_B_snpls_evdp_group
			/home/xinwen/cloud/yeastAse_redo/Z1_out/yps128_5_snpls_evdp_group

Goal: Union DNA snps counts only in evidence and pvalue>0.05 from snp count binom test
	Code: 	/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/Z5_merge_readsInGroup.sh
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/Z6_yps5rmB_gc.sh
	Input:	
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_rm11_B/D_*_[AH]_nosp.pileup.BPNAME
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/hst_yps128_5/D_*_[AH]_nosp.pileup.BPNAME

			/home/xinwen/cloud/yeastAse_redo/Z1_out/rm11_B_snpls_evdp_group
			/home/xinwen/cloud/yeastAse_redo/Z1_out/yps128_5_snpls_evdp_group
	Output:
			/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/3_5Bsnp_yr/D_*_[AH]_nosp.yr.gCount

Goal: Plot after-filter snp DNA sum reads, unioned
	Code:	/home/xinwen/cloud/yeastAse_redo/Scripts/Z1_yr_sum.R
	Input:	/home/xinwen/mnt/nnp/1_1_1_cali_hisat2/3_5Bsnp_yr/D_*_[AH]_nosp.yr.gCount
	Output:	/home/xinwen/cloud/yeastAse_redo/Z1_out/GroupSum_alleleCounts.pdf

Goal: Get rm11_B; yps128_5 annotation, and yps128 snp list with correspond gene name
	Code:	[HPC] && 
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/Z9_label_snp_gene.sh

	Input: 	/home/xinwen/mnt/nnp/0_3_0_pilon_busco/rm11_B.fasta   
			/home/xinwen/mnt/nnp/0_3_0_pilon_busco/yps128_5.fasta

			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/GCF_000146045.2_R64_genomic.gtf
			???/scer.gff
			???/scer.fasta
	
	Output:	/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/rm11_B_scer/rm11_B.gff
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_scer/yps128_5.gff

			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/rm11_B_scer/rm11_B.gtf				
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_scer/yps128_5.gtf

			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/rm11_B_scer/rm11_B.exons
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_scer/yps128_5.exons

			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/rm11_B_scer/rm11_B.splice_sites
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_scer/yps128_5.splice_sites

			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_snp_gene_nofilter.txt

Goal: Get filtered snp list from DNA data, and add gene Name for snp list
	Code: 	/home/xinwen/cloud/yeastAse_redo/Scripts/Z2_pre_snplist_gene.R
	Input: 	/home/xinwen/cloud/yeastAse_redo/Z1_out/cal_evd_p_allinfo.RData
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_snp_gene_nofilter.txt
	Output:	
			/home/xinwen/cloud/yeastAse_redo/Z2_out/rm11_B_snpls_evdp_gene
			/home/xinwen/cloud/yeastAse_redo/Z2_out/yps128_5_snpls_evdp_gene

			/home/xinwen/cloud/yeastAse_redo/Z2_out/rm11_B_snpls_evdp
			/home/xinwen/cloud/yeastAse_redo/Z2_out/yps128_5_snpls_evdp


Goal: Map RNA expression to both genomes
	Code: 	/home/xinwen/mnt/nnp/2_0_2_express_hisat2/Z1_hisat-build.sh
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/Z2_hisat2_express_rm11_B.sh
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/Z2_hisat2_express_yps128_5.sh
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/Z6_yps5rmB_gc.sh

			HPC~/software/count_pileup.py 
			# count by snp
			HPC~/software/group_reads.py	
			# union by group/gene
			HPC~/software/yps5rmB_gc.py
			# match corresponding snps in two alleles
	
	Input:	/home/xinwen/mnt/nnp/0_3_0_pilon_busco/rm11_B.fasta   
			/home/xinwen/mnt/nnp/0_3_0_pilon_busco/yps128_5.fasta
			
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/rm11_B_scer/rm11_B.exons
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_scer/yps128_5.exons

			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/rm11_B_scer/rm11_B.splice_sites
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_scer/yps128_5.splice_sites

			/home/xinwen/mnt/nnp/2_0_0_express/deMx_fixed_zero/*_[AH]_r[12].fq.gz
	
			/home/xinwen/cloud/yeastAse_redo/Z2_out/rm11_B_snpls_evdp_gene
			/home/xinwen/cloud/yeastAse_redo/Z2_out/yps128_5_snpls_evdp_gene

			/home/xinwen/cloud/yeastAse_redo/Z2_out/rm11_B_snpls_evdp
			/home/xinwen/cloud/yeastAse_redo/Z2_out/yps128_5_snpls_evdp


	
	Output: /home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_yps128_5/*_[AH].bam
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_rm11_B/*_[AH].bam
			
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_yps128_5/*_[AH].snpCount
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_rm11_B/*_[AH].snpCount

			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_yps128_5/*_[AH].geneCount
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_rm11_B/*_[AH].geneCount

			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/yps5rmB_geneCount/*_[AH].geneCount
	
Goal: Make a snp count file RNA with Gene Name, not unioned
	Code: 	/home/xinwen/cloud/yeastAse_redo/Scripts/Z2_snpCount_RNA.R

	Input:	/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_yps128_5/*_[AH].snpCount
			/home/xinwen/mnt/nnp/2_0_2_express_hisat2/hst_rm11_B/*_[AH].snpCount

			/home/xinwen/cloud/yeastAse_redo/Z1_out/posmap52.RData
			
			/home/xinwen/mnt/nnp/0_3_1_assembly_crossMap/yps128_5_snp_gene_nofilter.txt
	Output:
			/home/xinwen/cloud/yeastAse_redo/Z2_out/exp_snp_gene.RData

Goal: Make a read count data frame per gene, unioned. 
	Code: 	/home/xinwen/cloud/yeastAse_redo/Scripts/Z2_snpCount_RNA.R
	Input: 	/homt/xinwen/mnt/nnp/2_0_2_express_hisat2/yps5rmB_geneCount/*.geneCount
	Output: /home/xinwen/cloud/yeastAse_redo/Z2_out/exph.RData
			/home/xinwen/cloud/yeastAse_redo/Z2_out/Expressio_reads_sum.pdf	

Goal: Find cis-affected genes in exph, by randomly chose a subset.
	Code:	/home/xinwen/cloud/yeastAse_redo/Programs/cis_prmtr_func.R
			/home/xinwen/cloud/yeastAse_redo/Scripts/Z3_pack_express_function.R
			/home/xinwen/cloud/yeastAse_redo/Scripts/Z3_qsub_exph_redo.Rscript.sh

	Input:
			/home/xinwen/cloud/yeastAse_redo/Z3_out/exph_functions.RData
			/home/xinwen/cloud/yeastAse_redo/Z3_out/exph_random_subset.txt
									(both from Code Z3_pack_express_function.R)
	
	Output: /home/xinwen/mnt/nnp/3_0_1_infer_ase/exph_redo/exph_ci_*.RData

Goal: Discovery rate graph for exph
	Code:	/home/xinwen/cloud/yeastAse_redo/Scripts/Z4_exph_disR.R
	Input:	/home/xinwen/mnt/nnp/3_0_1_infer_ase/exph_redo/exph_ci_*.RData
	Output:	/home/xinwen/cloud/yeastAse_redo/Z4_out/fig_exph.pdf

######################### #####################  #######################################
Goal: Now, the data is spliced as follows: For example, in replication level 3:  For lowest seq depth: Sum every 2 replicates in my hybrid data, then get a pool of 18*17 possible samples, randomly pick 3 replicates from the pool, repeat for 25 times;l  for median seq depth, sum every 6 replicates in my hybrid data, get a pool of Choose(18,6), randomly pick 3 replicates from the pool, repeat 25 times. For high seq depth, sum every 12 replicate, get a pool of choose(18,12), randomly pick 3  replicates from the pool, repeat 25 time. ...  The same for replicates level 6, 12 and 18. The seq depth will not be a fix number. But a range.
	Code: 	/home/xinwen/cloud/yeastAse_redo/Scripts/Z5_rearrange_exph.R
	Input:	/home/xinwen/cloud/yeastAse_redo/Z2_out/exph.RData
	Output:	/home/xinwen/cloud/yeastAse_redo/Z5_out/depth_sample.RData

Goal: run cis-infer machinary on the depth_sample data
	Code:	
