#-----------------------#
#		PREP			#
#-----------------------#

library("adegenet"); library("poppr"); library(hierfstat); library("Demerelate")

#Modification of existing function from adegenet
repool_new<- function(genind_obj,vect_pops){
	genind_obj_sep<-seppop(genind_obj)
	genind_obj_merge<-genind_obj_sep[[vect_pops[1]]]
	for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
	genind_obj_merge
}

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Quercus_collab/")
	#setwd(paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))	#setwd(paste("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
	#setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii")
	
pop_sizes<-read.csv("3_oak_pop_sizes.csv",header=F)
species_names<-c("Qgeorgiana","Qgeorgiana_EST","Qgeorgiana_g","Qoglethorpensis","Qboyntonii")
file_names<-c("Qg_wild_w_NC.gen", "Qg_wild_w_NC_ESTSSRs.gen", "Qg_wild_w_NC_gSSRs.gen", "Qo_wild_w_gos.gen", "Qb_wild_w_ALL.gen")
file_names2<-c("Qg_wild_rel_w_NC.csv", "Qg_wild_rel_w_NC_ESTSSRs.csv", "Qg_wild_rel_w_NC_gSSRs.csv", "Qo_wild_rel_w_gos.csv", "Qb_wild_rel_w_ALL.csv")

min_pop_size_keep<-c(2,5,10)

for (min_p in min_pop_size_keep){

par(mfrow=c(2,3))	#for IBD plots

for (sp in 1:length(species_names)){
	this_species<-species_names[sp]
		
		setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Quercus_collab/")
#------------------------------------------------------------------------------#
#				SECTION ONE
#------------------------------------------------------------------------------#
	
##############
#  CLONES	 #
##############

#import the data- note individuals with no data at all are dropped
ade_test<-read.genepop(paste("genetic_data/",file_names[sp],sep=""),ncode=3)
pop_names<-levels(ade_test@pop)

#First put into poppr format
popr_test <- as.genclone(ade_test)
strata(popr_test) <- other(popr_test)$population_hierarchy[-1]
list_a<-mlg.id(popr_test)
#Function to pull out individual indices where clone length greater than 1
clone_index<-which(sapply(list_a,function(x) length(x)>1))
write.csv(list_a[clone_index],file=paste0(species_names[sp],"clones.csv"))

#This removes clones and then saves as new file for Genealex if desired
popr_nocl<-clonecorrect(popr_test,strata=~Pop)
#genind2genalex(genclone2genind(popr_nocl),file="QH_clone_free.csv")

#Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
GI_nocl<-genclone2genind(popr_nocl); 	GP_nocl<-genind2genpop(GI_nocl)

						  
###################
#  BASIC STATS	  #
###################
	
#narrow down to populations with >10 individuals- call these GP_sub, GI_sub for subset
pop_keep<- which(as.vector(table(GI_nocl@pop)>=min_p))
GI_sub<-repool_new(GI_nocl,pop_keep);	GP_sub<-GP_nocl[pop_keep,]
samp_size<-table(GI_nocl@pop)[pop_keep]

#now we are down to no clones and decently sampled populations, with minimal missing data!
#Summary statistics
par(mfrow=c(2,1))
Qsp_sumstats<-summary(GI_sub)
Qsp_alleles<-Qsp_sumstats$pop.n.all/ length(GI_sub@loc.n.all)
Qsp_all_rich<-colSums(allelic.richness(GI_sub)$Ar)/ length(GI_sub@loc.n.all)	#remember to divide by number of loci!!
Qsp_N_v_MLG<-poppr(popr_test)[pop_keep,2:3]
Qsp_Hexp<-poppr(popr_test)[pop_keep,10]

###################
#  PAIRWISE FST	  #
###################

sm_fst_mat<-as.matrix(pairwise.fst(GI_sub))
rownames(sm_fst_mat)<-pop_names[pop_keep];	colnames(sm_fst_mat)<-pop_names[pop_keep]
sm_fst_mat[sm_fst_mat==0]<-NA
Qsp_pwfst<-apply(sm_fst_mat,2,mean,na.rm=T)
#write.csv(sm_fst_mat,file="fst_mat.csv")


 #####################################
 #	Isolation by distance
 #####################################
 
 Qsp_coords <- read.csv("pop_locations.csv", sep=",", header=T, stringsAsFactors = FALSE)[,c(6,7)] #Loading locations
#Go through each species including for georgiana do all markers, EST, and gSSRs
rows_to_do<-list(1:9,10:18,19:27,28:34,35:41)
r<-rows_to_do[[sp]]
 Qsp_dist <- matrix(nrow = length(r), ncol = length(r))

	for(first in 1:length(r)){
	  for(second in 1:length(r)){
		Qsp_dist[first,second] <-  distm(Qsp_coords[r[first],], Qsp_coords[r[second],], fun = distGeo)/1000
	  }
	}
	Qsp_dist[Qsp_dist==0]<-NA
IBD<-lm(sm_fst_mat[lower.tri(sm_fst_mat)]~log(Qsp_dist[lower.tri(Qsp_dist)]))
plot(sm_fst_mat[lower.tri(sm_fst_mat)]~log(Qsp_dist[lower.tri(Qsp_dist)]),ylab="FST",xlab="distance (km)")
abline(IBD,col="violet")
print(summary(IBD))
text(800,.05,paste0("p=",round(summary(IBD)[[4]][2,4],3),"  R2=",round(summary(IBD)[[9]],3)))

###################
#  RELATEDNESS	  #
###################

Qsp_rel<-read.csv(paste0("genetic_data/",file_names2[sp]))
wang.results<-Demerelate(Qsp_rel,object=T,value="wang")
#morans.results<-Demerelate(Qsp_rel,object=T,value="morans")
ritland.results<-Demerelate(Qsp_rel,object=T,value="ritland")
Mxy.results<-Demerelate(Qsp_rel,object=T,value="Mxy")
#means of pairwise relatedness within each population
wang.means<-unlist(lapply(wang.results$Empirical_Relatedness,mean))[pop_keep]
ritland.means<-unlist(lapply(ritland.results$Empirical_Relatedness,mean))[pop_keep]
#morans.means<-unlist(lapply(morans.results$Empirical_Relatedness,mean))[pop_keep]
Mxy.means<-unlist(lapply(Mxy.results$Empirical_Relatedness,mean))[pop_keep]
						  
#calculate proportion of population that are half sibs or greater, and plot this
half_sib<-function(x) (sum(x>0.25)/length(x))

##############################################
#	Observed heterozygosity from diversity	 #
##############################################

a<-divBasic(paste("genetic_data/",file_names[sp],sep=""), bootstrap=100)
Qsp_Ho<-(a$Ho[nrow(a$Ho),])[pop_keep]
Qsp_Fis<-unlist(lapply(a$fis, function(x) x[length(x[,1]),1])[pop_keep])

			
#####################################						  
#  OUTPUT SUMMARY TABLE, DO T TESTS #
#####################################			
setwd(paste0("pop_keep_",min_p))			  
write.csv(file=paste(species_names[sp],"_summ_stats.csv",sep=""), cbind(pop_names[pop_keep],Qsp_N_v_MLG,Qsp_Hexp,Qsp_Ho,Qsp_Fis,Qsp_alleles,Qsp_all_rich,Qsp_pwfst,wang.means,ritland.means,Mxy.means,pop_sizes[pop_sizes[,1]==species_names[sp],5:6][pop_keep,]))
#cbind(pop_names[pop_keep],Qsp_N_v_MLG,Qsp_Hexp,Qsp_alleles,Qsp_all_rich,Qsp_pwfst,wang.means,ritland.means,Mxy.means,pop_sizes[pop_sizes[,1]==species_names[sp],5:6][pop_keep,])
#,as.numeric(table(Qsp_rel[,2]))
			  
#################
#	CLUSTERING	#
#################

#############
#	DAPC	#
#############
#pdf(file="Qha_dapc1.pdf")
#for (cl in 2:10){
#grp<-find.clusters(GI_sub,max.n.clust=12,n.pca=50,n.clust=cl)
#dapc1 <- dapc(GI_sub, grp$grp,n.pca=50,n.da=50);	scatter(dapc1)
#}
#dev.off()
pdf(file=paste0(file_names[sp],"_str_like.pdf"),width=40,height=9)
for (cl in 2:10){
grp<-find.clusters(GI_sub,max.n.clust=12,n.pca=50,n.clust=cl)
dapc1 <- dapc(GI_sub, grp$grp,n.pca=50,n.da=50);	compoplot(dapc1,show.lab=T)
}
dev.off()
pdf(file=paste0(file_names[sp],"_dapc2.pdf"))
dapc2 <- dapc(GI_sub,n.pca=50,n.da=50);	scatter(dapc2)
dev.off()
pdf(paste0(file_names[sp],"_snapclust.pdf"), height=5, width=18)
for (cl in 2:10){
grp<-find.clusters(GI_sub,max.n.clust=12,n.pca=50,n.clust=cl)
res<-snapclust(GI_sub,cl,pop.ini=grp$grp);		compoplot(res)
}
dev.off()


#############
#	PCA/CA	#
#############
#pdf(file="Qha_pca_bw.pdf")
#pca1<-dudi.pca(df = tab(GI_sub, NA.method = "zero"), scannf = FALSE, nf = 50)
#s.class(pca1$li,pop(GI_sub),xax=1,yax=2,sub="PCA 1-2",csub=2)
#dev.off()
#pdf(file="Qha_pca_col.pdf")
#col <- funky(length(pop(GI_sub)))
#s.class(pca1$li, pop(GI_sub),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)
#dev.off()
pdf(file=paste0(file_names[sp],"_ca.pdf"))
ca1 <- dudi.coa(tab(GP_sub),scannf=FALSE,nf=3)
s.label(ca1$li, sub="CA 1-2",csub=2)
dev.off()

setwd("..")

}

}


dev.off()





#-----------------------------------#
#	HIGHER LEVEL CALCULATIONS		#
#-----------------------------------#


#####################
#Glue all species together into one CSV
#####################

sp<-1
stats_total<-read.csv(file=paste(species_names[sp],"_summ_stats.csv",sep=""))
stats_total<-cbind(species_names[sp],stats_total)
for (sp in 2:length(species_names)){
	stats_total<-rbind(stats_total,
		cbind(species_names[sp],read.csv(file=paste(species_names[sp],"_summ_stats.csv",sep=""))))
	}
write.csv(stats_total,file="stats_total.csv")


setwd("pop_keep_10/")
stats_total<-read.csv(file="stats_total.csv")[,-1]

#Which species has highest genetic diversity 
for (subs in 1:3) {
	if (subs==1) {stats_3_sp<-stats_total[19:41,];	pdf(file="compare_3_sp_nuc.pdf",height=4,width=8)}	#EST_remov... only nuclear genomic markers for georgiana
	if (subs==2) {stats_3_sp<-stats_total[c(1:9,28:41),];	pdf(file="compare_3_sp_all.pdf",height=4,width=8)}	#EST included... all loci for georgiana
	if (subs==3) {stats_3_sp<-stats_total[c(10:18,28:41),];	pdf(file="compare_3_sp_est.pdf",height=4,width=8)}	#EST only for georgiana
	stats_3_sp[,1]<-factor(stats_3_sp[,1])
	anova_res_species<-matrix(ncol=4,nrow=3); rownames(anova_res_species)<-c("exphet", "allrich", "fst")
	colnames(anova_res_species)<-c("p.overall", "p.Qg-Qb", "p.Qo-Qb", "p.Qo-Qg")

	par(mfrow=c(1,3),mar=c(3,4,2,2))	
	boxplot(stats_3_sp[,6]~stats_3_sp[,1], ylab="exp heterozygosity", names=c("Q bo", "Q geo", "Q og"),cex.axis=1.25,cex.lab=1.25)
	anova_res_species[1,1]<-unlist(summary(aov(stats_3_sp[,6]~stats_3_sp[,1]))[[1]][5])[1]
	anova_res_species[1,2:4]<-TukeyHSD(aov(stats_3_sp[,6]~stats_3_sp[,1]))[[1]][,4]
	boxplot(stats_3_sp[,8]~stats_3_sp[,1], ylab="allelic richness", names=c("Q bo", "Q geo", "Q og"),cex.axis=1.25,cex.lab=1.25)
	anova_res_species[2,1]<-unlist(summary(aov(stats_3_sp[,8]~stats_3_sp[,1]))[[1]][5])[1]
	anova_res_species[2,2:4]<-TukeyHSD(aov(stats_3_sp[,8]~stats_3_sp[,1]))[[1]][,4]
	boxplot(stats_3_sp[,9]~stats_3_sp[,1],ylab="pop pw FST", names=c("Q bo", "Q geo", "Q og"),cex.axis=1.25,cex.lab=1.25)
	anova_res_species[3,1]<-unlist(summary(aov(stats_3_sp[,9]~stats_3_sp[,1]))[[1]][5])[1]
	anova_res_species[3,2:4]<-TukeyHSD(aov(stats_3_sp[,9]~stats_3_sp[,1]))[[1]][,4]
	write.csv(anova_res_species, "anova_res_species.csv")
	dev.off()
}
###############################
#Compare the EST and genomic microsatellites for the diversity metrics
###############################

#par(mfrow=c(3,2))
#for (x in c(6,8,9,10,11,12)) {
pdf(file="g_vs_EST_SSR.pdf", height=4.5, width=5)
par(mfrow=c(1,2), mar=c(3,4,1,1)) 
 boxplot(stats_total[11:20,6],stats_total[21:30,6], names=c("EST","genomic"), ylab="heterozygosity" )
 #p=0.019
 boxplot(stats_total[11:20,8],stats_total[21:30,8], names=c("EST","genomic"), ylab="allelic richness")
 #p=0.002
dev.off()
t.test(stats_total[11:20,8],stats_total[21:30,8])
t.test(stats_total[11:20,6],stats_total[21:30,6])
 #print(wilcox.test(stats_total[11:20,x],stats_total[21:30,x]))
 

 
#########################################
# Test genetic statistics ~ population size estimates
###########################################

#Genetic Statistics
#Hexp (5), Ar (6), Fst 87), Wang (9), minPS (12), maxPS (13)
#poulation sizes
#low estimate (12), high estimate (13)


#REGRESSIONS
# for (sp in 1:length(species_names)){
	# these_sum_stats<-read.csv(file=paste(species_names[sp],"_summ_stats.csv"))
	# pdf(file=paste0(species_names[sp],"regr.pdf"))
	# par(mfrow=c(4,2),mar=c(2,2,1,1),oma=c(2,2,1,1))
		# for (x in c(5,7,8,9)) {
			# for (y in c(12,13)) {
				# lm1<-lm(these_sum_stats[,x]~these_sum_stats[,y])
				# plot(these_sum_stats[,x]~these_sum_stats[,y])
				# if (summary(lm1)[[4]][2,4]<0.05) {
					# abline(lm1); mtext(summary(lm1)[[9]],line=3)
					# }
			# }
		# }
	# dev.off()
# }

#BOXPLOTS
setwd("pop_keep_10/")
pdf(file="3_sp_boxplots.pdf",width=8,height=10)
par(mfrow=c(4,4),mar=c(2,4,1,1),oma=c(2,2,3,1))

for (sp in c(1:2,4:5)){
	these_sum_stats<-read.csv(file=paste(species_names[sp],"_summ_stats.csv"))
	colnames(these_sum_stats)<-c("X","pop names","N","MLG","exp heterozygosity","alleles","allelic richness", "pw pop FST", "pw indiv relatedness","ritland","Mxy","psize low","psize high")
		for (x in c(5,7,8,10)) {
		boxplot(these_sum_stats[these_sum_stats[,12]>50,x],(these_sum_stats[these_sum_stats[,12]<=50,x]), cex.main=1,ylab="", names=c("",""))
		pval<-round(unlist(t.test(these_sum_stats[these_sum_stats[,12]>50,x],(these_sum_stats[these_sum_stats[,12]<=50,x]))[3]),3)
		if (pval<0.07) mtext(side=3,paste0("p=",pval))
		if (sp==5) mtext("pop>=50   pop<50",side=1, line=1)
		}
		
	}
mtext("         exp_het                          allelic richness                   pw pop FST                   relatedness   ", side=3, line=0,outer=T)
mtext("          Q boyntonii                            Q oglethorpensis                          Q georgiana (nucl)                  Q georgiana (EST)  ", side=2, line=0,outer=T)
dev.off()
		
		#relatedness xy plot
		#plot(stats_total[,10],stats_total[,11])
		
	
#---------------------------------------#
#	CORRELATION TO CLIMATE				#
#---------------------------------------#
	
	library(stats); library(dismo); library(rgdal)
	
#Download some bioclim data
setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Downloads/wc2.0_2.5m_bio")
 files <- list.files(pattern='tif', full.names=TRUE) #Load climate files
 bioclim2.5 <- stack(files) #Create a raster stack

#load population lat/long
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Quercus_collab/")
 Qsp_loc <- read.csv("pop_locations.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations

clim_reg <- extract(bioclim2.5, Qsp_loc[,6:7]) 

#place to keep regression results (p values) of all regressions run
	reg_pval<-matrix(nrow=dim(clim_reg)[2],ncol=4)
	reg_r2<-matrix(nrow=dim(clim_reg)[2],ncol=4)
	
stats_total<-read.csv(file="stats_total.csv")[,-1]


#Go through each species including for georgiana do all markers, EST, and gSSRs
rows_to_do<-list(1:9,10:18,19:27,28:34,35:41)
for (r in rows_to_do){
stat_list<-c(6,8,9,10)
	for (ss in 1:length(stat_list)){
		for (c in 1:ncol(clim_reg)){
			reg_pval[c,ss]<-summary(lm(stats_total[r,stat_list[ss]]~clim_reg[r,c]))[4][[1]][8]
			reg_r2[c,ss]<-unlist(summary(lm(stats_total[r,stat_list[ss]]~clim_reg[r,c]))[8]	)
		}
}
print(r)
print(which(reg_pval<0.01,arr.ind=T))
#Several are less than 0.05, but none after adjustment for MCT (even a liberal MCT)
print(sum((p.adjust(reg_pval,"BY"))<0.05))
}



#place to keep regression results (p values) of all regressions run
	reg_pval<-matrix(nrow=2,ncol=4)
	reg_r2<-matrix(nrow=2,ncol=4)
	
rows_to_do<-list(1:9,10:18,19:27,28:34,35:41)
for (r in rows_to_do){
stat_list<-c(6,8,9,10)
	for (ss in 1:length(stat_list)){
		for (c in 6:7){
			reg_pval[c-5,ss]<-summary(lm(stats_total[r,stat_list[ss]]~Qsp_loc[r,c]))[4][[1]][8]
			reg_r2[c-5,ss]<-unlist(summary(lm(stats_total[r,stat_list[ss]]~Qsp_loc[r,c]))[8]	)
		}
}
print(r)
print(which(reg_pval<0.01,arr.ind=T))
#Several are less than 0.05, but none after adjustment for MCT (even a liberal MCT)
print(sum((p.adjust(reg_pval,"BY"))<0.05))
}

 
 
 
 
 
 
 

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Quercus_collab/genetic_data/")
source("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/Fa_sample_funcs.R")
file_names<-c("Qg_total.gen", "Qg_total_ESTSSRs.gen", "Qg_total_gSSRs.gen", "Qb_total.gen", "Qo_total.gen")

#This will run over a loop of "Include all alleles (n_to_drop=0)" and "Include only alleles present in more than two copies (n_to_drop=2)"
for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-""
	if (n_to_drop==0) n_drop_file<-"_dr_0"
	
	region_makeup_list<-list(list(1:2,3,4,5:10),list(1:2,3,4,5:10),list(1:2,3,4,5:10),list(1:2,3,4,5:9),list(1:2,3,4,5:8))
	set_garden_p<-c(rep(11,3),10,9) 
	wild_results<-matrix(nrow=length(file_names),ncol=9+1)
	alleles_existing_by_sp<-matrix(nrow=length(file_names),ncol=9)

	for (sp in 1:length(file_names)){
		this_species<-file_names[sp]
		Spp_tot_genind<-read.genepop(this_species,ncode=3)

		#This code compares the wild to various ex situ populations or all the ex situ merged
		#Just put in garden and wild population numbers... currently all gardens are merged into one "population"
		wild_p<-unlist(region_makeup_list[[sp]]); garden_p<-set_garden_p[sp]
		n_ind_W<-table(Spp_tot_genind@pop)[wild_p];  n_ind_G<-table(Spp_tot_genind@pop)[garden_p]; 
		Spp_tot_genpop<-genind2genpop(Spp_tot_genind)
		Spp_tot_genind_sep<-seppop(Spp_tot_genind)
		alleles_cap<-colSums(Spp_tot_genind_sep[[garden_p]]@tab,na.rm=T)

			#Allele categories based only on wild populations (can look at all wild pop'ns or only one if you want)
		allele_cat_tot<-get.allele.cat(Spp_tot_genpop[wild_p], region_makeup_list[[sp]], 2, n_ind_W, glob_only=T,n_drop=n_to_drop)
			#This goes through each allele category and divides the number captured ex situ (alleles_cap) by the number of alleles existing (allele_cat_tot)
			list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")
			
		for (i in 1:9) alleles_existing_by_sp[sp,i]<- (sum((allele_cat_tot[[i]])>0,na.rm=T))
		
		for (l in 1:length(allele_cat_tot)) wild_results[sp,l]<-round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
		
		wild_results[sp,10]<-n_ind_G
	}
	
	wild_results<-cbind(file_names,wild_results)
	write.csv(wild_results,file=paste("ex_vs_in_situ",n_drop_file,".csv",sep=""))
	
} #close num to drop loop		
		


