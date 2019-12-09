# Erin Gilbert
# Created: ~April 2017
# Updated: Oct 17 2019

#Usage: rscript  ALL_SNP_FILE.csv NIL_INFO_FILE.csv TRAIT


############################################
######    Define Environment  #############
############################################


############################################
######    Define Functions     #############
############################################
Corner_text <- function(text, location="topright"){  #usage Corner_text(text="Use Corner_text() function")
  legend(location,legend=text, bty ="n", pch=NA) #allows arguments to be passed 
}


############################################
######    Load All the Data    #############
############################################
#Table ofr SNP values from Song et al.
SNPs<- read.table(file="./Worktable/masters_publication/data/loess/NIL_SoySNP50k_cM_pos.txt", header=T, stringsAsFactors=F, sep="\t")#all NIL SNPs in one file downloaded from Soybase
SNPs$CHROM<-gsub("Gm0", "", SNPs$CHROM)
SNPs$CHROM<-gsub("Gm", "", SNPs$CHROM)


#Table of NILs, genes, pedigree, donor, recurrent, maturity, trait type according to Bernard, and subcollection
NIL_master_table<-read.table(file="./Worktable/masters_publication/mapping/NIL_info_101619_reduced.txt", header=T, stringsAsFactors=F)

NIL_master_table<-data.frame(NIL_master_table) #create usable dataframe
print("NILs loaded...")


############################################
######    Main Mapping Code    #############
############################################


#########################################
####   Select NILs for each trait    ####
#########################################
#Create table of NILs with trait z from the trait_status table. Includes Donor, Recurrent, and wt/mut status

genes<-unique(NIL_master_table$genes)
for(g in 1:length(genes))
  #for(g in 43:length(genes))
{
  trait<-genes[g]
  print(trait)
  #Trait_NILs<-cbind(NIL_master_table[grep(paste(">", trait, "<", sep=""), NIL_master_table$genes), c(1,4,5)],Status=rep(trait_status[z, 2],length.out = length(grep(paste(">", trait, "<", sep=""),NIL_master_table$genes))))
  Trait_NILs<-NIL_master_table[which(NIL_master_table$genes==trait),]
  rownames(Trait_NILs)<-NULL #remove rownames from Trait_NILs
  
  
  if(nrow(Trait_NILs)<=2) #Need at least 3 NILs to map a trait
  {
    write(file="NoNILs.txt", paste(trait, nrow(Trait_NILs), "available NIL(s)", sep=" "), append=T)
    next
  }
  
  
  
  #########################################
  #### Gathering Alleles for Each NIL  ####
  #########################################
  All_alleles<-NULL
  All_alleles<-data.frame(All_alleles)
  print(names(SNPs))
  All_alleles<-subset(SNPs, select=c(CHROM, POS, WP_linkage_pos_plus_5))
  targets<-subset(SNPs, select= as.character(Trait_NILs$NIL))
  targetnils<-Trait_NILs$NIL[which(Trait_NILs$NIL %in% names(SNPs))]
  All_alleles<-cbind(All_alleles,targets) #adds in SNPs for the lines from Trait_NILs
  print("NIL alleles loaded")
  
  #add in the recurrent parent SNPs. Recurrent here is going to also represent non-recurrent lines that are the background for receiving a target trait.
  recurrent<-unique(Trait_NILs$Recurrent)
  recurrent <- recurrent[which(recurrent %in% names(SNPs))]
  All_alleles<-cbind(All_alleles,subset(SNPs, select=as.character(recurrent))) #adding donor parents to table All_alleles
  print("RP alleles loaded")
  
  #add in donor parent SNPs
  donor<-unique(Trait_NILs$Donor)
  donor <- donor[which(donor %in% names(SNPs))]
  All_alleles<-cbind(All_alleles,subset(SNPs, select=as.character(donor)))
  All_alleles<- subset(All_alleles, CHROM %in% c(1:20))
  print("DP alleles loaded")
  
  #########################################
  ##### Scoring Alleles for Each NIL  #####
  #########################################
  Score_holder<-All_alleles #The data frame in which alleles are given scores. Also where all score manipulation happens. Use this as work datarame, then store desired results in appropriate place.
  Score_holder[]<-lapply(Score_holder, as.character)
  print("Score holder made...")
  
  for(i in 1:nrow(Trait_NILs)) #work with each NIL to score it
  {
    NILTarget<-as.character(Trait_NILs$NIL[i]) #NIL being scored
    ParentRec<-as.character(Trait_NILs$Recurrent[i]) #Recurrent parent of NILTarget
    ParentDonor<-as.character(Trait_NILs$Donor[i]) #Donor parent of NILTarget
    print("assign target and parents...")
    print(ParentDonor)
    print(ParentRec)
    print(NILTarget)
    
    if( (ParentDonor %in% names(Score_holder)) & (ParentRec %in% names(Score_holder)) & (NILTarget %in% names(Score_holder)))
    {
      #Each NIL is scored all at once.
      Score_holder[Score_holder[,ParentRec]==Score_holder[,ParentDonor],NILTarget]<- "-" #monomorphic loci receive -
      print("#monomorphic loci receive -")
      Score_holder[Score_holder[,ParentRec]=="H" | Score_holder[,ParentDonor]=="H" | Score_holder[,ParentDonor]=="U" | Score_holder[,ParentRec]=="U",NILTarget]<- "-" #Heterozygous or missing reads in parents give the NIL locus a "-"
      print("Heterozygous or missing reads in parents give the NIL locus a -")
      Score_holder[Score_holder[,NILTarget]==Score_holder[,ParentRec] & Score_holder[,NILTarget]%in% c("A","T","C","G"),NILTarget]<-0 #if NIL matches recurrent parent at a polymorphic locus
      Score_holder[Score_holder[,NILTarget]==Score_holder[,ParentRec] & Score_holder[,NILTarget]%in% c("A","T","C","G") & Trait_NILs$Status[i]==0,NILTarget]<-0 #if NIL matches recurrent parent at a polymorphic locus & carries wildtype
      Score_holder[Score_holder[,NILTarget]==Score_holder[,ParentDonor] & Score_holder[,NILTarget]%in% c("A","T","C","G") & Trait_NILs$Status[i]==1,NILTarget]<-1 #if NIL matches donor parent at a polymorphic locus & carries mutant
      Score_holder[Score_holder[,NILTarget]==Score_holder[,ParentDonor] & Score_holder[,NILTarget]%in% c("A","T","C","G"), NILTarget]<-1 #Remaining polymorphic SNPs must match donor parent NIL matches donor parent
      Score_holder[Score_holder[,NILTarget]=="H",NILTarget]<-0.5 #Heterozygous reads in the NIL receive 0.5
      Score_holder[Score_holder[,NILTarget]=="U",NILTarget]<-"-" #missing reads in the NIL receive "-"
      Score_holder[Score_holder[,NILTarget]%in% c("A","T","C","G"), NILTarget]<-"-" #remaining receive "-"
    }
    
  }
  
  print("Alleles Scored...")
  
  Score_holder<-Score_holder[,1:(3+length(targetnils))] #remove the parents
  Raw_score_backup<-Score_holder
  rows = apply(Raw_score_backup[, 4:(3+length(targetnils))], 1, function(i) !all(i=="-")) #remove rows with all "-"
  Score_holder<-Raw_score_backup[rows,]
  rows = apply(Score_holder[, 4:(3+length(targetnils))], 1, function(i) !all(i=="*")) #remove rows with all "*"
  Score_holder<-Score_holder[rows,]
  UnimputedScores <- Score_holder
  
  #########################################
  ##### Impute Scores for Each NIL  #####
  #########################################
  
  #The first and last scores are set to zero to allow linear interpolation at the ends of the chromosomes
  Score_holder[1, 4:ncol(Score_holder)]<-"0" 
  Score_holder[nrow(Score_holder), 4:ncol(Score_holder)]<-"0"
  
  for( z in 4:ncol(Score_holder)){
    
    x<-which(Score_holder[,z]!="-" & Score_holder[,z]!="*")
    y<-which(Score_holder[,z] =="-")
    
    
    while(length(y)>0)
    {
      down<-x[min(which(x>y[1]))] #Finds the closest score downstream
      up<-x[max(which(x<y[1]))] #Finds the closest score upstream
      
      #Calculate the rate of change for the scores to be imputed
      slope=(as.numeric(Score_holder[down,z])- as.numeric(Score_holder[up,z]))/(down-up)
      
      if(as.numeric(Score_holder[down,z])!=as.numeric(Score_holder[up,z]))
      {
        impute<-seq(as.numeric(Score_holder[up,z]), as.numeric(Score_holder[down,z]), by=slope)
        impute<-impute[-1]
        impute<-as.character(impute)
      }else{
        impute<-rep(Score_holder[up,z], down-up)
      }
      
      Score_holder[(up+1):down,z]<- impute
      
      
      y<-which(Score_holder[,z]=="-")
    }
    
  }
  print("Scores imputed")
  
  
  #########################################
  ######## Calculating Significance #######
  ######## of Each Locus Score Set  #######
  #########################################
  
  ##The probability (P) of each locus allele originating from the donor parent can be represented by a binomial equation which takees into account the number of reccurent scores (0) versus donor scires (1) for that locus in the NILs selected for the trait of interest. The equation characterizes the deviation of the score results from the expected 1-0.5^6 (or about 98.6%) proportion of recurrent alleles according to the backcross scheme of the pedigrees. -log10(P) is calculated and stored here. ##
  
  Score_holder["NegLogP"] <-NA
  NIL_count<-length(targetnils)
  for(c in 1:nrow(Score_holder)) #Binomial calculations in this loop
  {
    Recurrent_scores_per_SNP<-sum(Score_holder[c,4:ncol(Score_holder)-1]==0)
    Donor_scores_per_SNP<- sum(Score_holder[c,4:ncol(Score_holder)-1]==1)
    Total_included_scores <- Recurrent_scores_per_SNP + Donor_scores_per_SNP
    Score_holder$NegLogP[c]<- -log10(choose(Total_included_scores, Recurrent_scores_per_SNP )*(0.984375^Recurrent_scores_per_SNP)* (0.015625^Donor_scores_per_SNP))
  }
  print("Significance calculated...")
  
  
  #####################################
  ######## Choosing interval  #########
  ######## based on -Log10(P)  ########
  #####################################
  
  All_scores<-Score_holder
  All_scores.keyrows<-rownames(All_scores)
  rownames(All_scores) <- seq(length=nrow(All_scores))
  intervals.raw<- data.frame(matrix(ncol=6, nrow=0), stringsAsFactors = F)
  colnames(intervals.raw)<- c("Chromosome", "Start", "End", "NegLog10Pmax", "StartcM", "EndcM")
  intervals<- data.frame(matrix(ncol=6, nrow=0), stringsAsFactors = F)
  colnames(intervals)<- c("Chromosome", "Start", "End", "NegLog10Pmax", "StartcM", "EndcM")
  all_max_chr<-as.numeric(unique(All_scores[which(All_scores$NegLogP>4),1]))
  
  for(max_chr in all_max_chr)
  {
    Score_holder<-subset(All_scores, CHROM==max_chr)
    Score_holder["ROWS"]<-row.names(Score_holder)
    rownames(Score_holder) <- seq(length=nrow(Score_holder))
    max_index<-which(Score_holder$NegLogP>4)
    #write.table(Score_holder, file="T_test.txt", quote=F, row.names = F, sep="\t")
    
    for(j in max_index)
    {
      start<-end<-j
      start_found<-F
      while(!start_found & start>1)
      {
        start<-start-1
        if(sum(Score_holder[start,4:(ncol(Score_holder))-1]==0)>0) ##if there are more than one zero
        {
          start_found<-T
          break
        }
      }
      
      
      end_found<-F
      while(!end_found & end<nrow(Score_holder))
      {
        end<-end+1
        if(sum(Score_holder[end,4:(ncol(Score_holder)-1)]==0)>0)
        {
          end_found<-T
          break
        }
      }
      
      
      if(start<=20)
      {
        start<-21
      }
      if(end==nrow(Score_holder))
      {
        end<-(nrow(Score_holder)-20)
      }
      intervals.raw<-rbind(intervals.raw, data.frame(Locus=trait,Chromosome=max_chr, Start=Score_holder$POS[start], End=Score_holder$POS[end], NegLog10Pmax=max(Score_holder$NegLogP[start:end]), StartcM=Score_holder$WP_linkage_pos_plus_5[start], EndcM=Score_holder$WP_linkage_pos_plus_5[end]))
      if(end-start>=2)
      {intervals<-rbind(intervals, data.frame(Locus=trait,Chromosome=max_chr, Start=Score_holder$POS[start], End=Score_holder$POS[end], NegLog10Pmax= max(Score_holder$NegLogP[start:end]), StartcM=Score_holder$WP_linkage_pos_plus_5[start], EndcM=Score_holder$WP_linkage_pos_plus_5[end]))}
    }
  }
  
  intervals<-intervals[!duplicated(intervals),]
  
  All_scores_backup<-All_scores
  All_scores<-All_scores[order(as.numeric(All_scores$CHROM), All_scores$POS),]
  
  if(nrow(intervals)==0)
  {
    write(paste(trait, " could not be mapped due to small intervals"),file="NoNILs.txt", append=T)
    
    #Plotting stuff
    chromosome_start_indices<-integer(0)
    for(i in c(1:20))
    {
      chromosome_start_indices<-c(chromosome_start_indices,min(which(All_scores$CHROM==i)))
    }
    chromosome_start_indices<-sort(as.numeric(chromosome_start_indices))
    
    chromosome_start_indices<- c(chromosome_start_indices, nrow(All_scores))
    chromosome_midpoint_indices<-integer(0)
    for(i in 2:21)
    {
      chromosome_midpoint_indices<-c(chromosome_midpoint_indices, (chromosome_start_indices[i] + chromosome_start_indices[i-1])/2)
    }
    
    jpeg(paste(trait,"-logP_unmapped.jpeg", sep="_"),  width= 1000, height=480, units="px")
    plot(All_scores$NegLogP, xaxt='n',  xlab= "Chromosome", ylab="-Log(P)", main= paste(trait, " Imputed NegLogP Scores"), col=ifelse((as.numeric(All_scores$CHROM) %%2)==0, "black","grey"))
    
    points(All_scores$NegLogP, col=ifelse((as.numeric(All_scores$CHROM) %%2)==0, "black","steelblue"))
    axis(side=1, labels=c(1:20), at = chromosome_midpoint_indices)
    Corner_text(text=paste( "Number of NILs:", nrow(Trait_NILs)))
    dev.off()
    
    
  }else{ 
    print("Intervals selected")
    
    #####################################
    ##### Final Interval Clean-up  ######
    #####################################
    interval_chr<-unique(intervals$Chromosome)
    intervals.final<- data.frame(matrix(ncol=5, nrow=0))
    colnames(intervals.final)<- c("Trait","Chromosome", "Start", "End", "NegLog10Pmax")
    if(nrow(intervals)==1)
    {
      intervals.final<-data.frame(Trait= trait, Chromosome=intervals$Chromosome[1], Start=intervals$Start[1], End=intervals$End[1], NegLog10Pmax=intervals$NegLog10Pmax[1],)
    }else{
      i<-1
      while(i<=nrow(intervals)) #stepping through intervals and looking for continuous intervals
      {
        final<-intervals[i,]
        done<-F
        while(!done)
        {
          if (i<nrow(intervals))
          {
            #if (as.numeric(intervals$Chromosome[i])==as.numeric(intervals$Chromosome[i+1]) & (as.numeric(intervals$Start[i+1])-as.numeric(intervals$End[i])<=2000000))
            if (as.numeric(intervals$Chromosome[i])==as.numeric(intervals$Chromosome[i+1]) & (as.numeric(as.character(intervals$StartcM[i+1]))-as.numeric(as.character(intervals$EndcM[i])))<=1.6)
            {final[4]<-intervals$End[i+1]
            if(final[5]<intervals$NegLog10Pmax[i+1])
            {final[5]<-intervals$NegLog10Pmax[i+1]}
            }
            else {done<-T}
          }else {done<-T}
          i<-i+1
        }
        
        intervals.final<-rbind(intervals.final, data.frame(Trait= final[1], Chromosome=final[2], Start=final[3], End=final[4], NegLog10Pmax= final[5]))
      }
    }
    write.table(intervals.final,"intervals.txt", sep="\t", row.names=F, col.names=F, append=T, quote=F)
    
    
    
    #####################################
    ###### Create Mapping Visuals  ######
    #####################################
    chromosome_start_indices<-integer(0)
    for(i in c(1:20))
    {
      chromosome_start_indices<-c(chromosome_start_indices,min(which(All_scores$CHROM==i)))
    }
    chromosome_start_indices<-sort(as.numeric(chromosome_start_indices))
    
    chromosome_start_indices<- c(chromosome_start_indices, nrow(All_scores))
    chromosome_midpoint_indices<-integer(0)
    for(i in 2:21)
    {
      chromosome_midpoint_indices<-c(chromosome_midpoint_indices, (chromosome_start_indices[i] + chromosome_start_indices[i-1]) /2)
    }
    
    
    jpeg(paste(trait,"logP.jpeg", sep="_"),  width= 1000, height=480, units="px")
    plot(All_scores$NegLogP, xaxt='n',  xlab= "Chromosome", ylab="-Log(P)", main= paste(trait, " -Log(P) Scores"), col=ifelse((as.numeric(All_scores$CHROM) %%2)==0, "green","brown"))
    points(All_scores$NegLogP, col=ifelse((as.numeric(All_scores$CHROM) %%2)==0, "black","steelblue"))
    axis(side=1, labels=c(1:20), at = chromosome_midpoint_indices)
   
    abline(h=4, col="red")
    Corner_text(text=paste( "Number of NILs:", length(targetnils)))
    dev.off()
    
    
  }
  
  All_scores<-All_scores_backup
  
  
  #####################################
  ##### Create Excel Output File  #####
  #####################################
  
  output.scores<-NULL
  output.alleles<-NULL
  

  output.scores<-cbind(All_scores,UnimputedScores)
  output.alleles<-All_alleles[All_scores.keyrows,]
  spacer<-c("-","-")
  spacer_a<-append(spacer, "AlleleType")
  spacer_r<-append(spacer, "Recurrent")
  spacer_d<-append(spacer, "Donor")
  allelic_labels<-NULL
  recurrent_labels<-NULL
  donor_labels<-NULL

  recurrent_labels<-append(spacer_r,Trait_NILs$Recurrent)
  recurrent_labels<-append(recurrent_labels, rep("+", length(unique(Trait_NILs$Recurrent))))
  recurrent_labels<-append(recurrent_labels, rep("-", length(unique(Trait_NILs$Donor))))
  recurrent_labels<-append(recurrent_labels,Trait_NILs$Recurrent)
  donor_labels<-append(spacer_d,Trait_NILs$Donor)
  donor_labels<-append(donor_labels, rep("-", length(unique(Trait_NILs$Recurrent))))
  donor_labels<-append(donor_labels, rep("+", length(unique(Trait_NILs$Donor))))
  donor_labels<-append(donor_labels,Trait_NILs$Donor)
  header<-do.call("rbind", list(recurrent_labels,donor_labels))
  header.final<- header
  output.scores<-subset(output.scores, !is.na(WP_linkage_pos_plus_5 )) 
  output.scores[output.scores==0.0]<-0
  

  output.final<-cbind(output.alleles, output.scores[,4:(ncol(output.scores))])
  output_file<-paste(trait, "_results.txt", sep="")
  write.table(header.final, file=output_file, append=T, row.names=F, col.names = F, quote=F, sep="\t")
  write.table(output.final, file=output_file, append=T, row.names=F, quote=F, sep="\t")
  
  
}
