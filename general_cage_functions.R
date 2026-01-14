###############
## Libraries
###############
suppressWarnings(suppressMessages(library(multcomp)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(doMC)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(itertools)))
suppressWarnings(suppressMessages(library(RColorBrewer)))

#suppressWarnings(library(lme4))

###############
## Utilities
##############

## redefine list function so you can return multiple variables from a function with one assignment call
#########################
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

## chop a string and return a specified field
#########################
chop=function(x,splitchar,field) {sapply(strsplit(x,splitchar),"[",field)}

## get sizes of all objects in memory, sorted
get_obj_sizes=function(){print("sort( sapply(ls(),function(x){format(object.size(get(x)), units = 'Mb')}))")}

## parse file paths and sync with remote teamdrive
#########################
parse_filepaths=function(filepaths){
  filedir=paste0(localDir,"/",dirname(filepaths[1])); 
  drivedir=paste(ifelse(length(filepaths)==3,filepaths[3],teamdriveDir),
                 ifelse(length(filepaths)==1,dirname(filepaths[1]),dirname(filepaths[2])),
                 sep="/"
  ); 
  include=ifelse(length(filepaths)==1,basename(filepaths[1]),basename(filepaths[2]))
  return(list(filedir,drivedir,include))
}

download_files=function(filepaths,dryrun=FALSE){
  list[filedir,drivedir,include] <- parse_filepaths(filepaths)
  system(paste("rclone copy",drivedir,filedir,paste0("--include=",include),ifelse(dryrun,"--dry-run","")))
  return(system(paste0("ls -lh ",localDir,"/",filepaths[1]),intern = TRUE))
}

upload_files=function(filepaths,dryrun=FALSE){
  list[filedir,drivedir,include] <- parse_filepaths(filepaths)
  system(paste("rclone copy",filedir,drivedir,paste0("--include=",include),ifelse(dryrun,"--dry-run","")))
  return(system(paste0("rclone ls ",drivedir," --include=",include),intern = TRUE))
}

check_files=function(filepaths,report=FALSE){
  list[filedir,drivedir,include] <- parse_filepaths(filepaths)
  if(!dir.exists(filedir)){dir.create(filedir,recursive = TRUE)}
  if(report){return(suppressWarnings(system(paste("rclone check",drivedir,filedir,paste0("--include=",include),"2>&1"),intern=T)))}
  return(length(grep("Failed",suppressWarnings(system(paste("rclone check",drivedir,filedir,paste0("--include=",include),"2>&1"),intern=T))))>0)
}

##############
## Data wrangling
##########

## load data
load_sampKey=function(sampleKeyFile){
  if(check_files(sampleKeyFile[1])){download_files(sampleKeyFile)}
  sampKey=fread(sampleKeyFile[1]);
  bacteria=sampKey$treatment=="Aceto" | sampKey$treatment=="Lacto"
  sampKey=sampKey[!bacteria,]
  sOrder=order(sampKey$tpt,sampKey$treatment,sampKey$cage)
  return(sampKey[sOrder,])
}
load_orchard_af_data=function(year,type="HAFs",set="cages",treatment=NULL){
  load(paste0("/mnt/cages/orchard_",year,"/Rdata/",type,".orch",substr(year,3,4),"_",set,".Rdata"))
  samps$year=year; samps$set=set; samps$type=type
  if(!is.null(treatment)){mask=samps$treatment==treatment; afmat=afmat[,mask];samps=samps[mask,]}
  return(list(sites,samps,afmat))
}
load_LIsea_af_data=function(){
  load(file="/mnt/cages/ref_data/heather_sea_cli/mel_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf.LIsamps_FREQ.Rdata")
  return(list(sites,samps,afmat))
}

load_cli_af_data=function(){
  load(file="/mnt/cages/ref_data/heather_sea_cli/mel_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf.clinal_samps_FREQ.Rdata")
  return(list(sites,samps,afmat))
}

load_glm_data=function(dataset,shuff=FALSE,renameCols=TRUE){
  filename=paste0("/mnt/cages/orchard_20",dataset,"/Rdata/glm.Ecages",ifelse(shuff,"_shuffle",""),".Rdata")
  switch(dataset,
         "sea"={
           load("/mnt/cages/ref_data/heather_sea_cli/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm.Rdata")
           df.glm=df.sea
         },
         "cli"={
           load("/mnt/cages/ref_data/heather_sea_cli/mel_clinal_spring.glm.noheader.Rdata")
           df.glm=df.cli
         },
         {
           load(filename)
           if(renameCols){
            coefCols=colnames(df.glm)[grepl("coef\\.",colnames(df.glm))]
            pCols=colnames(df.glm)[grepl("p\\.",colnames(df.glm))]
            df.glm<- df.glm %>% 
             rename_at(vars(coefCols), ~ paste0("coef.orch",dataset,".",chop(coefCols,"\\.",2))) %>% 
             rename_at(vars(pCols), ~ paste0("p.orch",dataset,".",chop(pCols,"\\.",2))) 
           }
        }
  )
  for(ii in grep("p\\.",colnames(df.glm))){df.glm[[ii]][df.glm[[ii]]==0]=NA;ii}
  return(df.glm)
}

load_meanaf_data=function(dataset,shuff=FALSE,renameCols=TRUE){
  filename=paste0("/mnt/cages/orchard_20",dataset,"/Rdata/meanaf.",ifelse(shuff,"Ecages_shuffle",paste0("orch",dataset)),".Rdata")
  load(filename)
  if(renameCols){colnames(meanaf)=paste0("orch",dataset,".",colnames(meanaf))}
  return(meanaf)
}

load_genes=function(){
  load("/mnt/cages/ref_data/d_mel/genes.Rdata")
  return(genes)
}

## make af data object
make_afData=function(thresh.altrd,thresh.minrd){
  if(check_files(afFiles)){download_files(afFiles)}
  if(check_files(repeatMaskerFiles)){download_files(repeatMaskerFiles)}
  
  samps=load_sampKey(sampleKeyFile) 
  afsamps=chop(basename(fread(system(paste0("ls ",afFiles[1]," | grep .af.meanrd$"),intern=T)[1])$V1),"_",1)
  sMask=match(samps$sampID,afsamps)
  
  cat("**Filtering sites from orchard_2017 dataset \n",file=siteFilterLog,append=FALSE)  
  #read af and filter
  df.af=foreach(chrom=chroms,.combine=rbind)%do%{
    rm(contammask,roundmask,foundermask, repmask,rdmask,refmask)
    
    af=fread(system(paste0("ls ",afFiles[1]," | grep ",chrom,".af$"),intern=T))
    sites=af[,1:5,with=F]
    rd=as.matrix(af[,seq(7,ncol(af),2),with=F])
    af=as.matrix(af[,seq(6,ncol(af),2),with=F])
    colnames(rd)=afsamps
    colnames(af)=afsamps
    cat(nrow(sites),",chrom",chrom,",chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## filter sites with potential human contamination from barcode swapping
    contam=fread(humanContamFilter[1])
    contammask=!(sites$pos%in%contam$pos[contam$chrom==chrom])
    sites=sites[contammask,];af=af[contammask,];rd=rd[contammask,]
    cat(sum(!contammask),",humanreads,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## filter out sites with non-overlapping allele freqs between round1 samps and round2 samps
    R1samps=afsamps%in%samps$sampID[samps$sequencingRd=="R1"]
    R1max=apply(af[,R1samps],1,max)
    R1min=apply(af[,R1samps],1,min)
    R2max=apply(af[,!R1samps],1,max)
    R2min=apply(af[,!R1samps],1,min)
    batchmask=!((R1min>R2max) | (R2min>R1max))
    sites=sites[batchmask,];af=af[batchmask,];rd=rd[batchmask,];
    cat(sum(!batchmask),",batcheffects,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## subset to relevant samps
    af=af[,sMask];rd=rd[,sMask];
    
    ## filter out sites not present in the baseline samples
    foundermask=apply(af[,samps$treatment=="Founder"],1,min)>thresh.founder
    sites=sites[foundermask,];af=af[foundermask,];rd=rd[foundermask,]
    cat(sum(!foundermask),",notstanding variation,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## filter repeat regions
    df.repeat<-fread(system(paste0("ls ",repeatMaskerFiles[1]," | grep ",chrom,".txt$"),intern=T))
    pos=sites$pos
    cc=cut(df.repeat$genoStart,c(0,pos),labels=FALSE)
    dd=cut(df.repeat$genoEnd,c(0,pos),labels=FALSE)
    cc[is.na(cc)]=max(cc,na.rm=TRUE)+1
    dd[is.na(dd)]=max(dd,na.rm=TRUE)+1
    pos.filter=eval(parse(text=paste0("c(",paste0(cc[dd>cc],":",dd[dd>cc]-1,collapse=","),")")))
    repmask=!(1:length(pos))%in%pos.filter
    sites=sites[repmask,];af=af[repmask,];rd=rd[repmask,]
    cat(sum(!repmask),",repeatregions,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    # filter out low + high rd sites
    ### since there are a few samps with low overall rd, only examine samps with mean overall rd>thresh.sampmeanrd
    sampsToExamine=colMeans(rd)>thresh.sampmeanrd
    rdmask=rowSums(rd[,sampsToExamine]<thresh.minrd)==0 & rowSums(rd>thresh.maxrd)==0 
    sites=sites[rdmask,];af=af[rdmask,];rd=rd[rdmask,]
    cat(sum(!rdmask),",lowreaddepth,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## keep only sites where ref allele is present
    refmask=(sites$ref==sites$major) | (sites$ref==sites$minor)
    sites=sites[refmask,];af=af[refmask,];rd=rd[refmask,]
    cat(sum(!refmask),",refallelenotfound,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    cbind(sites,af,rd)
  }
  
  #save as Rdata
  sites=df.af[,1:5,with=F]
  afmat=as.matrix(df.af[,5+1:nrow(samps),with=F])
  rd=as.matrix(df.af[,5+nrow(samps)+1:nrow(samps),with=F])
  flipMe=(sites$ref != sites$major)
  afmat[flipMe,]=1-afmat[flipMe,]
  sites$minor[flipMe]=sites$major[flipMe]
  colnames(sites)[5]="alt";sites=sites[,-4,with=F]
  save(sites,samps,afmat,rd,file=afData)
  upload_files(afData)
}

### make pooldata object for use with poolfstat
make_pooldata=function(afmat,rd,sites,sampIDs,poolSize){
  make.example.files(writing.dir=tempdir())
  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15),)
  pooldata@refallele.readcount=round(rd*(1-afmat))
  pooldata@readcoverage=rd
  pooldata@snp.info=as.matrix(sites)
  pooldata@poolsizes=rep(poolSize,ncol(afmat))
  pooldata@nsnp=nrow(afmat)
  pooldata@npools=ncol(afmat)
  pooldata@poolnames=sampIDs
  return(pooldata)
}

## merge two dataframes by chromosome and position
merge.by.pos.all <- function(a, b) {  
  if("chrom"%in%colnames(a) && "chrom"%in%colnames(b)){
    merge(a, b, by=c('chrom','pos'), all=TRUE,suffixes=c('', ncol(a)))
  } else { merge(a, b, by='pos', all=TRUE,suffixes=c('', ncol(a)))}
}
merge.by.pos <- function(a, b) {  
  if("chrom"%in%colnames(a) && "chrom"%in%colnames(b)){
    merge(a, b, by=c('chrom','pos'),suffixes=c('', ncol(a)))
  } else {merge(a, b, by='pos',suffixes=c('', ncol(a)))}
}

## apply a function to all pairs of columns in a matrix
fun.mat=function (x,myfun) {
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  
  ncy <- ncx <- ncol(x)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(i)) {
      x2 <- x[, i]
      y2 <- x[, j]
      r[i, j] <- myfun( x2, y2)
    }
  }
  r <- r + t(r) - diag(diag(r))
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  r
}


###############
## Transformations/Calculations
###############

## calculate N-effective (Feder 2012), the variance of frequency estimates,
## approximated by calculating the effective number of observations at a given locus, conditional on read depth and number of chromosomes in the sample
n_eff=function(rd,pooledChroms){
  round((rd*pooledChroms-1)/(rd+pooledChroms))
}

### apply FDR adjustment to treatment divergence p-values in a data frame
fdr_adjust=function(pvals){
  rank(pvals)*pvals
}

fdr_adjust_all=function(df,comparisons){
  pcols=match(paste0("p.",comparisons),colnames(df))
  fdr=foreach(pp=pcols,.combine=cbind)%do%{
    fdr=fdr_adjust(df[[pp]])
    fdr[fdr>1]=1
    fdr
  }
  colnames(fdr)=colnames(df)[pcols]
  fdr
}

### get effect size for treatment/tpt pairs
effect_size=function(afmeans,comparisons,treatment=NULL,tpt=NULL){
  groups=strsplit(comparisons,'_')
  if(!is.null(treatment)){groups=lapply(groups,function(x){paste0(treatment,".t",x)})}
  if(!is.null(tpt)){groups=lapply(groups,function(x){paste0(x,".t",tpt)})}
  #print(groups)
  es=do.call(cbind,lapply(groups,function(gg){abs(afmeans[,gg[1]]-afmeans[,gg[2]])}))
  colnames(es)=comparisons
  es
}

### simulate panmixia
create_panmixia=function(afmeans,rd){
  afpan=apply(rd,2,function(rd_i){
    rbinom(length(rd_i),rd_i,afmeans)/rd_i
  }
  )
}

## Recode allele freqs from 'ref-major-minor' to 'ref-alt'
recode_AF=function(df,afcols=NULL,coefcols=NULL){
  dropMe= (df$ref != df$major) & (df$ref != df$minor)
  df=df[!dropMe,]
  flipMe=df$ref == df$minor
  if(!is.null(afcols)){
    df[flipMe,afcols]=1-df[flipMe,afcols,with=F]
  }
  if(!is.null(coefcols)){
    df[flipMe,coefcols]=-1*df[flipMe,coefcols,with=F]
  }
  colnames(df)[grep("minor",colnames(df))]="alt"
  df$alt[flipMe]=df$major[flipMe]
  df=df[,-grep("major",colnames(df)),with=F]
  return(df)
}

## get enrichment of significant sites (pval< pThresh) in each 250-SNP window
## for each column of df.glm with column header "p.*"
get_win_enrichment=function(df.glm,meanaf,thresh.effsize=.02,pThresh=.05,nSNPs=250,stepSize=50,ncores=detectCores()-5){
  
  comparisons=gsub("coef.","",colnames(df.glm)[grep("coef.._.",colnames(df.glm))])
  ESmask=effect_size(meanaf,comparisons,treatment="E")>thresh.effsize
  
  divSig=df.glm[,colnames(df.glm)%in%paste0("p.",comparisons)]<=pThresh & ESmask
  divSig[is.na(divSig)]=FALSE
  divSucc=colSums(divSig);divFail=colSums(!divSig);
  divConc=apply(df.glm[,colnames(df.glm)%in%paste0("coef.",comparisons)],2,function(x){(x*df.glm$coef.cli)>0})
  
  chromStarts=aggregate(1:nrow(df.glm),by=list(df.glm$chrom),FUN=min)$x
  chromEnds=c(chromStarts[-1]-1,nrow(df.glm))
  winStarts=unlist(mapply(function(chrStart,chrStop){seq(chrStart,chrStop-nSNPs,stepSize)},chromStarts,chromEnds))
  
  cliSig=!is.na(df.glm$p.cli) & df.glm$p.cli<=pThresh; cliSucc=sum(cliSig);cliFail=sum(!cliSig)
  seaSig=!is.na(df.glm$p.sea) & df.glm$p.sea<=pThresh; seaSucc=sum(seaSig);seaFail=sum(!seaSig)
  cat(cliSucc," sites with clinal pval<",pThresh,"\n")
  cat(seaSucc," sites with seasonal pval<",pThresh,"\n")
  
  sites=df.glm[,c("chrom","pos")]
  rm(ESmask)
  
  system.time({
    df.winEnrich=do.call(rbind,mclapply(1:length(winStarts),function(ww){
        if(ww%%100==0){cat(ww-1,"windows analysed of ",length(winStarts),"\n")}
      winStart=winStarts[ww]
      mask=winStart:(winStart+nSNPs-1);
      nCli=sum(cliSig[mask])
      nSea=sum(seaSig[mask])
      clinalEnrich=phyper(nCli,cliSucc,cliFail,nSNPs,lower.tail=FALSE)
      seaEnrich=phyper(nSea,seaSucc,seaFail,nSNPs,lower.tail=FALSE)
      divEnrich=foreach(ii=1:length(comparisons),.combine=rbind)%do%{
        nDiv=sum(divSig[,ii][mask])
        c(nDiv,phyper(nDiv,divSucc[ii],divFail[ii],nSNPs,lower.tail=FALSE)) #phyper(sampSuccess,popSuccess,popFail,sampSize,lower.tail=FALSE)
      }
      data.frame(ww,sigType=c("clinal","seasonal",comparisons),sigCt=c(nCli,nSea,divEnrich[,1]),enrichment=c(clinalEnrich,seaEnrich,divEnrich[,2]))
  },mc.cores = ncores))})
  df.winEnrich=cbind(sites[winStarts[df.winEnrich$ww],],sites$pos[winStarts[df.winEnrich$ww]+(nSNPs-1)],df.winEnrich)
  colnames(df.winEnrich)[2:3]=c("winStart","winStop")
  return(df.winEnrich)
}

match_sites=function(sitesToMatch,by){
  siteCounts=aggregate(as.numeric(sitesToMatch),by=by,FUN=sum)
  matchedSitesIX=do.call(c,lapply(1:nrow(siteCounts),function(ii){
    matched=Reduce("&",lapply(grep("Group",colnames(siteCounts)),function(gg){by[[gg]]==siteCounts[[gg]][ii]}))
    sample(which(!sitesToMatch & matched),siteCounts$x[ii])
    }))
  matchedSitesMask=rep(FALSE,length(sitesToMatch))
  matchedSitesMask[matchedSitesIX]=TRUE
  return(matchedSitesMask)
}


###############
### Statistics
################

## get standard error
se <- function(x) sqrt(var(x)/length(x))

## calculate Fst between two (equi-length vectors of) allele freqs
Fst=function(p1, p2) {
  fhat <- p1/2 + p2/2  # avg freq across both pops ie (p1+p2)/2
  Htot <- 2*fhat*(1-fhat) #heterozygosity of the avg freq ie 2pq 
  Hwith <- p1*(1-p1) + p2*(1-p2) #avg heterozygosity of indivdual pop freqs  ie. (2pq + 2pq)/2
  fst=(Htot-Hwith)/Htot # how different are they? scaled by total
  fst[fhat==0 | fhat==1]=0 # set fst for fixed sites to zero
  return(fst)
}

## calculate Fst between all pairs of columns in an allele freqs matrix
Fst.mat=function (x) {
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  
  ncy <- ncx <- ncol(x)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(i)) {
      x2 <- x[, i]
      y2 <- x[, j]
      r[i, j] <- mean(Fst( x2, y2),na.rm=T)
    }
  }
  r <- r + t(r) - diag(diag(r))
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  r
}

## calculate Fst between more than two (equi-length vectors of) allele freqs
Fstmulti <- function(pAll) { #pAll should be a matrix with a separate column for each population
  fhat <- rowMeans(pAll)  # avg freq across pops 
  Htot <- 2*fhat*(1-fhat) #heterozygosity of the avg freq ie 2pq 
  Hwith <- rowMeans(2*pAll*(1-pAll)) #avg heterozygosity of indivdual pop freqs  
  (Htot-Hwith)/Htot # how different are they? scaled by total
}

## get zscore
zscore <- function(x,m=mean(x),s=sd(x)){
  z<- (x - mean(x)) / sd(x)
  return(z)
}

## calculate Heterozygosity from an allele frequency vector
hetz <- function(p) { 2*p*(1-p) }

## function to get drift-effective population size, following the method of Reed et al https://www.genetics.org/content/197/2/781.full
##--> based on expectation that the variance in AF over time is proportional to the frequency it started at and the number of gens since starting
get_drift_Ne=function(pInit,pEv,gen){
  
  flipMask=pInit>.5
  pInit[flipMask]=1-pInit[flipMask]
  pEv[flipMask]=1-pEv[flipMask]
  
  snpAFBinID=cut(pInit,breaks = seq(0,.5,.05),labels = FALSE)
  do.call(rbind,
          lapply(1:10,function(bin){
            mask=snpAFBinID==bin & pInit>.01 & pEv>.01 & pEv<.99
            nSites=sum(mask)
            if(nSites==0){
              return(data.frame(bin,initialMAF=NA,nSites,medianShift=NA,Ne=NA))
            } else{
              medianShift=median(abs(pEv[mask]-pInit[mask]))
              sqDev=(pEv[mask] - pInit[mask])^2/(pEv[mask]*(1-pEv[mask]))
              Ne=-gen/(2*log(1-1/nSites)*sum(sqDev)) 
              return(data.frame(bin,initialMAF=median(pInit[mask]),nSites,medianShift,Ne))
            }
          }))
}

### do PCA for sites on each chromosome and for genome
do_pca_perChrom=function(afmat,sites,centerMe=TRUE,scaleMe=TRUE,siteMask=NULL,chroms=NULL){
    if(is.null(chroms)){chroms=c("2L","2R","3L","3R","X","genome")}
    if(is.null(siteMask)){siteMask=rep(TRUE,nrow(sites))}
    list.pca=lapply(chroms,function(chrom){
      chrMask=if(chrom=="genome"){rep(TRUE,nrow(sites))} else{sites$chrom==chrom}
      return(prcomp(t(afmat[siteMask & chrMask,]),center = centerMe,scale=scaleMe))
    })
    names(list.pca)=chroms
    return(list.pca)
}


do_Fst_perChrom=function(afmat,sites,siteMask=NULL,chroms=NULL){
    if(is.null(chroms)){chroms=c("2L","2R","3L","3R","X","genome")}
    if(is.null(siteMask)){siteMask=rep(TRUE,nrow(sites))}
    list.fst=lapply(chroms,function(chrom){
      chrMask=if(chrom=="genome"){rep(TRUE,nrow(sites))} else{sites$chrom==chrom}
      return(Fst.mat(afmat[siteMask & chrMask,]))
    })
    names(list.fst)=chroms
    return(list.fst)
}


## function to fit a generalzied linear model to allele frequencies at each site in siteIX (assumes afmat is already loaded)
## calculates the coefficient and pvalue for the contribution to the mdoel of each sample feature specified in model.vars
## can choose whether to fit model using 'glm' and quasibinomial error model, or 'glmer' with binomial error model, but the ability to specify some variables as random effects
## when using glmer can specify one or more sample features to be random effect variables (usually 'cage')
## specify the name of a sample feature in 'cmpAll' to convert that feature to a factor( if not already) and calculate coeffs/pvals for every pairwise comparison of factor levels
## specify the name of a sample feature in 'dont report' to omit coeffs/pvals from that feature from the returned results
fit_GLM=function(siteIX,sampIX,samps,model.vars,poolCt=100,glmType="glm",randEff=NULL,cmpAll=NULL,dontReport=NULL){
  df=samps[sampIX,colnames(samps)%in%model.vars,with=F];
  if(!is.null(cmpAll)){df[[cmpAll]]=factor(df[[cmpAll]])}
  
  formulaString=paste0(colnames(df),collapse=" + ")
  formulaString=paste0("cts ~ ",formulaString," + 1")
  
  if(glmType=="glmer"){
    sapply(randEff,function(rE){
      formulaString=gsub(rE,paste0('(1|',rE,')'),formulaString)
    })
  }
  cat("Model Forumla is: \n",formulaString)
  
  mclapply(siteIX,function(ix){
    
    Neff=((poolCt*2*rd[ix,sampIX])-1)/(poolCt*2+rd[ix,sampIX])
    cts=cbind(round(Neff*afmat[ix,sampIX]),round(Neff*(1-afmat[ix,sampIX])))
    
    model=switch(glmType,
                 glm=glm(as.formula(formulaString),family="quasibinomial",data=df),
                 glmer=glmer(as.formula(formulaString),family="binomial",nAGQ = 0,data=df)
    )
    if(!is.null(cmpAll)){
      model.multcomp=summary(eval(parse(text=paste0("glht(model, mcp(",cmpAll,"='Tukey'))")))) 
      cp=cbind(coefficients(model.multcomp),model.multcomp$test$pvalues)
      row.names(cp)=gsub(" - ","_",row.names(cp))
    } else{cp=summary(model)$coefficients[-1,c(1,4)]}
    colnames(cp)=c("coef","p")
    cp=cp[grep(dontReport,row.names(cp),invert = TRUE),,drop=F]
    results=c(cp[,1],cp[,2]);names(results)=c(paste0("coef.",row.names(cp)),paste0("p.",row.names(cp)))
    return(results)
  },mc.preschedule=TRUE)
  
}

set_up_sampData=function(sampIX,samps,model.vars,cmpAll){
  sampData=samps[sampIX,colnames(samps)%in%model.vars,drop=FALSE];
  sampData=sampData[,apply(sampData,2,function(x){length(unique(x))>1}),drop=F]
  if(!is.null(cmpAll)){sampData[[cmpAll]]=factor(sampData[[cmpAll]])}
  return(sampData)
}  

extract_coef_pval=function(model,cmpAll=NULL,dontReport=NULL){
  if(!is.null(cmpAll)){
    model.multcomp=summary(eval(parse(text=paste0("glht(model, mcp(",cmpAll,"='Tukey'))")))) 
    cp=cbind(coefficients(model.multcomp),model.multcomp$test$pvalues)
    row.names(cp)=gsub(" - ","_",row.names(cp))
  } else{cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];}
  if(!is.null(dontReport)){cp=cp[grep(paste0("(",paste0(dontReport,collapse="|"),")"),row.names(cp),invert = TRUE),,drop=FALSE];}
  return(cp)
}

fit_GLM_one = function(af.site,rd.site,sampData,formulaString,cmpAll=NULL,dontReport=NULL){
  Neff=((poolCt*2*rd.site)-1)/(poolCt*2+rd.site);
  cts=cbind(round(Neff*af.site),round(Neff*(1-af.site)))
  model=glm(as.formula(formulaString),family="quasibinomial",data=sampData) 
  cp=extract_coef_pval(model,cmpAll,dontReport)      
}

fit_GLM_all=function(siteIX,sampIX,samps,model.vars,poolCt=100,cmpAll=NULL,dontReport=NULL){
  
  sampData=set_up_sampData(sampIX,samps,model.vars,cmpAll)
  
  formulaString=paste0(colnames(sampData),collapse=" + ")
  formulaString=paste0("cts ~ ",formulaString," + 1")
  
  cat("Model Formula is: \n",formulaString,"\n")
  nSites=length(siteIX);cat("there are ",nSites," sites\n")
  do.call(rbind,mclapply(siteIX,function(ix){
    
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    
    cp=fit_GLM_one(afmat[ix,sampIX],rd$EC[rd$chrom==sites$chrom[ix]],sampData,formulaString,cmpAll,dontReport)
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.",row.names(cp)),paste0("p.",row.names(cp)))
    return(results)},mc.cores=ncores))
}

###################
## Annotations
##################
write_vcf=function(df,filename){
  colnames(df)[match(c("chrom","pos","ref","alt"),colnames(df))]=c("#CHROM","POS","REF","ALT")
  df$ID=1:nrow(df)
  df$QUAL=NA
  df$FILTER="PASS"
  df$INFO=NA
  cat("##fileformat=VCFv4.1\n",file = filename)
  suppressWarnings(write.table(df[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"),with=F],sep="\t",file=filename,row.names=F,quote=F,append = TRUE))
}

annotate_vcf=function(vcf,mergeByPos=TRUE){
  system(paste0(annotationScript[1]," ",vcf," ",gffFeaturesList[1],' "',gffFiles[1],'"'));
  ann=fread(gsub(".vcf$",".genes",vcf))
  if(mergeByPos){
    geneRows=ann$feature=="gene"
    ann.genenameagg=aggregate(ann$Name[geneRows],by=list(ann$chrom[geneRows],ann$pos[geneRows]),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann.genefullagg=aggregate(ann$fullname[geneRows],by=list(ann$chrom[geneRows],ann$pos[geneRows]),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann.flybaseagg=aggregate(ann$ID[geneRows],by=list(ann$chrom[geneRows],ann$pos[geneRows]),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann.featureagg=aggregate(paste(ann$Name,ann$feature,sep="|"),by=list(ann$chrom,ann$pos),FUN=function(x){paste(unique(x),collapse=";")})
    ann.entrezagg=aggregate(ann$EntrezGene,by=list(ann$chrom,ann$pos),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann=suppressWarnings(Reduce(function(x, y) merge(x, y, all=TRUE,by=c("Group.1","Group.2")), 
                                list(ann.genenameagg, ann.genefullagg,ann.flybaseagg,ann.featureagg,ann.entrezagg)
    ))
    colnames(ann)=c("chrom","pos","name","fullname","fbID","features","entrezID")
  }
  return(ann[order(ann$chrom,as.numeric(ann$pos)),])
}

add_snpEff=function(df,fields="effect"){
  snpEff=fread(snpEffFile[1])
  df=cbind(df,snpEff[match(paste0(df$chrom,df$pos),paste0(snpEff$chrom,snpEff$pos)),fields,with=F])
  return(df)
}

add_goTerms=function(df){
  load(gotermData[1])
  goMask=df.go$fbID%in%unlist(strsplit(df$fbID,";"))
  df.go=aggregate(df.go$goterm[goMask],by=list(df.go$fbID[goMask]),FUN=function(x){paste(x,collapse="|")})
  df$goterms=sapply(df$fbID,function(x){goMask=df.go$Group.1%in%strsplit(x,";");paste(df.go$x[goMask],collapse="|")})
  return(df)
}
add_clinal=function(df,clinalityFile="/mnt/cages/ref_data/heather_sea_cli/mel_clinal_spring.glm.noheader"){
  df.cli=fread(clinalityFile[1])[,1:4,with=F]
  colnames(df.cli)=c("chrom","pos","coef.cli","p.cli")
  df=merge(df,df.cli,by=c("chrom","pos"),all.x=T)
  return(df)
}
add_seasonal=function(df,seasonalFile="/mnt/cages/ref_data/heather_sea_cli/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm"){
  df.sea=fread(seasonalFile[1])[,1:4,with=F]
  colnames(df.sea)=c("chrom","pos","coef.sea","p.sea")
  df=merge(df,df.sea,by=c("chrom","pos"),all.x=T)
  return(df)
}
mark_centro_telo_meres=function(df){  
  #from Comeron https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002905#s4
  #alternative: see Kauer et al 2003 http://www.genetics.org/content/165/3/1137 : loci with recombination rates <0.0001% recombination per kilobase after adjusting for zero recombination in males (i.e., multiplying by 0.67 for the X chromosome and by 0.5 for the third chromosome)
  chroms=c("2L","2R","3L","3R","X")
  starts=c(.5,5.2,0.7,9.7,2.3)*1000000
  ends=c(17.4,20.8,19.9,26.9,20.8)*1000000
  
  mask=do.call(c,mapply(function(chr,s,e){which(df$chrom==chr & (df$pos<s | df$pos>e))},chroms,starts,ends))
  return(seq_along(df$pos) %in% mask)
}

################
## Plots
################
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 5, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

gg_color_hue=function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_site=function(af,samps){
  samps$af=af
  ggplot(samps,aes(x=tpt,y=af,color=cage)) + geom_point() + facet_wrap(facets=vars(treatment),scales = "free_x") + 
    theme_minimal()
}

plot_MDS=function(df.fst,samps,treatments,tpts,drawArrows=FALSE,mytitle){
  sampMask=samps$treatment%in%treatments & samps$tpt%in%tpts
  fit<-cmdscale(as.matrix(df.fst),k=2,eig=TRUE)
  
  df=cbind(as.data.frame(fit$points[sampMask,]),samps[sampMask,])
  df$group=paste0(df$treatment,df$cage)
  
  df2=df[order(df$treatment,df$cage,df$tpt),]
  df2$V1.end=c(df2$V1[-1],0)
  df2$V2.end=c(df2$V2[-1],0)
  df2=df2[df2$tpt<4,]
  df2[df2$treatment=="FOUND",c("V1.end","V2.end")]=df2[df2$treatment=="FOUND",c("V1","V2")]
  df2$shift.y=df2$V2.end-df2$V2
  df2$shift.x=df2$V1.end-df2$V1
  
  ### color by treatment,alpha by timept, lines within cage
  ggMDS<-ggplot(df,aes(V1,V2,label=sampID,color=treatment)) + 
    labs(x="MDS1", y="MDS2", color="Treatment", alpha="Timepoint") + 
    geom_point(size=2,aes(alpha=tpt)) +
    ggtitle(mytitle) +
    scale_alpha(range = c(0.4, 1),limits=c(1,4),breaks=sort(unique(df$tpt)),labels=sort(unique(df$tpt))) +
    scale_color_manual(labels=c('Founder' ="Founder","E"="E","F"="N","ED"="ED","B"="B","S"="S"),
                       values = c('Founder' = "black", 
                                  'E' = gg_color_hue(3)[1],
                                  'F' = gg_color_hue(3)[2],
                                  'ED' = gg_color_hue(3)[3],
                                  'B' = gg_color_hue(5)[2],
                                  'S' = gg_color_hue(5)[5])) + 
    theme_minimal() +  
    theme(plot.title = element_text(size=12, face="bold"),axis.text = element_blank())
  
  if(drawArrows){
    ggMDS<- ggMDS + geom_segment(data=df2,
                                 aes(x=V1,y=V2,xend=V1.end,yend=V2.end,col=treatment,alpha=tpt+1,linetype="Cage"),
                                 arrow=arrow(type = "closed",length = unit(0.3,"cm")),
                                 show.legend = FALSE) +
      scale_linetype_manual(values=c("Cage" = 4))
  }
  
  return(ggMDS)
}

plot_sites_manhattan=function(chrom,pos,score,group=NULL,ylab="score",mytitle="",pointSize=.7,pointShape=16,pointAlpha=.7){
  df=data.frame(chrom,pos,score);
  if(!is.null(group)){df$group=group}
  gg<-ggplot(df,aes(x=pos/1000000,y=score)) + facet_wrap(vars(chrom),scales = "free_x",nrow=1) + 
    theme_minimal() + labs(y=ylab,x="position (Mb)") + ggtitle(mytitle)
  if(!is.null(group)){
    gg<- gg + geom_point(aes(color=group),size=pointSize,shape=pointShape,alpha=pointAlpha)
  } else{
    gg<-gg + geom_point(size=pointSize,shape=pointShape,pointAlpha=.7)
  }
  gg
}

make_window_GLM_plot=function(df,scoreCol,colorCol,myalpha=.7,myshape=1){
  ### plots scores for window-based enrichment as a point in the center of the window
  df$score=df[[scoreCol]]
  df$color=df[[colorCol]]
  gg<-ggplot(df,aes((winStart+winStop)/2000000,score)) + 
    geom_point(aes(color=color),alpha=myalpha,shape=myshape) + 
    theme_minimal() + facet_wrap(~chrom,nrow = 1,scales="free_x") + 
    labs(x="genomic window (Mb)",y=expression("-log"[10] * "(p"["enrichment"] *")"))
  gg
}





###MCB additions

calc_Neff=function(rd,poolCt){
  ((poolCt*2*rd)-1)/(poolCt*2+rd)
}


fit_GLM=function(afMatrix,rdMatrix,sampleData,model.vars,poolCt=100,ncores){
  ## afMatrix: matrix of allele frequencies; rows are sites, columns are samples
  ## rdMatrix: matrix of read depths corresponding to allele freqs in afMatrix
  ## sampleData: data frame of metadata about each sample; rows are samples, columns are metadata items - samples be in the same order as columns in afMatrix and rdMatrix
  ## model.vars: vector of strings with names of the columns of sampleData to be included in the model
  ## poolCt: number of individuals that were pooled to create each sample
## this function requires the library doMC

  registerDoMC(ncores)
  ## set up sample data
  df=as.data.frame(sampleData[,colnames(sampleData)%in%model.vars])  # subsets the sample data object to just the variables to be included in the model
  colnames(df) <- model.vars

  ## set up formula for GLM
  formulaString=paste0(colnames(df),collapse=" + ")  ## assumes all variables are additive - ie no interactions
  formulaString=paste0("cts ~ ",formulaString," + 1")  ## add a constant
  cat("Model Formula is: \n",formulaString,"\n")

  ## run GLM over all sites in parallel
  Neff=calc_Neff(rdMatrix,poolCt)  ## see calc_Neff function above
  do.call(rbind,mclapply(1:nrow(afMatrix),function(ix){  ## mclapply is another way to run things in parallel
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    cts=cbind(t(round(Neff[ix,]*afMatrix[ix,])),t(round(Neff[ix,]*(1-afMatrix[ix,]))))  ## makes a 2 column matrix of read counts for alt reads and ref reads for each sample
    df$cts=cts
    df=as.data.frame(df)
    model=glm(formulaString,family="quasibinomial",data=df)
    cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];  ## grab the model coefficients and p-values
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.", TPA, "_", TPB),paste0("p.",TPA, "_", TPB))
    return(results)
  },mc.cores=ncores))}


fit_GLM_ContinuousTime=function(afMatrix,rdMatrix,sampleData, vec, model.vars,poolCt=100,ncores){
  ## afMatrix: matrix of allele frequencies; rows are sites, columns are samples
  ## rdMatrix: matrix of read depths corresponding to allele freqs in afMatrix
  ## sampleData: data frame of metadata about each sample; rows are samples, columns are metadata items - samples be in the same order as columns in afMatrix and rdMatrix
  ## model.vars: vector of strings with names of the columns of sampleData to be included in the model
  ## poolCt: number of individuals that were pooled to create each sample
## this function requires the library doMC

  registerDoMC(ncores)
  ## set up sample data
  df=as.data.frame(sampleData[,colnames(sampleData)%in%model.vars])  # subsets the sample data object to just the variables to be included in the model
  colnames(df) <- model.vars

  ## set up formula for GLM
  formulaString=paste0(colnames(df),collapse=" + ")  ## assumes all variables are additive - ie no interactions
  formulaString=paste0("cts ~ ",formulaString," + 1")  ## add a constant
  cat("Model Formula is: \n",formulaString,"\n")

  ## run GLM over all sites in parallel
  Neff=calc_Neff(rdMatrix,poolCt)  ## see calc_Neff function above
  do.call(rbind,mclapply(1:nrow(afMatrix),function(ix){  ## mclapply is another way to run things in parallel
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    cts=cbind(t(round(Neff[ix,]*afMatrix[ix,])),t(round(Neff[ix,]*(1-afMatrix[ix,]))))  ## makes a 2 column matrix of read counts for alt reads and ref reads for each sample
    df$cts=cts
    df=as.data.frame(df)
    model=glm(formulaString,family="quasibinomial",data=df)
    cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];  ## grab the model coefficients and p-values
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.", vec[1], "_", vec[length(vec)]), paste0("p.",vec[1], "_", vec[length(vec)]))
    return(results)
  },mc.cores=ncores))}


fit_GLM_reps=function(afMatrix,rdMatrix,sampleData,model.vars,poolCt=100,ncores){
  ## afMatrix: matrix of allele frequencies; rows are sites, columns are samples
  ## rdMatrix: matrix of read depths corresponding to allele freqs in afMatrix
  ## sampleData: data frame of metadata about each sample; rows are samples, columns are metadata items - samples be in the same order as columns in afMatrix and rdMatrix
  ## model.vars: vector of strings with names of the columns of sampleData to be included in the model
  ## poolCt: number of individuals that were pooled to create each sample
## this function requires the library doMC

  registerDoMC(ncores)
  ## set up sample data
  df=as.data.frame(sampleData[,colnames(sampleData)%in%model.vars])  # subsets the sample data object to just the variables to be included in the model
  colnames(df) <- model.vars

  ## set up formula for GLM
  formulaString=paste0(colnames(df),collapse=" + ")  ## assumes all variables are additive - ie no interactions
  formulaString=paste0("cts ~ ",formulaString," + 1")  ## add a constant
  cat("Model Formula is: \n",formulaString,"\n")

  ## run GLM over all sites in parallel
  Neff=calc_Neff(rdMatrix,poolCt)  ## see calc_Neff function above
  do.call(rbind,mclapply(1:nrow(afMatrix),function(ix){  ## mclapply is another way to run things in parallel
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    cts=cbind(t(round(Neff[ix,]*afMatrix[ix,])),t(round(Neff[ix,]*(1-afMatrix[ix,]))))  ## makes a 2 column matrix of read counts for alt reads and ref reads for each sample
    df$cts=cts
    df=as.data.frame(df)
    model=glm(formulaString,family="quasibinomial",data=df)
    cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];  ## grab the model coefficients and p-values
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.biolreps.tp", tp),paste0("p..biolreps.tp",tp))
    return(results)
  },mc.cores=ncores))}


fit_GLM_Side=function(afMatrix,rdMatrix,sampleData,model.vars,poolCt=100,ncores){
  ## afMatrix: matrix of allele frequencies; rows are sites, columns are samples
  ## rdMatrix: matrix of read depths corresponding to allele freqs in afMatrix
  ## sampleData: data frame of metadata about each sample; rows are samples, columns are metadata items - samples be in the same order as columns in afMatrix and rdMatrix
  ## model.vars: vector of strings with names of the columns of sampleData to be included in the model
  ## poolCt: number of individuals that were pooled to create each sample
## this function requires the library doMC

  registerDoMC(ncores)
  ## set up sample data
  df=as.data.frame(sampleData[,colnames(sampleData)%in%model.vars])  # subsets the sample data object to just the variables to be included in the model
  colnames(df) <- model.vars

  ## set up formula for GLM
  formulaString=paste0(colnames(df),collapse=" + ")  ## assumes all variables are additive - ie no interactions
  formulaString=paste0("cts ~ ",formulaString," + 1")  ## add a constant
  cat("Model Formula is: \n",formulaString,"\n")

  ## run GLM over all sites in parallel
  Neff=calc_Neff(rdMatrix,poolCt)  ## see calc_Neff function above
  do.call(rbind,mclapply(1:nrow(afMatrix),function(ix){  ## mclapply is another way to run things in parallel
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    cts=cbind(t(round(Neff[ix,]*afMatrix[ix,])),t(round(Neff[ix,]*(1-afMatrix[ix,]))))  ## makes a 2 column matrix of read counts for alt reads and ref reads for each sample
    df$cts=cts
    df=as.data.frame(df)
    model=glm(formulaString,family="quasibinomial",data=df)
    cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];  ## grab the model coefficients and p-values
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.tpt", tp),paste0("p.tpt",tp))
    return(results)
  },mc.cores=ncores))}




############function to get allele frequency shifts in subsequent time segments - used for Resolving Evolutionary Phases Orchard 2021
get.NextSegShifts = function(x, base.data, sites, af.shifts){
    chr = x$chrom
    p = x$pos
    sign = x$sign.shift
    seg = as.character(x$comparison)
    next.seg = strsplit(seg, split = "_")[[1]][2]
    next.seg = paste0(next.seg, "_", as.numeric(next.seg) + 1)
    min.freq = x$af.base.mean - 0.05
    max.freq = x$af.base.mean + 0.05
    pos.matched = base.data %>% filter(chrom == chr & pos != p & between(af.base.mean, min.freq, max.freq)) %>%
                                sample_n(1) %>% pull(pos)         
    matchedShift.NextSeg = cbind(sites, af.shifts) %>% filter(chrom == chr & pos == pos.matched) %>% pull(paste0("dAF.", next.seg))
    afShift.NextSeg = (cbind(sites, af.shifts) %>% filter(chrom == chr & pos == p) %>% pull(paste0("dAF.", next.seg))) * sign
    matchedShift.NextSeg = cbind(sites, af.shifts) %>% filter(chrom == chr & pos == pos.matched) %>% pull(paste0("dAF.", next.seg))
    x = cbind(x, pos.matched, next.seg, afShift.NextSeg, matchedShift.NextSeg)
    return(x)
}



###Code to take in sig sites for a specific time segment, produce dataframe of afShifts in subsequent windows and also run stats comparing target sites to matched controls
get.NextSegShifts.ConcatData.Stats = function(RData, baseData, afShifts, sigSites, outputFile){
    load(RData)
    load(baseData)
    base.data = cbind(sites.base, afmat.base.mean)
    af.shifts = read.csv(afShifts)
    df.sig = read.csv(sigSites)
    seg = strsplit(sigSites, "ResolvePhases.")[[1]][2]
    seg = strsplit(seg, "Start")[[1]][1]
    df = df.sig %>%
        dplyr::select(chrom, pos, afShift, comparison) %>%
        mutate(sign.shift = sign(afShift),
              phased.afShift = abs(afShift)) 
    df = left_join(df, base.data)
    df.new = data.frame()
    for (i in 1:nrow(df)){
        vec = get.NextSegShifts(df[i,], base.data, sites, af.shifts)
        df.new = rbind(df.new, vec)
    }
    write.csv(df.new, paste0("./", outputFile, "csv"))
    
    stats = data.frame()
    comparisons = unique(df.new$comparison)
    for (comp in comparisons){
        df.comp = df.new %>% filter(comparison == comp)
        segment = comp
        test.seg = strsplit(segment, split = "_")[[1]][2]
        test.seg = paste0(test.seg, "_", as.numeric(test.seg) + 1)
        medianShift.target = median(df.comp$phased.afShift)
        medianShift.target.nextseg = median(df.comp$afShift.NextSeg)
        medianShift.matched.nextseg = median(df.comp$matchedShift.NextSeg)
        pval = t.test(df.comp$afShift.NextSeg, df.comp$matchedShift.NextSeg)$p.val
        data = cbind(segment, test.seg, medianShift.target, medianShift.target.nextseg, medianShift.matched.nextseg, pval)
        stats = rbind(stats, data)   
    }
    write.csv(stats, paste0("./", outputFile, "Stats.csv"))
}




##################################Same as above, but with new parameter that allows variable next seg length#############################
get.NextSegShifts.ConcatData.Stats.V2 = function(RData, baseData, afShifts, sigSites, outputFile, next.seg.length){
    load(RData)
    load(baseData)
    base.data = cbind(sites.base, afmat.base.mean)
    af.shifts = read.csv(afShifts)
    df.sig = read.csv(sigSites)
    seg = strsplit(sigSites, "ResolvePhases.")[[1]][2]
    seg = strsplit(seg, "Start")[[1]][1]
    df = df.sig %>%
        dplyr::select(chrom, pos, afShift, comparison) %>%
        mutate(sign.shift = sign(afShift),
              phased.afShift = abs(afShift))
    df = left_join(df, base.data)
    df.new = data.frame()
    for (i in 1:nrow(df)){
        vec = get.NextSegShifts.V2(df[i,], base.data, sites, af.shifts, next.seg.length)
        df.new = rbind(df.new, vec)
    }
    write.csv(df.new, paste0("./", outputFile, ".csv"))

    stats = data.frame()
    comparisons = unique(df.new$comparison)
    for (comp in comparisons){
        df.comp = df.new %>% filter(comparison == comp)
        segment = comp
        test.seg = strsplit(segment, split = "_")[[1]][2]
        test.seg = paste0(test.seg, "_", as.numeric(test.seg) + next.seg.length)
        medianShift.target = median(df.comp$phased.afShift)
        medianShift.target.nextseg = median(df.comp$afShift.NextSeg)
        medianShift.matched.nextseg = median(df.comp$matchedShift.NextSeg)
        pval = t.test(df.comp$afShift.NextSeg, df.comp$matchedShift.NextSeg)$p.val
        data = cbind(segment, test.seg, medianShift.target, medianShift.target.nextseg, medianShift.matched.nextseg, pval)
        stats = rbind(stats, data)
    }
    write.csv(stats, paste0("./", outputFile, "Stats.csv"))
}


get.NextSegShifts.V2 = function(x, base.data, sites, af.shifts, next.seg.length){
    chr = x$chrom
    p = x$pos
    sign = x$sign.shift
    seg = as.character(x$comparison)
    next.seg = strsplit(seg, split = "_")[[1]][2]
    next.seg = paste0(next.seg, "_", as.numeric(next.seg) + next.seg.length)
    min.freq = x$af.base.mean - 0.05
    max.freq = x$af.base.mean + 0.05
    pos.matched = base.data %>% filter(chrom == chr & pos != p & between(af.base.mean, min.freq, max.freq)) %>%
                                sample_n(1) %>% pull(pos)         
    matchedShift.NextSeg = cbind(sites, af.shifts) %>% filter(chrom == chr & pos == pos.matched) %>% pull(paste0("dAF.", next.seg))
    afShift.NextSeg = (cbind(sites, af.shifts) %>% filter(chrom == chr & pos == p) %>% pull(paste0("dAF.", next.seg))) * sign
    matchedShift.NextSeg = cbind(sites, af.shifts) %>% filter(chrom == chr & pos == pos.matched) %>% pull(paste0("dAF.", next.seg))
    x = cbind(x, pos.matched, next.seg, afShift.NextSeg, matchedShift.NextSeg)
    return(x)
}


###################Get Next Seg Shifts for LOA
get.NextSegShifts.loa = function(df.sig.matched, Comp, NextSeg){
    df = df.sig.matched %>% filter(comparison == Comp)
    df.shifts = data.frame()
    for (i in 1:nrow(df)){
        chr = df[i, "chrom"]
        p = df[i, "pos"]
        pos.matched = df[i, "pos.matched"]
        sign = df[i, "sign.shift"]
        seg = Comp
        next.seg = NextSeg
        matchedShift.NextSeg = cbind(sites, afShifts.loc) %>% filter(chrom == chr & pos == pos.matched) %>% pull(paste0("dAF.", next.seg))
        afShift.NextSeg = (cbind(sites, afShifts.loc) %>% filter(chrom == chr & pos == p) %>% pull(paste0("dAF.", next.seg))) * sign
        x = cbind((df[i,] %>% dplyr::select(chrom, pos, comparison, sigLevel, sign.shift, phased.afShift, af.base.mean, pos.matched)),
                 next.seg, afShift.NextSeg, matchedShift.NextSeg)
        df.shifts = rbind(df.shifts, x)
    }
    return(df.shifts)
}


#############Get Matched SNP LOA

get.MatchedSNP.loa = function(x, base.data){
    chr = x$chrom
    p = x$pos
    seg = as.character(x$comparison)
    min.freq = x$af.base.mean - 0.05
    max.freq = x$af.base.mean + 0.05
    pos.matched = base.data %>% filter(chrom == chr & pos != p & between(af.base.mean, min.freq, max.freq)) %>%
                                sample_n(1) %>% pull(pos)
    x = cbind(x, pos.matched)
    return(x)
    }

#####New version of Get Next Seg Shifts that does not loop through each row, but just selects all pertinent rows from af.shifts dataframe simultaneously and binds to df.sig
get.NextSegShifts.V4 = function(df.sig.matched, RData, Comp, next.segs){
    load(RData)
    df = df.sig.matched %>% filter(comparison == Comp) %>%
        mutate(sign.shift = sign(afShift), phased.afShift = abs(afShift))
    af.shifts = get_af_shifts(afmat, samps, cage_set = NULL, next.segs)
    af.shifts = cbind(sites, af.shifts)
    af.shifts.target = af.shifts
    colnames(af.shifts.target) = paste(colnames(af.shifts),"target",sep=".")
    colnames(af.shifts.target)[1:2] = c('chrom', 'pos')
    af.shifts.matched = af.shifts
    colnames(af.shifts.matched) = paste(colnames(af.shifts),"matched",sep=".")
    colnames(af.shifts.matched)[1:2] = c('chrom', 'pos')
    df.target = left_join((df %>% dplyr::select(chrom, pos)), af.shifts.target) %>% dplyr::select(-chrom, -pos)
    signs = df$sign.shift
    df.target.phased = data.frame(nrow = nrow(df))
    for (cols in 1:ncol(df.target)){
        col = df.target[,cols]
        col.new = signs * col
        df.target.phased = cbind(df.target.phased, col.new)
    }
    df.target.phased = as.data.frame(df.target.phased[,-1])
    colnames(df.target.phased) = colnames(df.target)    
    df.matched = df %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)
    df.matched = left_join(df.matched, af.shifts.matched) %>% dplyr::select(-chrom, -pos)
    df.shifts = cbind(df, df.target.phased, df.matched)
    return(df.shifts)
    }



###Get Matched SNP - #########################
########takes in vector containing chrom, pos, and comparison; as well as a base.data file with a column "af.base.mean" indicating mean af acrose baseline samples
get.MatchedSNP = function (x, base.data) 
{
    chr = x$chrom
    p = x$pos
    seg = as.character(x$comparison)
    min.freq = x$af.base.mean - 0.05
    max.freq = x$af.base.mean + 0.05
    pos.matched = base.data %>% filter(chrom == chr & pos != 
        p & between(af.base.mean, min.freq, max.freq)) %>% sample_n(1) %>% 
        pull(pos)
    x = cbind(x, pos.matched)
    return(x)
}


##########################Get Weekly AF shifts -###################
#####takes in a df.sig.matched file (which contains info on comparison, chrom/pos of target and matched site); RData (for sites info); and weekly.shifts dataframe
get.WeeklyShifts = function (df.Sig.Matched, RData, weekly.shifts) {
    load(RData)
    df = df.Sig.Matched %>% mutate(sign.shift = sign(afShift), 
        phased.afShift = abs(afShift))
    af.shifts = read.csv(weekly.shifts)
    af.shifts = cbind(sites, af.shifts)
    af.shifts.target = af.shifts
    colnames(af.shifts.target) = paste(colnames(af.shifts), "target", 
            sep = ".")
    colnames(af.shifts.target)[1:2] = c("chrom", "pos")
    af.shifts.matched = af.shifts
    colnames(af.shifts.matched) = paste(colnames(af.shifts), 
        "matched", sep = ".")
    colnames(af.shifts.matched)[1:2] = c("chrom", "pos")
    df.target = as.data.frame(left_join((df %>% dplyr::select(chrom, pos)), 
        af.shifts.target) %>% dplyr::select(-chrom, -pos))
    df.matched = as.data.frame(left_join((df %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)), 
        af.shifts.matched) %>% dplyr::select(-chrom, -pos))
    signs = df$sign.shift
    df.target.phased = data.frame(nrow = nrow(df))
    for (cols in 1:ncol(df.target)) {
        col = df.target[, cols]
        col.new = signs * as.numeric(col)
        df.target.phased = cbind(df.target.phased, col.new)
        }
        df.target.phased = as.data.frame(df.target.phased[, -1])
        colnames(df.target.phased) = colnames(df.target)
    df.matched.phased = data.frame(nrow = nrow(df))
    for (cols in 1:ncol(df.matched)) {
            col = df.matched[, cols]
            col.new = 1 * col
            df.matched.phased = cbind(df.matched.phased, col.new)
        }
    df.matched.phased = as.data.frame(df.matched.phased[, -1])
    colnames(df.matched.phased) = colnames(df.matched)
    df.shifts = cbind(df, df.target.phased, df.matched.phased)
    return(df.shifts)
    }

############################################Get sig sites and matched sites
Get.Sig.Match = function(glmres, HAFsFile, BaseData){
    load(HAFsFile)
    load(glmres)

    seg.length = strsplit(glmres, '[_]')[[1]][2]
    sites=df.glm %>% dplyr::select(chrom,pos)

    comparisons = grep("coef.", names(df.glm), value = TRUE)  
    for(i in (1:length(comparisons))){
        comparisons[i] = strsplit(comparisons[i], '[.]')[[1]][2]
    }
#####Get Sig Sites
    af.shifts <- get_af_shifts(afmat, samps,cage_set=NULL, comparisons)
    FDR=get_glm_FDR.V2(df.glm)
    df.sig <- get_sig_sites(df.glm, comparisons, FDR, afShifts = af.shifts, fdrThreshs, esThreshs)
    write.csv(df.sig, paste0('df.sig.ExpandingPhases.', seg.length,'.orch21.csv'), row.names = FALSE)
    df.sig = df.sig %>% filter(FDR < 0.05) #filter sig sites by FDR 0.05
    comparisons = as.character(unique(df.sig$comparison)) ##In case any comparisons have no sig sites, re-generate comparisons vector

##############Find matched snp
    load(BaseData) 
    base.data = cbind(sites.base, afmat.base.mean)
    sites.sub = sites %>% mutate(snp = paste0(chrom, "_", pos))
    base.data = base.data %>% mutate(snp = paste0(chrom, "_", pos))
    base.data = base.data %>% filter(snp %in% sites.sub$snp) %>% dplyr::select(-snp)
    head(base.data)
    df.sig = df.sig %>% filter(FDR < 0.05) #filter sig sites by FDR 0.05
    comparisons = as.character(unique(df.sig$comparison)) ##In case any comparisons have no sig sites, re-generate comparisons vector
#############Find matched snps
    load(BaseData) 
    base.data = cbind(sites.base, afmat.base.mean)
    sites.sub = sites %>% mutate(snp = paste0(chrom, "_", pos))
    base.data = base.data %>% mutate(snp = paste0(chrom, "_", pos))
    base.data = base.data %>% filter(snp %in% sites.sub$snp) %>% dplyr::select(-snp)

    df.sig = left_join(df.sig, base.data)
    df.sig.matched = data.frame()
        for (i in 1:nrow(df.sig)){
            vec = get.MatchedSNP.loa(df.sig[i,], base.data)
            df.sig.matched = rbind(df.sig.matched, vec)
        }
    write.csv(df.sig.matched, paste0("df.sig.Matched.ExpandingPhases.", seg.length, ".csv") , row.names = FALSE)

}

#######################################Get shift stats for target and matched control sites
##Takes in df.shifts dataframe produced by get.WeeklyShifts
get.shift.stats = function(df.shifts){
    comps = as.character(unique(df.shifts$comparison))
    df.stat.meta = data.frame()
    for (comp in comps){
        df.c = df.shifts %>% filter(comparison == comp)
        segs = grep(".target", names(df.c), value = TRUE)  
        for(i in (1:length(segs))){
            segs[i] = strsplit(segs[i], '[.]')[[1]][2]
        }
    for (seg in segs){
            df.stat = data.frame()
                ID.Seg = comp
                Test.Seg = seg
                matched.vec = df.c %>% dplyr::select(paste0("dAF.", 
                    seg, ".matched"))
                matched.vec = matched.vec[, 1]
                median.matched = median(matched.vec)
                target.vec = df.c %>% dplyr::select(paste0("dAF.", 
                    seg, ".target"))
                target.vec = target.vec[, 1]
                median.target = median(target.vec)
                if (length(target.vec) > 10) {
                    pvalue = t.test(target.vec, matched.vec)$p.value
                }
                else {
                    pvalue = NA
                }
                df.stat = cbind(ID.Seg, Test.Seg, median.matched, median.target, 
                    pvalue)
             df.stat.meta = rbind(df.stat.meta, df.stat)

        }
     df.stat.meta = rbind(df.stat.meta, df.stat)
    }
        return(df.stat.meta)
}



#######Function for Expanding Phases get weekly shifts/stats
ExpandPhase.Shifts.Stats = function(df.sig.matched, HAFsFile, WeeklyShifts){
################Get shifts at target and matched sites for all weekly windows
    df = read.csv(df.sig.matched)
    comparisons = as.character(unique(df$comparison))
    df.shifts = data.frame()
    for (comp in comparisons){
        df.sig.matched.comp = df %>% filter(comparison == comp)
        df.shifts.comp = get.WeeklyShifts(df.sig.matched.comp, HAFsFile, WeeklyShifts)
        df.shifts = rbind(df.shifts, df.shifts.comp)
    }
    write.csv(df.shifts, paste0("df.shifts.ExpandingPhases.", seg.length, ".csv") , row.names = FALSE)

####Get median shift and stats for each weekly segment
    df.stats = get.shift.stats(df.shifts)
    write.csv(df.stats, paste0("df.stats.ExpandingPhases.", seg.length, ".csv") , row.names = FALSE)
    
}

#############################Calculated expected covereage (EEC) from haf-pipe generated data##################################

calc_expected_ec=function(rd,gen,pct_missing,nof_snps,chrom_length,recomb_rate){

#example input: rd=5,gen=20,pct_missing=2,nof_snps=283438
# rd: actual read depth
# gen: nof generations since population founding
# pct_missing: percent of founder genotype calls that are missing (%, not a fraction)
# nof_snps: average number of snps per chromosome
# chrom_length: average length of a chromosome
# recomb_rate: average recombination rate

    q=18
    mycoeffs=data.frame(a=0.5199118,b=-0.6909052,c=0.3553630)
    winSize=round(qexp(q/100,1/((chrom_length)/((recomb_rate)*(chrom_length)*(gen)+1)))/1000)
    nReadsPerWin=(rd)*(nof_snps)*(winSize)*1000/(chrom_length)
    ec = 10^(mycoeffs$a * log10( nReadsPerWin ) + mycoeffs$b * log10(1+pct_missing) + mycoeffs$c )
    ec
    }

######################



##################Get shift stats, combining all comparisons in df.shifts dataframes
get.shift.stats.ComparisonsCombined = function(df.shifts){
    df.stat.meta = data.frame()
    segs = grep(".target", names(df.shifts), value = TRUE)
        for (i in (1:length(segs))) {
            segs[i] = strsplit(segs[i], "[.]")[[1]][2]
        }
    for (seg in segs) {
                df.stat = data.frame()
                ID.Seg = 'AllSigSites'
                Test.Seg = seg
                matched.vec = df.shifts %>% dplyr::select(paste0("dAF.", 
                    seg, ".matched"))
                matched.vec = matched.vec[, 1]
                median.matched = median(matched.vec)
                target.vec = df.shifts %>% dplyr::select(paste0("dAF.", 
                    seg, ".target"))
                target.vec = target.vec[, 1]
                median.target = median(target.vec)
                if (length(target.vec) > 10) {
                    pvalue = t.test(target.vec, matched.vec)$p.value
                }
                else {
                    pvalue = NA
                }
                df.stat = cbind(ID.Seg, Test.Seg, median.matched, 
                    median.target, pvalue)
                df.stat.meta = rbind(df.stat.meta, df.stat)
            }
            df.stat.meta = rbind(df.stat.meta, df.stat)
        }

########################################

#######Shuffle Sites - using a df.base with starting frequency percentile info
######
shuffle.sites = function(df.base, d.sites, ncores){
    registerDoMC(ncores)
    do.call(rbind, mclapply(1:nrow(d.sites), function(ix) {
        if (ix%%1000 == 0) {
            cat("working on site ", ix, "\n")
        }
	s = d.sites[ix,]$snp
        c = as.character((df.base %>% filter(snp == s))$chrom)
        perc = as.character((df.base %>% filter(snp == s))$percentile)
        d.site.shuff = df.base %>% filter(chrom == c) %>% filter(snp != s & percentile == perc) %>% sample_n(1)
        return(d.site.shuff)
    }, mc.cores = ncores))  
}
##################
#################

###Shuffle sites by ensuring they are within 2 percent of starting frequency

shuffle.sites.V2 = function (df.base, d.sites, ncores) 
{
    registerDoMC(ncores)
    do.call(rbind, mclapply(1:nrow(d.sites), function(ix) {
        if (ix%%1000 == 0) {
            cat("working on site ", ix, "\n")
        }
        s = d.sites[ix, ]$snp
        c = as.character((df.base %>% filter(snp == s))$chrom)
        min.freq = (df.base %>% filter(snp == s))$af.base.mean - 0.02
        max.freq = (df.base %>% filter(snp == s))$af.base.mean + 0.02
        d.site.shuff = df.base %>% filter(chrom == c) %>% filter(snp != 
            s & between(af.base.mean, min.freq, max.freq)) %>% sample_n(1)
        return(d.site.shuff)
    }, mc.cores = ncores))
}  

######################
#######################



####Shuffle Sites w/ in 5% starting frequency
get.MatchedSNP.ShuffleSites = function (x, base.data) 
{
    chr = x$chrom
    p = x$pos
    min.freq = x$af.base.mean - 0.05
    max.freq = x$af.base.mean + 0.05
    pos.matched = base.data %>% filter(chrom == chr & pos != 
        p & between(af.base.mean, min.freq, max.freq)) %>% sample_n(1) %>% 
        pull(pos)
    x = cbind(x, pos.matched)
    return(x)
}

##Shuffle sites w/in 1% starting frequency and get 10 matched sites
get.MatchedSNP.ShuffleSites.10x = function (x, base.data) 
{
    chr = x$chrom
    p = x$pos
    min.freq = x$af.base.mean - 0.01
    max.freq = x$af.base.mean + 0.01
    pos.matched = base.data %>% filter(chrom == chr & pos != 
        p & between(af.base.mean, min.freq, max.freq)) %>% sample_n(10) 
    matched.vec = pos.matched$pos
    d.matched = data.frame()
    d.matched = rbind(d.matched, matched.vec)
    names(d.matched) = c('pos.matched.1', 'pos.matched.2', 'pos.matched.3','pos.matched.4', 'pos.matched.5', 'pos.matched.6',
                    'pos.matched.7', 'pos.matched.8', 'pos.matched.9', 'pos.matched.10')
    d.new = cbind(x, d.matched)
    return(d.new)
}
########################################


####New get shift stats function specifically when running stats on each chromosomal arm
###This version only requires 5 observed snps to test
get.shift.stats.ByChrom  = function (df.shifts) 
{
    comps = as.character(unique(df.shifts$comparison))
    df.stat.meta = data.frame()
    for (comp in comps) {
        df.c = df.shifts %>% filter(comparison == comp)
        segs = grep(".target", names(df.c), value = TRUE)
        for (i in (1:length(segs))) {
            segs[i] = strsplit(segs[i], "[.]")[[1]][2]
        }
        for (seg in segs) {
            df.stat = data.frame()
            ID.Seg = comp
            Test.Seg = seg
            matched.vec = df.c %>% dplyr::select(paste0("dAF.", 
                seg, ".matched"))
            matched.vec = matched.vec[, 1]
            median.matched = median(matched.vec)
            target.vec = df.c %>% dplyr::select(paste0("dAF.", 
                seg, ".target"))
            target.vec = target.vec[, 1]
            median.target = median(target.vec)
            if (length(target.vec) > 5) {
                pvalue = t.test(target.vec, matched.vec)$p.value
            }
            else {
                pvalue = NA
            }
            df.stat = cbind(ID.Seg, Test.Seg, median.matched, 
                median.target, pvalue)
            df.stat.meta = rbind(df.stat.meta, df.stat)
        }
        df.stat.meta = rbind(df.stat.meta, df.stat)
    }
    return(df.stat.meta)
}

##############################

############Get Segment Shifts for Ecol Phases LOA - Orchard 2021
get.NextSegShifts.EcolPhasesLOA = function (df.sig.matched, RData, Comp, next.segs) 
{
    load(RData)
    df = df.sig.matched %>% filter(comparison == Comp) %>% mutate(sign.shift = sign(afShift), 
        phased.afShift = abs(afShift))
    af.shifts = get_af_shifts(afmat.lo, samps.lo, cage_set = NULL, 
        next.segs)
    af.shifts = cbind(sites, af.shifts)
    af.shifts.target = af.shifts
    colnames(af.shifts.target) = paste(colnames(af.shifts), "target", 
        sep = ".")
    colnames(af.shifts.target)[1:2] = c("chrom", "pos")
    af.shifts.matched = af.shifts
    colnames(af.shifts.matched) = paste(colnames(af.shifts), 
        "matched", sep = ".")
    colnames(af.shifts.matched)[1:2] = c("chrom", "pos")
    df.target = left_join((df %>% dplyr::select(chrom, pos)), 
        af.shifts.target) %>% dplyr::select(-chrom, -pos)
    signs = df$sign.shift
    df.target.phased = data.frame(nrow = nrow(df))
    for (cols in 1:ncol(df.target)) {
        col = df.target[, cols]
        col.new = signs * col
        df.target.phased = cbind(df.target.phased, col.new)
    }
    df.target.phased = as.data.frame(df.target.phased[, -1])
    colnames(df.target.phased) = colnames(df.target)
    df.matched = df %>% dplyr::select(chrom, pos.matched) %>% 
        rename(pos = pos.matched)
    df.matched = left_join(df.matched, af.shifts.matched) %>% 
        dplyr::select(-chrom, -pos)
    df.shifts = cbind(df, df.target.phased, df.matched)
    return(df.shifts)
}

############################
#########


####Function to get Sig Sites and Matched sites from orchard 2021 e cage data
Get.Sig.Match.orch21 = function(glmres, HAFsFile) {
    load(HAFsFile)
    load(glmres)
    seg.length = strsplit(glmres, "[_]")[[1]][2]
    sites = df.glm %>% dplyr::select(chrom, pos)
    comparisons = grep("coef.", names(df.glm), value = TRUE)
        for (i in (1:length(comparisons))) {
            comparisons[i] = strsplit(comparisons[i], "[.]")[[1]][2]
        }
    af.shifts <- get_af_shifts(afmat, samps, cage_set = NULL, 
                comparisons)
    FDR = get_glm_FDR.V2(df.glm)
    df.sig <- get_sig_sites(df.glm, comparisons, FDR, afShifts = af.shifts, 
        fdrThreshs, esThreshs)
    write.csv(df.sig, paste0("df.sig.ExpandingPhases.", seg.length, 
        ".orch21.csv"), row.names = FALSE)
    matched.sites = read.csv('~/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/ShuffledSites10X/Orch2021.ShuffledSites10x.csv')
    #Get one of the 10 random columns of shuffled sites
    rand.col = floor(runif(1, min=5, max=14))
    matched.sites.loc = matched.sites[,c(2,3,4,rand.col)]
    names(matched.sites.loc) = c('chrom', 'pos', 'af.base.mean', 'pos.matched')

    ##join df.sig to matched.sites dataframe
    df.sig.matched = left_join(df.sig, matched.sites.loc)

    write.csv(df.sig.matched, paste0("df.sig.Matched.ExpandingPhases.", seg.length, 
        ".orch21.csv"), row.names = FALSE)
    
    }
###########
########


###################New get next seg shift function that simply uses the loa RData file (accomodate samps.lo, afmat.lo (instead of samps, afmat))

get.NextSegShifts.V5 = function (df.sig.matched, RData, Comp, next.segs) 
{
    load(RData)
    df = df.sig.matched %>% filter(comparison == Comp) %>% mutate(sign.shift = sign(afShift), 
        phased.afShift = abs(afShift))
    af.shifts = get_af_shifts(afmat.lo, samps.lo, cage_set = NULL, 
        next.segs)
    af.shifts = cbind(sites, af.shifts)
    af.shifts.target = af.shifts
    colnames(af.shifts.target) = paste(colnames(af.shifts), "target", 
        sep = ".")
    colnames(af.shifts.target)[1:2] = c("chrom", "pos")
    af.shifts.matched = af.shifts
    colnames(af.shifts.matched) = paste(colnames(af.shifts), 
        "matched", sep = ".")
    colnames(af.shifts.matched)[1:2] = c("chrom", "pos")
    df.target = left_join((df %>% dplyr::select(chrom, pos)), 
        af.shifts.target) %>% dplyr::select(-chrom, -pos)
    signs = df$sign.shift
    df.target.phased = data.frame(nrow = nrow(df))
    for (cols in 1:ncol(df.target)) {
        col = df.target[, cols]
        col.new = signs * col
        df.target.phased = cbind(df.target.phased, col.new)
    }
    df.target.phased = as.data.frame(df.target.phased[, -1])
    colnames(df.target.phased) = colnames(df.target)
    df.matched = df %>% dplyr::select(chrom, pos.matched) %>% 
        rename(pos = pos.matched)
    df.matched = left_join(df.matched, af.shifts.matched) %>% 
        dplyr::select(-chrom, -pos)
    df.shifts = cbind(df, df.target.phased, df.matched)
    return(df.shifts)
}

##################################
################################

##Get intra segment phased allele frequency shifts at target and matched controls for all comparisons in df.sig.matched file
get.IntraSegShifts = function(df.sig.matched,HAFsFile, Comp, next.segs) {
    load(HAFsFile)
    df = df.sig.matched %>% filter(comparison == Comp) %>% mutate(sign.shift = sign(afShift), 
        phased.afShift = abs(afShift))
    af.shifts = get_af_shifts(afmat.lo, samps.lo, cage_set = NULL, 
            next.segs)
    af.shifts = cbind(sites, af.shifts)
    af.shifts.target = af.shifts
    load(HAFsFile)
    af.shifts.matched = get_af_shifts(afmat.lo, samps.lo, cage_set = NULL, 
       next.segs)

    af.shifts.matched = cbind(sites, af.shifts.matched)


    if (ncol(af.shifts) > 2) {
        colnames(af.shifts.target) = paste(colnames(af.shifts), 
            "target", sep = ".")
        colnames(af.shifts.target)[1:2] = c("chrom", "pos")
        af.shifts.matched = af.shifts
        colnames(af.shifts.matched) = paste(colnames(af.shifts), 
            "matched", sep = ".")
        colnames(af.shifts.matched)[1:2] = c("chrom", "pos")
        df.target = left_join((df %>% dplyr::select(chrom, pos)), 
            af.shifts.target) %>% dplyr::select(-chrom, -pos)
        signs = df$sign.shift
        df.target.phased = data.frame(nrow = nrow(df))
        for (cols in 1:ncol(df.target)) {
            col = df.target[, cols]
            col.new = signs * col
            df.target.phased = cbind(df.target.phased, col.new)
        }
        df.target.phased = as.data.frame(df.target.phased[, -1])
        colnames(df.target.phased) = colnames(df.target)
        df.matched = left_join((df %>% dplyr::select(chrom, pos.matched) %>% 
            rename(pos = pos.matched)), af.shifts.matched) %>% dplyr::select(-chrom, -pos)
        df.matched.phased = data.frame(nrow = nrow(df))
        for (cols in 1:ncol(df.matched)) {
            col = df.matched[, cols]
            col.new = 1 * col
            df.matched.phased = cbind(df.matched.phased, col.new)
        }
        df.matched.phased = as.data.frame(df.matched.phased[, -1])
        colnames(df.matched.phased) = colnames(df.matched)
         df.shifts = cbind(df, df.target.phased, df.matched.phased)
        names(df.shifts) = c('ix', 'chrom','pos', 'coef.div', 'p.div','sigLevel', 'FDR','afShift','comparison','snp','pos.matched','sign.shift',
			   'phased.afShift',"dAF.target",  "dAF.matched")
     }else {
        df.shifts = matrix(nrow = 1, ncol = 15)
        df.shifts = as.data.frame(df.shifts)
        names(df.shifts) = c('ix', 'chrom','pos', 'coef.div', 'p.div','sigLevel', 'FDR','afShift','comparison','snp','pos.matched','sign.shift',
			 'phased.afShift',"dAF.target",  "dAF.matched")
        df.shifts$comparison = next.segs
    }
     
    return(df.shifts)
}

###################
##################
#########


####Get stats from a df.shifts file generated with get.IntraSegShifts
get.shift.stats.IntraSeg = function (df.shifts) {
    comps = as.character(unique(df.shifts$comparison))
    df.stat.meta = data.frame()
    for (comp in comps) {
        df.c = df.shifts %>% filter(comparison == comp)
        df.stat = data.frame()
        ID.Seg = comp
        Test.Seg = comp
        matched.vec = df.c %>% dplyr::select(dAF.matched)
            matched.vec = matched.vec[, 1]
            median.matched = median(matched.vec)
            target.vec = df.c %>% dplyr::select(dAF.target)
            target.vec = target.vec[, 1]
            median.target = median(target.vec)
            if (length(target.vec) > 10) {
                pvalue = t.test(target.vec, matched.vec)$p.value
            }
            else {
                pvalue = NA
            }
            df.stat = cbind(ID.Seg, Test.Seg, median.matched, 
                median.target, pvalue)
            df.stat.meta = rbind(df.stat.meta, df.stat)
        }
        df.stat.meta = rbind(df.stat.meta, df.stat)
        return(df.stat.meta)

    }
##################


##############Get matched snps from base data, matching on recomb rate, inversion status, and distance from focal snp
get.MatchedSNP.BaseData = function(x, base.data) 
{
    chr = x$CHROM
    p = x$POS
    min.freq = x$af.base.mean - 0.025
    max.freq = x$af.base.mean + 0.025
    min.recomb = x$RECOM - 0.5
    max.recomb = x$RECOM + 0.5
    p.matched = as.data.frame(base.data %>% mutate(snp = paste0(CHROM, POS)) %>%
        filter(CHROM == chr & POS != p & between(af.base.mean, min.freq, max.freq)) %>%  #filter on chrom and starting freq.
        filter(Inv == x$Inv & between(RECOM, min.recomb, max.recomb)) %>% ##filter on inv status and recomb rate
        filter(POS > p + 50000  | POS < p - 50000)) %>% #filter on position and ensure 
        sample_n(15) %>% rename(pos.matched = POS) 
    p.matched = p.matched$snp
    matched = data.frame()
    matched = rbind(matched, p.matched)
    colnames(matched) = c('match1', 'match2', 'match3', 'match4', 'match5', 'match6', 'match7', 'match8', 'match9', 'match10',
                         'match11', 'match12', 'match13', 'match14', 'match15')
    d.new = cbind(x, matched)
    return(d.new)
    }

##################
#################


#####Select a single matched site from the matched site dataframe
select.matched = function (x, matched.sites, df.c) 
{
    p.target = x$snp
    matched = matched.sites %>% filter(snp == p.target)
    pos.matched = as.character(as.data.frame(t(matched[1, 8:22])) %>% 
        filter(!V1 %in% df.c$snp) %>% sample_n(1) %>% pull(1))
    d.new = cbind(x, pos.matched)
    return(d.new)
}

############
############