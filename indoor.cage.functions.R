###Functions for analysis of indoor cage data

####Fit continuous time glm
fit_GLM_ContinuousTime = function (afMatrix, rdMatrix, sampleData, vec, model.vars, poolCt = 100, 
    ncores) 
{
    registerDoMC(ncores)
    df = as.data.frame(sampleData[, colnames(sampleData) %in% 
        model.vars])
    colnames(df) <- model.vars
    formulaString = paste0(colnames(df), collapse = " + ")
    formulaString = paste0("cts ~ ", formulaString)
    cat("Model Formula is: \n", formulaString, "\n")
    Neff = calc_Neff(rdMatrix, poolCt)
    do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
        if (ix%%10000 == 0) {
            cat("working on site ", ix, "\n")
        }
        cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
            ] * (1 - afMatrix[ix, ]))))
        df$cts = cts
        df = as.data.frame(df)
        model = glm(formulaString, family = "quasibinomial", 
            data = df)
        cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
        results = c(cp[, 1], cp[, 2])
        names(results) = c(paste0("coef.", vec[1], "_", vec[length(vec)]), 
            paste0("p.", vec[1], "_", vec[length(vec)]))
        return(results)
    }, mc.cores = ncores))
}


####get N effective - written by S. Greenblum
calc_Neff = function (rd, poolCt) 
{
    ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}

########
########

###Get af shifts - written by S. Greenblum
get_af_shifts = function (afmat, samps, cage_set = NULL, comparisons) 
{
    df.shifts = do.call(cbind, lapply(comparisons, function(cc) {
        tt = rev(chop(cc, "_", 1:2))
        if (is.null(cage_set)) {
            cage_set = unique(samps$cage)
        }
        cageMask = samps$cage %in% cage_set
        t1Mask = cageMask & samps$tpt == tt[1]
        if (sum(t1Mask) > 1) {
            af1 = rowMeans(afmat[, t1Mask])
        }
        else {
            af1 = afmat[, t1Mask]
        }
        t2Mask = cageMask & samps$tpt == tt[2]
        if (sum(t2Mask) > 1) {
            af2 = rowMeans(afmat[, t2Mask])
        }
        else {
            af2 = afmat[, t2Mask]
        }
        ss = data.frame(dAF = af2 - af1) %>% rename_all(.funs = function(x) {
            paste0("dAF.", cc)
        })
        return(ss)
    }))
    return(df.shifts)
}
#######
#######


##get glm FDR - written by SG, modified by MCB
get_glm_FDR.V2 = function (df.glm) 
{
    pvals = df.glm %>% dplyr::select(contains("p."))
    pval.corr = data.frame(nrow = nrow(pvals))
    for (i in 1:ncol(pvals)) {
        pval.vec = p.adjust(pvals[, i], method = "BH")
        pval.corr = cbind(pval.corr, pval.vec)
    }
    FDR = as.data.frame(pval.corr[, -1])
    colnames(FDR) = colnames(pvals)
    FDR = as.matrix(FDR)
    return(FDR)
}

########
########


###Get significant sites from glm results generated with fit_GLM_ContinuousTime function

get.sig.sites = function(glm.file, rdata, comps, fdrThreshs, esThreshs){
                    chroms = c('2L', '2R', '3L', '3R', 'X')
                    load(rdata)
                    af.shifts = get_af_shifts(afmat, samps, cage_set = NULL, comparisons = comps)
                    #Get fdr 
                    load(glm.file)
                    FDR = get_glm_FDR.V2(df.glm = df.glm)
                    load(glm.file)
                    df.sig = get_sig_sites(df.glm, comparisons = comps, FDR  , af.shifts, fdrThreshs,
                                          esThreshs)
            return(df.sig)
                       }
#####

##another function to get singificant sites, but this is specifically geting the significant site dataframe
get_sig_sites = function (df.glm, comparisons, FDR, afShifts, fdrThreshs, esThreshs) 
{
    pSig = Reduce("+", lapply(1:length(fdrThreshs), function(ii) {
        (0 + (FDR <= fdrThreshs[ii] & abs(afShifts) >= esThreshs[ii]))
    }))
    do.call(rbind, lapply(comparisons, function(cc) {
        cc.ix = match(cc, comparisons)
        df.glm %>% mutate(ix = 1:nrow(df.glm)) %>% rename(coef.div = paste0("coef.", 
            cc), p.div = paste0("p.", cc)) %>% dplyr::select(ix, 
            chrom, pos, coef.div, p.div) %>% mutate(sigLevel = pSig[, 
            cc.ix], FDR = FDR[, cc.ix], afShift = afShifts[, 
            cc.ix], comparison = cc) %>% filter(sigLevel > 0, 
            !is.na(p.div))
    })) %>% mutate(comparison = factor(comparison, comparisons))
}

####Get matched snp - takes in df.sig generated with get_sig_sites and baseline allele frequency data for a particular cage

get.matchedSNP = function (x, snps.c, base.data)
{
    chr = as.character(x$chrom)
    p = x$pos
    min.freq = x$af.base - 0.025
    max.freq = x$af.base + 0.025
    min.recomb = x$RECOM - 0.5
    max.recomb = x$RECOM + 0.5
    i.s = x$Inv.Status
    p.matched = as.data.frame(base.data) %>% filter(chrom ==
        chr & pos != p & between(af.base.mean, min.freq, max.freq)) %>%
        filter(between(RECOM, min.recomb, max.recomb) & Inv.Status == i.s) %>% filter(pos >
        p + 10000 | pos < p - 10000) %>% filter(!snp %in% snps.c) %>%
        sample_n(20)
    p.matched = as.data.frame(t((p.matched)$pos))
    d.new = cbind(x, p.matched)
    return(d.new)
}




###Get significant sites using observed distributions of FDR and afShift for each chromosome - for Clustering
get.sig.ClusterSites = function(glm.file, rdata, comps){
                    chroms = c('2L', '2R', '3L', '3R', 'X')
                    load(rdata)
                    af.shifts = get_af_shifts(afmat, samps, cage_set = NULL, comparisons = comps)
                    #Get fdr 
                    load(glm.file)
                    FDR = get_glm_FDR.V2(df.glm = df.glm)

                    #Determine quantiles of FDR and afShift for sigLevel filtering by chromosomal arm
                    af = as.data.frame(cbind(sites, af.shifts))
                    fdr = as.data.frame(cbind(sites, FDR))

                    df.scores = data.frame()
                    for (c in chroms){
                        af.c = af %>% filter(chrom == c)
                        fdr.c = fdr %>% filter(chrom == c)
                        afs = sort(as.vector(quantile(abs(af.c[,3]), probs = c(0.85, 0.75, 0.5))))
                        ps = sort(as.vector(quantile(fdr.c[,3], probs = c(0.2, 0.1, 0.05))), decreasing = TRUE)
                        df = as.data.frame(cbind(afs, ps))
                        df$score = c(1, 2, 3)
                        df$chrom = c
                        df.scores = rbind(df.scores, df)
                    }

                    thresh = df.scores
                    names(thresh) = c('af.shift', 'FDR', 'score', 'chrom')

                    #Generate sig files separately for each chromosome w/ specific sigLevel parameters, then bind to a meta df.sig
                    df.sig.meta = data.frame()
                    for (chr in chroms){
                        load(glm.file)
                        fdr.c = as.data.frame(cbind(sites, FDR) %>% filter(chrom == chr) %>% dplyr::select(-chrom, -pos))
                        df.glm.c = df.glm %>% filter(chrom == chr)
                        af.shifts.c = as.data.frame(cbind(sites, af.shifts) %>% filter(chrom == chr) %>% dplyr::select(-chrom, -pos))
                        thresh.c = thresh %>% filter(chrom == chr)
                        fdrThreshs=sort(as.numeric(thresh.c$FDR),decreasing = TRUE) 
                        esThreshs=sort(as.numeric(thresh.c$af.shift))
                        df.sig = get_sig_sites(df.glm = df.glm.c, comparisons = comps, FDR = fdr.c, afShifts = af.shifts.c, fdrThreshs,
                                          esThreshs)
                        df.sig.meta = rbind(df.sig.meta, df.sig)

                        }
                    return(df.sig.meta)
                    }

################
###############

###function to project fst-mds points onto linear axis derived from other time points/phases
project_points=function(m1,b1,myx,myy){
  m2=-1*1/m1
  b2=myy-m2*myx
  project.x=(b2-b1)/(m1-m2)
  project.y=b1+m1*project.x
  mydist=sqrt(project.x^2 + (project.y-b1)^2)*c(-1,1)[1+as.numeric(project.x>0)];
  return(cbind(project.x,project.y,mydist))
}


#########
#########

##function to get p-value from hypergeometric test (used to compare overlap in SNPs identified in indoor and outdoor mesocosm)
hyperg.p = function(test.snps, comp.snps, num.sites){
    overlap = length(intersect(test.snps, comp.snps))
    test.length = length(test.snps)
    comp.length = length(comp.snps)
    h.p = (phyper(overlap - 1, test.length, num.sites - test.length, comp.length, lower.tail = FALSE, log.p = FALSE))
    return(h.p)
}