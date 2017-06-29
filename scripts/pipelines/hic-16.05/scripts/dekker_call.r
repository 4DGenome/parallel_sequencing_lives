#!/usr/bin/env Rscript

# dependencies

library("magrittr")
library("dplyr")
library("parallel")
library("hicutils")
library("Matrix")
library("matrixStats")
library("mgcv")
library("tidyr")
library("methods")
library("MASS")

options(stringsAsFactors = F, scipen = 999)

# function to get data

get_matrix <- function(infile, region, bins){

    paste("/software/mb/bin/tabix",
          infile,
          region) %>%
        pipe %>%
        read.delim(head = F) %>%
        setNames(c("chr", "p1", "p2", "n")) %>%
        mutate(bin1 = paste0(chr, ":", p1) %>% factor(levels = bins),
               bin2 = paste0(chr, ":", p2) %>% factor(levels = bins)) %>%
        xtabs(n ~ bin1 + bin2, ., sparse = T)

}


    
# function to storecompressed and indexed

write_tabix <- function(obj, filename, FUN = write.table, ...){
    
    zz <- paste("/software/mb/bin/bgzip -c >", filename) %>%
        pipe("w")

    on.exit(close(zz))
    
    FUN(obj, zz, ...)
    
}

# function to compute insulation index and get insulation TADs

get_insulation <- function(mat_file, p_is = 500000, p_ids = 250000, p_nt = 0.1, scratch = tmpdir){
    
    tmp <- tempfile(tmpdir = scratch)
    dir.create(tmp)
    oldir <- getwd()
    setwd(tmp)

    paste("perl -I /software/mb/el7.2/perl5/lib/perl5/",
          "/software/mb/el7.2/cworld-dekker/scripts/perl/matrix2insulation.pl",
          "--is ", p_is,
          "--ids ", p_ids,
          "--nt ", p_nt,
          "-i ", mat_file) %>%
        system

    borders <- list.files(patt = "insulation.boundaries.bed$") %>%
        read.delim(, skip = 1, head = F) %>%
        setNames(c("chr", "start", "end", "binname", "score")) %>%
        mutate(end = c(start[-1], NA))
    
    index <- list.files(patt = "insulation.bedGraph") %>%
        read.delim(skip = 1, head = F) %>%
        setNames(c("chr", "start", "end", "score"))
    
    list.files(tmp, patt = "is.*ids") %>% unlink
    unlink(tmp)

    setwd(oldir)

    list(borders = borders, index = index)

}

# funtion to output messages

catn <- function(...) cat(..., sep = "\n")


# settings

## reso <- 50e3
## infile <- paste0("/nfs/users/project/4DGenome_no_backup/data/hic/merged/ES_rep1/",
##                 "ES_rep1_raw_50kb.tsv.gz")
## ncores <- 10
## outfile <- "~/prueba.pdf"
## tempdir <- ""

# get arguments

input <- commandArgs(trailingOnly = TRUE)
infile <- input[1]
inbam <- input[2]
reso <- as.numeric(input[3])
ncores <- as.numeric(input[4])
outfile <- input[5]
p_is <- as.numeric(input[6])
p_ids <- as.numeric(input[7])
p_nt <- as.numeric(input[8])
tmpdir <- input[9]
 
if(tmpdir == "") tmpdir <- tempfile()

dir.create(tmpdir, F)

catn()
catn("Input parameters")
paste("Input file", infile) %>% catn
paste("Input bam", inbam) %>% catn
paste("Resolution", reso) %>% catn
paste("N cores", ncores) %>% catn
paste("Output suffix", outfile) %>% catn
paste("parameter 'IS'", p_is) %>% catn
paste("parameter 'IDS'", p_ids) %>% catn
paste("parameter 'NT'", p_nt) %>% catn
paste("Temp dir", tmpdir) %>% catn
catn()

# make bins

paste("Creating bins") %>% catn

bins <- make_bins(inbam, reso) %>%
    mutate(binid = paste(bin, pos + reso, sep = "-"))

chromosomes <- unique(bins$chr)

#chromosomes <- chromosomes[!(gsub("^chr", "", chromosomes) %>% as.numeric %>% is.na)]

# get contact matrix

mat_file <- tempfile(tmpdir = tmpdir)

paste("Getting matrices") %>% catn

mclapply(chromosomes, try(function(chrom){
    
    mat <- get_matrix(infile, chrom, filter(bins, chr == chrom) %$% bin)
    
    colnames(mat) <- rownames(mat) <- filter(bins, chr == chrom) %$% binid

    out <- mat %>%
        simetrize_matrix %>%
        band(-100, 100) %>%
        as.matrix

    ids <- paste0("bin",
                  1:nrow(out),
                  "|hg38|",
                  rownames(out))

    rownames(out) <- colnames(out) <- ids

    outd <- out
    aux <- outd
    diag(aux) <- 0
    diag(outd) <- max(aux)
    rm(aux)

    write.table(out, paste(mat_file, chrom, sep = "_"), quote = F, sep =  "\t", col.names = NA)
    
}, silent = T), mc.cores = ncores, mc.preschedule = F)


# for chromosome

paste("Getting insulation indexes") %>% catn

# compute insulation score and call TAD borders
    
aux <- mclapply(chromosomes, function(chrom){
    paste(mat_file, chrom, sep = "_") %>%
        (function(nn) try(get_insulation(nn,
                                         p_is = p_is,
                                         p_ids = p_ids,
                                         p_nt = p_nt),
                          silent = T))
}, mc.cores = ncores, mc.preschedule = F)

print(aux[sapply(aux, class) == "try-error"])

aux <- aux[sapply(aux, class) != "try-error"]

# arrange to store

out <- list()
out$borders <- lapply(aux, "[[", "borders") %>% do.call(rbind, .)
out$index <- lapply(aux, "[[", "index") %>% do.call(rbind, .)

# add last bin to last TAD border (per chromosome)

out$borders$end[is.na(out$borders$end)] <- with(group_by(bins, chr) %>%
                                                summarize(end = max(pos)),
                                                end[match(out$borders$chr[is.na(out$borders$end)], chr)])

names(out$borders)[1] <- names(out$index)[1] <- "#chr"

# store and index
    
write_tabix(out$borders,
            paste(outfile, "is", p_is, "ids", p_ids,
                  "nt", p_nt, "borders.bed.gz", sep ="_"),
            row.names = F, sep = "\t", quote = F)

paste("/software/mb/bin/tabix -f -p bed",
      paste(outfile, "is", p_is, "ids", p_ids,
            "nt", p_nt, "borders.bed.gz", sep ="_")) %>%
    system

write_tabix(out$index,
            paste(outfile, "is", p_is, "ids", p_ids,
                  "nt", p_nt, "index.bed.gz", sep ="_"),
            row.names = F, sep = "\t", quote = F)

paste("/software/mb/bin/tabix -f -p bed",
      paste(outfile, "is", p_is, "ids", p_ids,
            "nt", p_nt, "index.bed.gz", sep ="_")) %>%
    system

catn("Finished!")
