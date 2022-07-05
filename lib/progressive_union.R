.libPaths(c(.libPaths(),"/home/jbonnie1/rlibs/4.0.2/","/home/jbonnie1/R/4.0.4","/usr/local/lib/R/site-library","/home/jbonnie1/rlibs/4.0.2/gcc/9.3.0","/usr/local/lib/R/library"))
require(tidyr)
require(ggplot2)
require(data.table)
require(openssl)
require(dplyr)
require(hash)
# require(digest)
source('/home/jbonnie1/scr16_blangme2/jessica/dandd/dev-dandD/lib/plot_progressive.R')


library("optparse")

option_list = list(
  make_option(c("-x", "--species"), type="character", default=NULL, 
              help="species tagname used to locate data directory", metavar="character"),
              # option_list=c('HVSVC2','human','salmonella','ecoli')),
  
  make_option(c("-o", "--out"), type="character", default=".", 
              help="output directory path [default= %default]", metavar="character"),
  
  make_option(c("-K", "--maxk"), type="numeric", default=20, 
              help="maximum value of k to sketch [default= %default]", metavar="numeric"),
  
  make_option(c("-k", "--mink"), type="numeric", default=14, 
              help="minimum value of k to sketch [default= %default]", metavar="numeric"),
  
  make_option(c("-g", "--ngenomes"), type="numeric", default=NULL, 
              help="number of genomes to include in the random orderings [default= %default]", metavar="numeric"),
  
  make_option(c("-n", "--norder"), type="numeric", default=30, 
              help="Number of random orderings to sketch [default= %default]", metavar="numeric"),
  
  make_option(c("-d", "--dashing"), type="character", default="/home/jessica/lib/dashing/dashing", 
              help="location of dashing program [default= %default]", metavar="character"),
  
  make_option(c("-a", "--append"), type="logical", default=FALSE, action="store_true",
              help="append new random orderings to previously orderings [default= %default]", metavar="logical"),
  
  make_option(c("-s", "--sortorder"), type="character", default=NULL, 
              help="location of file with 2 columns (no titles): 1) name of fasta file and 2) expected position in the sorted order. If not provided it will be sought in outdir or created there."),
  
  make_option(c("-i", "--inputdir"), type="character", default='/home/jbonnie1/scr16_blangme2/jessica/data', 
              help="absolute path to directory of genomic data with species names as subdirectories containing fasta files."),
  
  make_option(c("-r", "--registers"), type="numeric", default=20, 
              help="number of registers to use during sketching[default= %default]", metavar="numeric")
  
); 
# print(getwd())
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

norder <- opt$norder
tag <- opt$species
dashingcmd <- opt$dashing
mink <- opt$mink
maxk <- opt$maxk
outdir <- R.utils::getAbsolutePath(opt$out)
print(outdir)


sketchdir<-file.path(outdir,"sketches",tag)
dir.create(sketchdir, showWarnings = FALSE)
gdir_base<-opt$inputdir
nregister=opt$registers

hashkey_loc=file.path(sketchdir,paste0(tag,"_hashkey.rda"))
oldresultsfile=file.path(outdir,paste0(tag,"_results.rda"))
resultsfile=file.path(sketchdir, paste0(tag,"_results.rda"))

if(file.exists(oldresultsfile)){ file.rename(oldresultsfile, resultsfile)}

parval=20

# dashingcmd<- '~/lib/dashing/dashing'
# tag <- 'salmonella'
species=tag
genomedir <- file.path(gdir_base,tag)
if (tag == 'HVSVC2'){ 
  genomedir <- file.path(gdir_base,'HVSVC2','consensus')
  species='HVSVC2'}

if (tag == 'HVSVC2_2'){ 
  genomedir <- file.path(gdir_base,'HVSVC2','consensus_snv_sv')
  species='HVSVC2'}
if (tag == 'HVSVC2_snv'){ 
  genomedir <- file.path(gdir_base,'HVSVC2','consensus_snv')
  species='HVSVC2'}
if (tag == 'HVSVC2_indel'){ 
  genomedir <- file.path(gdir_base,'HVSVC2','consensus_indel')
  species='HVSVC2'}
if (tag == 'HVSVC2_sv'){ 
  genomedir <- file.path(gdir_base,'HVSVC2','consensus_sv')
  species='HVSVC2'}
if (tag == 'HVSVC2_snv_indel'){ 
  genomedir <- file.path(gdir_base,'HVSVC2','consensus_snv_indel')
  species='HVSVC2'}


# setwd(opt$out)
# setwd("~/scr16_blangme2/jessica/dandd/progressive_union/")
# results.save = list()
# results.save<-list(salmonella=list(),human=list(),ecoli=list())
if (file.exists(resultsfile)){
load(resultsfile)} else{ results.save=list()}
if (file.exists(hashkey_loc)){
  load(hashkey_loc)} else{ hashkey=hash()}
# load('orderings.rda')

fasta_files <- sort(list.files(genomedir,pattern = "*\\.(fa.gz|fasta.gz|fna.gz|fasta|fa)"))

old_expected_sort<-file.path(outdir,paste0(tag,"_sortorder.csv"))
expected_sort<-file.path(sketchdir,paste0(tag,"_sortorder.csv"))
if (file.exists(old_expected_sort)){file.rename(old_expected_sort, expected_sort)}

# print(expected_sort)
if (is.null(opt$sortorder)){
  if (! file.exists(expected_sort)){
    gfiles <- fasta_files
    # write.csv(file=expected_sort, data.frame(gfiles,seq(1,length(gfiles))),col.names = FALSE, row.names=FALSE, quote=FALSE)
  } else {
    order.df <- read.csv(expected_sort, header = TRUE, strip.white = TRUE,stringsAsFactors = FALSE) %>%
      arrange(position)
    complement <- fasta_files[!fasta_files %in% order.df$name]
    if (length(complement) > 0){
    gfiles <- append(order.df$names,complement)
    } else {
      gfiles <- order.df$name
      # print(gfiles)
      }
    
  }
  order.df <- data.frame(name=gfiles,position=seq(1,length(gfiles)))
  write.csv(file=expected_sort, x=order.df, row.names=FALSE, quote=FALSE)
} else {
  provided_order <- read.csv(expected_sort, header = TRUE, strip.white = TRUE,stringsAsFactors = FALSE) %>%
    arrange(position)
  if (file.exists(expected_sort)){
    order.df <- read.csv(expected_sort, header = TRUE, strip.white = TRUE,stringsAsFactors = FALSE) %>%
      arrange(position)
    sanity <- inner_join(provided_order,order.df,by='position')
    if (any(sanity$name.x != sanity$name.y)){ stop("OH NO, THE PROVIDED ORDER DOES NOT MATCH WHAT WE HAVE BEEN USING!! 
                                                             You need to bomb all of you ngen folders greater than 1. Or else, run a function that
                                                             isn't here yet to fix them.")
      error(1)}
    gfiles <- provided_order$name
    
  }
}

get_hash_string <-function(filenames){
  outvect=character(length=length(filenames))
  for (i in seq_along(filenames)){
    outvect[i] <- hashkey[[filenames[i]]]
  }
  return(outvect)}

assign_hash_string <- function(filename, len=6){
  print(paste0("%%%",filename,"%%%"))
  # newalphanum=digest::digest(file=filename, algo='md5', serialize = FALSE)
  # newtrunc=stringr::str_trunc(alphanum, len, "right",ellipsis = '')
  
  alphanum=digest::digest(filename, algo='md5', serialize = FALSE)
  trunc=stringr::str_trunc(alphanum, len, "right",ellipsis = '')
  print(paste0(filename, ' file keyed to hash ', trunc))
  if (has.key(key = filename, hashkey)){
    return(hashkey[[filename]])
  }
  if (trunc %in% values(hashkey)){
    warning(paste0('The hashvalue ', trunc, ' has 2 keys! ', filename, ' will be assigned to a longerhash.'))
    assign_hash_string(filename, len=len+1)
  }
  else { hashkey[[filename]] = trunc}
}

# print(gfiles)
gcount<-ifelse(is.null(opt$ngenomes), length(gfiles), opt$ngenomes)
print(paste0("I KNOW GCOUNT AND IT IS: ", gcount))
# print(paste("gcount:",gcount))
gfiles<-gfiles[1:gcount]
lapply(gfiles, function(x){
  assign_hash_string(x)
})
save(hashkey, file=hashkey_loc)

orderings_loc<-file.path(outdir,paste0(species,"_",gcount,"_orderings.rda"))
new_orderings <- list()
if(opt$append){
new_orderings <- lapply(1:norder,function(c){sample(1:gcount,replace = FALSE)})
}
if(file.exists(orderings_loc)){
  load(orderings_loc)
  orderings <- unique(c(orderings,new_orderings))
}else {
 orderings <- new_orderings 
}

save(orderings, file=orderings_loc)

nameSketch2<- function(indices, gfiles, kval, registers){
  # Function to return a hash independent of the ordering of the input files
  if (length(indices) == 1){
    return(paste0(basename(gfiles[indices[[1]]]),".w.",kval,".spacing.",registers,".hll"))
    
  } else {
    filename <-paste0(tag,"_",paste0(sort(indices),collapse="_"),"k",kval,"_r",registers,".hll")
    # string<- as.character(paste0(paste0(sort(input_files),collapse=""),"k",kval))
    
    # filename=paste0(digest::digest(string, algo='md5', serialize = FALSE),kval,length(input_files),".hll")
    return(filename)
  }
}

nameSketch<-function(input_files,kval,registers){
  # Function to return a hash independent of the ordering of the input files
  if (length(input_files) == 1){
    return(paste0(basename(input_files[[1]]),".w.",kval,".spacing.",registers,".hll"))
  } else {
    hash_keys<-get_hash_string(input_files)
    sorted_keys <- sort(hash_keys)
    filename<- as.character(paste0(tag,"_",paste0(sorted_keys,collapse="_"),"_k",kval,"_r",registers,".hll"))
    return(filename)
  }
}

nameSketch3<-function(input_files,kval, gfiles, registers){
  # Function to return a hash independent of the ordering of the input files
  if (length(input_files) == 1){
    return(paste0(basename(input_files[[1]]),".w.",kval,".spacing.",registers,".hll"))
  } else {
    recapitulate=sort(which(gfiles %in% input_files))
    # oldname=nameSketch(input_files,kval)
    filename<- as.character(paste0(tag,"_",paste0(recapitulate,collapse="_"),"_k",kval,"_r",registers,".hll"))
    
    return(filename)
  }
}

reconstituteGenomeList <- function(filename){
  parts <- unlist(stringr::str_split(filename, "_"))
  return(array(parts[2:(length(parts)-1)]))
}

# 
# nameSketch<-function(input_files,kval){
#   # Function to return a hash independent of the ordering of the input files
#   if (length(input_files) == 1){
#     return(paste0(basename(input_files[[1]]),".w.",kval,".spacing.10.hll"))
#   }
#   else {
#     string<- as.character(paste0(paste0(sort(input_files),collapse=""),"k",kval))
#     
#     filename=paste0(digest::digest(string, algo='md5', serialize = FALSE),kval,length(input_files),".hll")
#     return(filename)
#   }
# }



omni_list<-list()

for (i in 1:length(orderings)){
  reorder <- gfiles[orderings[[i]]]
  reorder_paths <- unlist(lapply(reorder,function(x){file.path(genomedir,x)}))
  
  for (kval in mink:maxk){
    sketchkdir=file.path(sketchdir, paste0("k",kval))
    dir.create(sketchkdir, showWarnings = FALSE)
    if (i ==1){
      sketchkndir=file.path(sketchkdir,paste0("ngen",1))
      dir.create(sketchkndir, showWarnings = FALSE)
      for (fname in gfiles){
        newsketchprefix <- nameSketch(c(fname),kval, registers=nregister)
        sketchprefix <- nameSketch3(c(fname),kval,gfiles,registers=nregister)
        sketch_loc <- file.path(sketchkndir, sketchprefix)
        newsketch_loc <- file.path(sketchkndir, newsketchprefix)
        if (! file.exists(sketch_loc) | file.size(sketch_loc) == 0L){
          # print(file.path(sketchdir,sketchprefix))
          command=paste0("~/lib/dashing/dashing sketch -k", kval," -p", parval,
                         " --prefix ", sketchkndir, " -S ",nregister," " ,file.path(genomedir, fname))
          print(command)
          system(command, ignore.stdout = FALSE, ignore.stderr = FALSE)
       }
      } }
    # reorder_sketches <- unlist(lapply(reorder, function(x){
    #   # oldpath=file.path(sketchdir,nameSketch(x,kval))
    #   # newpath=file.path(sketchkndir,nameSketch3(x,kval,gfiles, registers=nregister))
    #   # if (file.exists(oldpath)){
    #   #   print(paste("Old Name:",oldpath))
    #   #   print(paste("New Name:", newpath))
    #   #   file.rename(oldpath, newpath)
    #   # }
    #   return(newpath)
    # }))
    # print(reorder_sketches)
    for (g in 1:gcount){
      sketchkndir <- file.path(sketchkdir,paste0("ngen",g))
      dir.create(sketchkndir, showWarnings = FALSE)
      oldunionprefix <- file.path(sketchkndir,nameSketch3(reorder[1:g], kval,gfiles, registers=nregister))
      unionprefix <- file.path(sketchkndir,nameSketch(reorder[1:g], kval,registers=nregister))
      if (file.exists(oldunionprefix)){
        print(paste("#### I rename",oldunionprefix, "to", unionprefix))
        file.rename(oldunionprefix, unionprefix)
      }
      
      
      omni_list=append(omni_list,list(list(ordering=i, ngenomes=g, kval=kval, union_loc=unionprefix)))
      # print(paste("ngenomes",g," kval ",kval))
      if (g>1 ){
        # print(paste( reorder_paths[1:g], collapse=" "))
        old_alt_input1= file.path(sketchkdir,paste0("ngen",g-1),nameSketch3(reorder[1:(g-1)], kval, gfiles,registers=nregister))
        old_input2= file.path(sketchkdir,paste0("ngen",1),nameSketch3(reorder[g], kval, gfiles,registers=nregister))
        
        alt_input1= file.path(sketchkdir,paste0("ngen",g-1),nameSketch(reorder[1:(g-1)], kval, registers=nregister))
        alt_input2= file.path(sketchkdir,paste0("ngen",1),nameSketch(reorder[g], kval, registers=nregister))
        
        # command <- paste0("~/lib/dashing/dashing union -p ", parval," -z -o ", unionprefix, " ", paste( reorder_sketches[1:g], collapse=" "))
        command <- paste0("~/lib/dashing/dashing union -p ", parval," -z -o ", unionprefix, " ", alt_input1, " ", alt_input2)
        
        print(command)
        if (! file.exists(unionprefix) | file.size(unionprefix) == 0L){
          system(command, ignore.stdout = FALSE)
        }
      }
    }
  }
}



omni.table<-rbindlist(omni_list)
old_index_loc=paste0('/scratch16/blangme2/jessica/dandd/progressive_union/',tag,'_',length(orderings),'_index.tsv')
index_loc=file.path(sketchdir,paste0(tag,'_',length(orderings),'_index.tsv'))
if(file.exists(old_index_loc)){file.rename(old_index_loc,index_loc)}

old_card_loc=paste0("/scratch16/blangme2/jessica/dandd/progressive_union/",tag,"_",length(orderings),"_cardinalities.txt")
card_loc=file.path(sketchdir,paste0(tag,"_",length(orderings),"_cardinalities.txt"))
if(file.exists(old_card_loc)){file.rename(old_card_loc,card_loc)}



write.table(unique(omni.table$union_loc),index_loc, row.names = FALSE, col.names = FALSE, quote = FALSE)

command <- paste0("~/lib/dashing/dashing card --presketched -F ", index_loc," -o ",card_loc)
print(command)
system(command,wait = TRUE)
# cmd_output <-fread(cmd = paste0("~/lib/dashing/dashing card --presketched -F ", index_loc))

card <- read.table(card_loc, header = F, stringsAsFactors = F, skip=1) %>% 
  rename(union_loc=V1, UniqXMatch=V2)
# cbind(omni.table,card) %>% mutate(union_loc != union_loc2) %>% View()
delta.table <- left_join(omni.table,card) %>%
  mutate(ngenomes=as.numeric(ngenomes),delta_pos=UniqXMatch/kval) %>%
  group_by(ordering,ngenomes) %>%
  slice_max(delta_pos) %>% 
  ungroup() %>%  
  filter(UniqXMatch > 0) 


# if (opt$append){
#   results.save[[length(results.save)]] <- list(
#     card = card,
#     omni = omni.table,
#     delta = delta.table,
#     # tp = tp,
#     orderings = orderings
#   )
# } else {

results.save[[paste0("ng",gcount)]] = list(
                               card = card,
                               omni = omni.table,
                               delta = delta.table,
                               # tp = tp,
                               orderings = orderings,
                               gfiles=gfiles,
                               sort=order.df
                             )
# }
save(results.save,file=resultsfile)
print("SAVED RESULTS FILE!")

tp <- plotProgressiveUnion(species=tag, out = sketchdir, ngenome=gcount, delta.table)

# tp <-
#   delta.table %>%
#   mutate(Indicator=ordering <= 6) %>%
#   # filter(ordering %in% c(1:20)) %>%
#   mutate(ngenomes=as.integer(ngenomes),kval=as.factor(kval))  %>%
#   ggplot(.) +
#   geom_smooth(aes(y=delta_pos, x=ngenomes, linetype="Loess"), method='loess',formula=y~x) +
#   geom_point(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes, shape=kval), size=2) +
#   scale_linetype_manual(name="Fitted Curve",values=c(2)) +
#   geom_line(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes, group=ordering,color=as.factor(ordering)), size=.5) +
#   theme_bw() +
#   scale_x_continuous(breaks= scales::pretty_breaks(10)) +
#   scale_shape_discrete(name="argmax(k)")+
#   # scale_color_discrete(name = "Random Genome Ordering") +
#   xlab("Number of Genomes in Set") +
#   ylab("Value of \u03b4") +
#   ggtitle(label =paste0("Values of \u03b4 over orderings of ",gcount," ",stringr::str_to_title(tag)," Genomes"),
#           subtitle=paste0("Fit over ",length(orderings)," Orderings")) +
#   guides(color = 'none')
# # guides(color=guide_legend(ncol=3)) +
# # theme(legend.position = "bottom")
# 
# 
# 
# ggplot2::ggsave(filename = file.path(outdir,paste0("plot_",tag,"_ng",gcount,"_",length(orderings),".pdf")), 
#                 plot = tp, 
#                 device = cairo_pdf, 
#                 dpi = 1200, 
#                 width = 15,
#                 height = 10, 
#                 units = "cm")
