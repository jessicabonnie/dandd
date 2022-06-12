.libPaths(c("/home/jbonnie1/rlibs/4.0.2/","/home/jbonnie1/R/4.0.4","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()))
require(tidyr)
require(ggplot2)
require(data.table)
require(openssl)
require(dplyr)


library("optparse")

option_list = list(
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="species tagname used to locate data directory", metavar="character"),
              # option_list=c('HVSVC2','human','salmonella','ecoli')),
  
  make_option(c("-o", "--out"), type="character", default=".", 
              help="output directory path [default= %default]", metavar="character"),
  
  make_option(c("-K", "--maxk"), type="numeric", default=20, 
              help="maximum value of k to sketch [default= %default]", metavar="numeric"),
  
  make_option(c("-k", "--mink"), type="numeric", default=14, 
              help="minimum value of k to sketch [default= %default]", metavar="numeric"),
  
  make_option(c("-g", "--ngenomes"), type="numeric", default=10, 
              help="number of genomes to include in the random orderings [default= %default]", metavar="numeric"),
  
  make_option(c("-n", "--norder"), type="numeric", default=30, 
              help="Number of random orderings to sketch [default= %default]", metavar="numeric"),
  
  make_option(c("-d", "--dashing"), type="character", default="/home/jessica/lib/dashing/dashing", 
              help="location of dashing program [default= %default]", metavar="character")
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
resultsfile=file.path(outdir,paste0(tag,"_results.rda"))
sketchdir<-file.path(outdir,"sketches",tag)
dir.create(sketchdir, showWarnings = FALSE)


parval=20
gdir_base='/home/jbonnie1/scr16_blangme2/jessica/data'

# dashingcmd<- '~/lib/dashing/dashing'
# tag <- 'salmonella'

genomedir <- file.path(gdir_base,tag)
if (tag == 'HVSVC2'){ 
  genomedir <- file.path(gdir_base,tag,'consensus')}


# setwd(opt$out)
# setwd("~/scr16_blangme2/jessica/dandd/progressive_union/")
results.save = list()
# results.save<-list(salmonella=list(),human=list(),ecoli=list())
# load(resultsfile)
# load('orderings.rda')



gfiles<- list.files(genomedir,pattern = "*.fa.gz")
print(gfiles)
gcount<-unlist(ifelse(is.null(opt$ngenomes), length(gfiles), opt$ngenomes))
# print(paste("gcount:",gcount))
gfiles<-gfiles[1:gcount]
orderings <- lapply(1:norder,function(c){sample(1:gcount,replace = FALSE)})



nameSketch2<- function(indices, gfiles, kval){
  # Function to return a hash independent of the ordering of the input files
  if (length(indices) == 1){
    return(paste0(basename(gfiles[indices[[1]]]),".w.",kval,".spacing.10.hll"))
  }
  else {
    filename <-paste0(tag,"_",paste0(sort(indices),collapse="_"),"k",kval,".hll")
    # string<- as.character(paste0(paste0(sort(input_files),collapse=""),"k",kval))
    
    # filename=paste0(digest::digest(string, algo='md5', serialize = FALSE),kval,length(input_files),".hll")
    return(filename)
  }
}

nameSketch3<-function(input_files,kval, gfiles){
  # Function to return a hash independent of the ordering of the input files
  if (length(input_files) == 1){
    return(paste0(basename(input_files[[1]]),".w.",kval,".spacing.10.hll"))
  }
  else {
    recapitulate=sort(which(gfiles %in% input_files))
    # oldname=nameSketch(input_files,kval)
    filename<- as.character(paste0(tag,"_",paste0(recapitulate,collapse="_"),"_k",kval,".hll"))
    
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
        # sketchprefix <- nameSketch(c(fname),kval)
        sketchprefix <- nameSketch3(c(fname),kval,gfiles)
        sketch_loc <- file.path(sketchkndir, sketchprefix)
        if (! file.exists(sketch_loc) | file.size(sketch_loc) == 0L){
          # print(file.path(sketchdir,sketchprefix))
          command=paste0("~/lib/dashing/dashing sketch -k", kval," -p", parval,
                         " --prefix ", sketchkndir, " ", file.path(genomedir, fname))
          print(" ")
          print(command)
          print(" ")
          system(command, ignore.stdout = FALSE, ignore.stderr = FALSE)
       }
      } }
    reorder_sketches <- unlist(lapply(reorder, function(x){
      # oldpath=file.path(sketchdir,nameSketch(x,kval))
      newpath=file.path(sketchkndir,nameSketch3(x,kval,gfiles))
      # if (file.exists(oldpath)){
      #   print(paste("Old Name:",oldpath))
      #   print(paste("New Name:", newpath))
      #   file.rename(oldpath, newpath)
      # }
      return(newpath)
    }))
    print(reorder_sketches)
    for (g in 1:gcount){
      sketchkndir <- file.path(sketchkdir,paste0("ngen",g))
      dir.create(sketchkndir, showWarnings = FALSE)
      unionprefix <- file.path(sketchkndir,nameSketch3(reorder[1:g], kval,gfiles))
      omni_list=append(omni_list,list(list(ordering=i, ngenomes=g, kval=kval, union_loc=unionprefix)))
      # print(paste("ngenomes",g," kval ",kval))
      if (g>1 ){
        # print(paste( reorder_paths[1:g], collapse=" "))
        alt_input1= file.path(sketchkdir,paste0("ngen",g-1),nameSketch3(reorder[1:(g-1)], kval, gfiles))
        alt_input2= file.path(sketchkdir,paste0("ngen",1),nameSketch3(reorder[g], kval, gfiles))
        
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
index_loc=paste0('/scratch16/blangme2/jessica/dandd/progressive_union/',tag,'_index.tsv')
card_loc=paste0("/scratch16/blangme2/jessica/dandd/progressive_union/",tag,"_cardinalities.txt")
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




tp <- 
  delta.table %>% 
  mutate(Indicator=ordering <= 6) %>%
  # filter(ordering %in% c(1:20)) %>%
  mutate(ngenomes=as.integer(ngenomes),kval=as.factor(kval))  %>%
  ggplot(.) + 
  geom_smooth(aes(y=delta_pos, x=ngenomes, linetype="Loess"), method='loess',formula=y~x) +
  geom_point(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes, shape=kval), size=2) + 
  scale_linetype_manual(name="Fitted Curve",values=c(2)) + 
  geom_line(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes, group=ordering,color=as.factor(ordering)), size=.5) + 
  theme_bw() + 
  scale_x_continuous(breaks= scales::pretty_breaks(10)) +
  scale_shape_discrete(name="argmax(k)")+
  # scale_color_discrete(name = "Random Genome Ordering") + 
  xlab("Number of Genomes in Set") + 
  ylab("Value of \u03b4") + 
  ggtitle(label =paste0("Values of \u03b4 over orderings of ",gcount," ",stringr::str_to_title(tag)," Genomes")) +
  guides(color = FALSE)
# guides(color=guide_legend(ncol=3)) + 
# theme(legend.position = "bottom")




results.save = append(results.save,
                             list(list(
                               card = card,
                               omni = omni.table,
                               delta = delta.table,
                               tp = tp,
                               orderings = orderings
                             )))

save(results.save,file=resultsfile)

ggplot2::ggsave(filename = paste0("plot_",tag,length(orderings),".pdf"), 
                plot = tp, 
                device = cairo_pdf, 
                dpi = 1200, 
                width = 15,
                height = 10, 
                units = "cm")
