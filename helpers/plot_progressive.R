.libPaths(c("~/rlibs/4.0.2/gcc/9.3.0/","~/R/rstudio/4.0",.libPaths()))
#.libPaths(c(.libPaths(),"~/rlibs/4.0.2/","~/rlibs/4.0.2/gcc/9.3.0","~/R/4.0.4","/usr/local/lib/R/site-library","/usr/local/lib/R/library"))

require(tidyr)
require(ggplot2)
require(data.table)
require(openssl)
require(dplyr)


# tag="HVSVC2_2"
# outdir="/home/jbonnie1/scr16_blangme2/jessica/dandd/progressive_union"
# gcount=10
# load(resultsfile)
#library("optparse")

# option_list = list(
#   make_option(c("-s", "--species"), type="character", default=NULL, 
#               help="species tagname used to locate data directory", metavar="character"),
#   # option_list=c('HVSVC2','human','salmonella','ecoli')),
#   
#   make_option(c("-o", "--out"), type="character", default=".", 
#               help="output directory path [default= %default]", metavar="character"),
#   
#   make_option(c("-g", "--ngenomes"), type="numeric", default=10, 
#               help="number of genomes to include in the random orderings [default= %default]", metavar="numeric")
# ); 
# 
# resultsfile=file.path(outdir,paste0(tag,"_results.rda"))
# load(resultsfile)
# ngtag=paste0("ng",opt$ngenomes)
# 
# plot_progressive(delta=results.save[[ngtag]]$delta, out=opt$out, ngenomes=opt$ngenomes, species=opt$species)
# 
# 
# 
# opt_parser = OptionParser(option_list=option_list, resultsfile);
# opt = parse_args(opt_parser);

plotProgressiveUnion <- function(species, out, delta, nshow=Inf, method='loess'){
tag=species
outdir=out
gcount=max(delta$ngenomes)
item=paste0("ng",gcount)
norder=max(delta$ordering)

tp <-
  delta %>%
  mutate(Indicator=ordering <= nshow) %>%
  # filter(ordering %in% c(1:20)) %>%
  mutate(ngenomes=as.integer(ngenomes),kval=as.factor(kval))  %>%
  ggplot(.) +
  geom_smooth(aes(y=delta_pos, x=ngenomes, linetype=method), method=method,formula=y~x) +
  geom_point(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes, shape=kval), size=.5) +
  scale_linetype_manual(name=paste0("Fit (",norder," Orderings)"),values=c(2)) +
  geom_line(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes,
                                                     group=ordering, color=as.factor(ordering)), size=.5) +
  theme_bw() +
  scale_x_continuous(breaks= scales::pretty_breaks(10)) +
  scale_shape_discrete(name="argmax(k)")+
  # scale_color_discrete(name = "Random Genome Ordering") +
  xlab("Number of Genomes in Set") +
  ylab("Value of \u03b4") +
  ggtitle(label=paste0("Values of \u03b4* over orderings of ",gcount," ",stringr::str_to_title(tag)," Genomes")) +
  guides(color = 'none') +
  scale_y_continuous(labels = scales::label_number_auto())


ggplot2::ggsave(filename = file.path(outdir,paste0("plot_",tag,"_",item,"_",norder,".pdf")), 
                plot = tp, 
                device = cairo_pdf, 
                dpi = 1200, 
                width = 15,
                height = 10, 
                units = "cm")
return(tp)
}


plotAbba <- function(progu, step=3){
  # ordkey <- progu %>%
  #   select(ordering, step, flist4) %>%
  #   filter(step %in% c(3,4,5)) %>%
  #   pivot_wider(names_from=step,values_from=flist4, names_prefix = "step") %>%
  #   inner_join(progu)
  ggplot(abba.quad2) + 
    geom_point(aes(y=delta, x=step, shape=as.factor(kval)), size=.5) +
    geom_line(aes(y=delta, x=step, color=as.factor(fam_step5), group=as.factor(ordering)), size=.5) + 
    theme(legend.position = "bottom")
  
  
  
}

plotCumulativeUnion <- function(progu, title, summarize=TRUE, nshow=0, abba=FALSE){
  gcount=max(progu$ngen)
  norder=max(progu$ordering)
  datasets=unique(progu$dataset)
  
  alphas <- sapply(unique(progu$dataset),function(x){as.numeric(alpha(filter(progu,dataset==x)))}, USE.NAMES = TRUE)
print(alphas)
  summary <- summarize(group_by(progu, ngen, dataset), mean=mean(delta)) %>%
    mutate(alpha=alphas[dataset]) %>% 
    mutate(legend.name=paste0(dataset," (\u03b1=",round(alpha,3),")"))

  progu$opacity <- 1
  if (abba){
    most=max(progu$ngen)
    tmp <- progu %>% ungroup() %>%
      select(ordering, delta_pos, ngen) %>%
      group_by(ngen) %>%
      mutate(maxpos=max(delta_pos)) %>% ungroup() %>%
      mutate(lighten=ifelse(ngen< 4, ifelse(delta_pos == maxpos,ordering,NA ),NA))
    progu$opacity[progu$ordering %in% na.omit(tmp$lighten)] <- .5
  }
  
  tp <-
    ggplot() 
  
  if (nshow > 0){
    meanlinetype=c(2,3)
    progu <- progu %>%
      mutate(Indicator=ordering <= nshow) %>%
      # filter(ordering %in% c(1:20)) %>%
      mutate(ngen=as.integer(ngen),kval=as.factor(kval))
    print(names(progu))
    tp <- tp +
      geom_line(data=summary, aes(y=mean, x=ngen, linetype=legend.name)) +
    geom_point(data=filter(progu,Indicator),aes(y=delta, x=ngen, shape=kval), size=.5) +
    geom_line(data=filter(progu,Indicator),
              aes(y=delta, x=ngen,# linetype="Individual Ordering",
                  group=ordering, color=as.factor(ordering), alpha=opacity), size=.5) +
      scale_alpha_identity()+
      scale_shape_discrete(name="argmax(k)") +
      guides(color = 'none')
  }
  else{
    meanlinetype=c(1)
    print(names(summary))
    tp <- tp +
      geom_line(data=ungroup(summary), aes(y=mean, x=ngen, color=legend.name),linetype=c(1)) +
      scale_color_discrete(name = NULL) 
  }
  if (length(datasets) > 1){
    tp <- tp + 
      scale_linetype_manual(name=paste0("Mean (",norder," Orderings)")) 
  }
  
  tp <- tp + 
    #scale_linetype_manual(name=paste0("Mean (",norder," Orderings)")) +
    theme_bw() +
    scale_x_continuous(breaks= scales::pretty_breaks(10)) +
    # scale_color_discrete(name = "Random Genome Ordering") +
    xlab("Number of Genomes in Set") +
    ylab("Value of \u03b4*") +
    ggtitle(label=title)  +
    scale_y_continuous(labels = scales::label_number_auto()) +
    theme(legend.position = c(.75,.25))

  return(tp)
}






plot_scatter <- function (item, ytitle, filename, heap, alpha, param, plotlog=T) {
  param = latex2exp::TeX(param)
  if (abs(alpha) <= 1) {
    color_heap = "green"
  } else {
    color_heap = "red"
  }
  pl = ggplot(data = item, aes(x = x, y = y, color="Items")) +
    geom_point()+
    geom_line(data=heap,aes(x=x, y=y, color="Heap's law"), alpha=0.5)+
    theme_classic() +
    scale_color_manual(values=c(color_heap, "#000000"))+
    xlab("No. of genomes") +
    ylab(ytitle) + 
    ggtitle(param) +
    theme_bw(base_size = 16) + 
    theme(legend.text = element_text(size=5), plot.title =
            element_text(color="black", size=14))
  if (plotlog) {
    pl = pl + coord_trans(y='log10')#,x='log10')
  }
  plot(pl)
  #ggsave(filename=file.path(DIROUT, filename), plot=pl)
}



alphas <- function(progu){
  output <- sapply(unique(progu$dataset),
                   function(x){as.numeric(alpha(filter(progu,dataset==x)))}, USE.NAMES = TRUE)
return(output)
}


alpha <- function(progu.df){
  avg.progu <- ungroup(progu.df) 
  if (! "mean_delta" %in% colnames(progu.df)){
    avg.progu <- progu
  }
  # progu.df <- group_by (ngenmutate(progu.df, av)
  # if (! dataset %in% colnames(progu.df)){
  #   progu.df <- ungroup(mutate(progu.df, dataset="dataset"))
  # }
  avg.progu <- ungroup(progu.df) %>% 
    group_by(ordering) %>% 
    mutate(delta_delta=delta - lag(delta, default = delta[1])) %>%
    ungroup() %>% group_by(ngen) %>% 
    summarize(mean_delta=mean(delta), mean_delta_delta2 = mean(delta_delta)) %>%
    mutate(mean_delta_delta1 = mean_delta - lag(mean_delta, default = mean_delta[1]))
  new_item <- avg.progu$mean_delta_delta1
  N=length(new_item)
  # ALPHA
  x = 2:N
  model = lm(log(new_item[x])~log(x))
  alpha = abs(model$coefficients[2])
  return(alpha)
  }

graph_heaps <- function(){}

heaps <- function(progu.df){
  # progu.df <- group_by (ngenmutate(progu.df, av)
  avg.progu <- progu.df %>%
    group_by(ordering) %>% 
    mutate(delta_delta=delta_pos - lag(delta_pos, default = delta_pos[1])) %>%
    ungroup() %>% group_by(ngen) %>% 
    summarize(mean_delta=mean(delta_pos), mean_delta_delta2 = mean(delta_delta)) %>%
    mutate(mean_delta_delta1 = mean_delta - lag(mean_delta, default = mean_delta[1]))
  new_item <- avg.progu$mean_delta_delta1
  N=length(new_item)
  # ALPHA
  x = 2:N
  model = lm(log(new_item[x])~log(x))
  alpha = abs(model$coefficients[2])
  
  #str(summary(model))
  #summary(model)
  adj.r.sq = summary(model)$adj.r.squared
  x_h = seq(2,N,length.out=1000)
  y_h = exp(model$coefficients[1])*x_h^(model$coefficients[2])
  print(alpha)
  heap = data.frame(x=x_h,y=y_h)
  param = paste0("$\\alpha = $", signif(alpha, digits=4))
  plot_scatter(data.frame(x=2:N,y=new_item[2:N]), "No. of new items", 
               paste0("new_items.png"), heap, alpha, param)
  plot_scatter(data.frame(x=2:N,y=new_item[2:N]), "No. of new items", 
               paste0("new_items_log.png"), heap, alpha, param, plotlog=T)
  
}



clean_abba <- function(df, prefix="allvar_"){
  tmp <- df %>% tibble()
  tmp <- tmp %>% rowwise() %>%
    # make a column removing the brackets and spaces and quotes from the string(fasta list) column
    mutate(flist1=stringr::str_split(gsub("\\[|\\ |\\]|\\'","",fastas), pattern=",")) %>%
    # make a column stripping the file to the basename and removing expected extensions
    mutate(flist2=list(gsub(paste0(prefix,"|.fasta.gz|.fa|.fasta"),"",basename(flist1)))) %>%
    # mutate(flist3=list(lapply(flist1,function(z){gsub(paste0(prefix,"|.fasta.gz"),"",basename(z))}))) %>%
    # get the name of the last haplotype added
    mutate(last = unlist(flist2)[[ngen]]) %>%
    # get the sorted list of the previous ones (for combination in graph)
    mutate(prior = list(sort(flist2[1:(ngen-1)]))) %>%
    # make a sorted list of the whole list
    mutate(flist4=list(sort(flist2))) 
  
  # for each of the dataset names
  allsets <- unique(tmp$last)
  for (dataset in allsets){
    # assign all -1
    tmp[[dataset]] <- -1
    # reassign the ones that contain the dataset at issue to 1
    tmp[[dataset]][unlist(lapply(tmp$flist2, function(x){ dataset %in% x }))] <- 1
    # reassign (again!) whichever one is in last place
    tmp[[dataset]][tmp$last==dataset] <- 0
  }
  # rename and add desired columns
  tmp <- tmp %>% rename(step=ngen) %>% 
    group_by(ordering) %>% 
    mutate(deltadelta=delta-lag(delta,default=delta[2]),previous=lag(last))
  tmp$deltadelta[tmp$step == 1] <- NA
  return(tmp)
}
