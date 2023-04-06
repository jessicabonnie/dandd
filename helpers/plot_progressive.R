#.libPaths(c("~/rlibs/4.0.2/gcc/9.3.0/","~/R/rstudio/4.0",.libPaths()))

#.libPaths(c(.libPaths(),"~/R/rstudio/4.0"))
# .libPaths(c("~/R/rstudio/4.0", .libPaths()))
#.libPaths(c(.libPaths(),"~/R/rstudio/4.0","~/rlibs/4.0.2/gcc/9.3.0"))
#.libPaths(c(.libPaths(),"~/rlibs/4.0.2/","~/rlibs/4.0.2/gcc/9.3.0","~/R/4.0.4","/usr/local/lib/R/site-library","/usr/local/lib/R/library"))

require(tidyr)
require(ggplot2)
require(data.table)
require(dplyr)
require(stringr)
require(stringi)
#require(latex2exp)

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
  #geom_point(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes, shape=kval), size=.5) +
  scale_linetype_manual(name=paste0("Fit (",norder," Orderings)"),values=c(2)) +
  geom_line(data=function(x) subset(x,Indicator),aes(y=delta_pos, x=ngenomes,
                                                     group=ordering, color=as.factor(ordering)), size=.5) +
  theme_bw() +
  scale_x_continuous(breaks= scales::pretty_breaks(10)) +
  #scale_shape_discrete(name="argmax(k)")+
  # scale_color_discrete(name = "Random Genome Ordering") +
  xlab("Number of Genomes in Set") +
  ylab("Value of \u03b4") +
  ggtitle(label=paste0("Values of \u03b4* over orderings of ",gcount," ", stringr::str_to_title(tag)," Genomes")) +
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

fixme <- function(x){ifelse(x==0,1,ifelse(x==-1,0,x))}

plotAbba <- function(progu, title="Title Here", stepnum=3, sdmult=0, nshow=NA, child=NA, relationship="parent"){
  subtitletext<-NA
  if (sdmult > 0){
  subtitletext=paste0("Subsetted deltadelta >/- mean +/- ",sdmult,"SD")
  sub.orderings <- progu %>% group_by(step)  %>%
    mutate(sd=sd(deltadelta), mean=mean(deltadelta)) %>%
    filter(deltadelta < mean-sdmult*sd | deltadelta > mean+sdmult*sd) %>%
    pull(ordering) %>% unique()
  } else{
    subtitletext<-NA
    sub.orderings <- progu %>% group_by(step)  %>%
      mutate(sd=sd(deltadelta), mean=mean(deltadelta)) %>%
      # filter(deltadelta < mean-sdmult*sd | deltadelta > mean+sdmult*sd) %>%
      pull(ordering) %>% unique()
  }
  
  if (!is.na(nshow)){
    sub.orderings <- head(sub.orderings,nshow)
  }
  maxorder=max(progu$ordering)
  average <- progu %>% group_by(step)  %>%
    summarize(delta=mean(delta))
  
  if (! is.na(child)){
    childsplit<-strsplit(child, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]]
    childpre=paste0(childsplit[1:length(childsplit)-1])
    childend=childsplit[length(childsplit)]
    padding=nchar(childend)
    parent2 <- paste0(childpre, str_pad(str=as.numeric(childend)-1,pad="0",width=padding, side="left"))
    mom <- paste0(childpre, str_pad(str=as.numeric(childend)-2,pad="0",width=padding, side="left"))
    matchparent<-matchstring<-paste0(child,"*|",mom,"*")
    matchstring<-paste0(child,"*|",mom,"*|",parent2,"*")
    haplomom <- paste0(mom,c("_1","_2"))
    haplomom <- haplomom[haplomom %in% names(progu)]
    haplodad <- paste0(parent2,c("_1","_2"))
    haplodad <- haplodad[haplodad %in% names(progu)]
    haploparents=c(haplomom,haplodad)
    print(mom)
    ordkey <- progu %>%
      select(ordering, step, flist4, matches(matchstring)) %>%
      filter(step == stepnum) %>%
      rowwise() %>%
      mutate(across(matches(matchstring),.fns=~fixme(.x))) %>%
      mutate(flist4=stringr::str_flatten(flist4)) %>% 
      mutate(fam=rowSums(across(matches(matchstring)))) %>% 
    mutate(mom=rowSums(across(contains(mom)))) %>% 
      mutate(parent=rowSums(across(matches(matchparent)))) %>%
      select(-matches(matchstring)) %>% 
      pivot_wider(names_from=step,values_from=c(flist4,fam,mom,parent), values_fn = as.factor, names_prefix = "step")
    progu <- inner_join(progu,ordkey)
  }
  
  progu2 <- progu %>% #inner_join(ordkey) %>% 
    mutate(orderingf=ifelse(ordering %in% sub.orderings, "outlier", "inside"))
  mean_text=paste0("Mean Over ",maxorder," Orderings")
  
  # ggplot(progu2) + 
  #   geom_line(aes(y=delta, x=step, color=as.factor(fam_step5), group=as.factor(ordering)),
  #             size=.5) + 
  #   theme(legend.position = "bottom")
  
  gg <- ggplot(filter(progu2, orderingf=="outlier"))  + 
    geom_line(data=average,aes(y=delta,x=step,linetype="mean"))
  
  if (!is.na(child)){
    color.var=paste0(relationship,"_step",stepnum)
    gg <- gg + geom_line(aes(y=delta, x=step, color=as.factor(.data[[color.var]]), 
                  group=as.factor(ordering), linetype="Individual Ordering"), size=.5)
  }
  else{
    gg <- gg + geom_line(aes(y=delta, x=step, color=as.factor(ordering),#color=as.factor(.data[[color.var]]), 
                             group=as.factor(ordering), linetype="Individual Ordering"), size=.5)
  }
  #+ ggtitle(subtitle=titletext)
  
  #color.var=paste0("fam_step",step)
  #geom_point(aes(y=delta, x=step, shape=as.factor(kval)), size=.5) +
    gg <- gg + theme_bw() + 
    theme(legend.position = "bottom") + 
    ggtitle(label = title,subtitle=subtitletext) +
      #paste0("Colored By Number of Related Haplotypes at Step ",step),subtitle=titletext) +
    #guides(color="none") +
    scale_linetype_manual(values = c("mean" = "solid", "Individual Ordering" = "dotted"),
                          labels=c("mean"=mean_text, "Individual Ordering" = "Individual Ordering"),
                          name=NULL)
  
  #color.var=paste0("parent_step",step)
  # gg + #geom_point(aes(y=delta, x=step, shape=as.factor(kval)), size=.5) +
  #   geom_line(aes(y=delta, x=step, #color=as.factor(.data[[color.var]]), 
  #                 group=as.factor(ordering)), size=.5) + theme_bw() + 
  #   theme(legend.position = "bottom") + ggtitle(paste0("Colored By Number of Parental Haplotypes at Step ",step),subtitle=titletext)
  # 
  print(sub.orderings)
  return(gg)
}



plotCumulativeUnion <- function(progu, title, summarize=TRUE, nshow=0, abba=FALSE){
  gcount=max(progu$ngen)
  norder=max(progu$ordering)
  datasets=unique(progu$dataset)
  
  alphas <- sapply(unique(progu$dataset),function(x){as.numeric(alpha(filter(progu,dataset==x)))}, USE.NAMES = TRUE)
print(alphas)
  summary <- summarize(group_by(progu, ngen, dataset), mean=mean(delta)) %>%
    mutate(alpha=alphas[dataset]) %>% 
    mutate(legend.print=paste0(dataset," (\u03b1=",round(alpha,3),")"))

  progu$opacity <- 1
  # if (abba){
  #   most=max(progu$ngen)
  #   tmp <- progu %>% ungroup() %>%
  #     select(ordering, delta_pos, ngen) %>%
  #     group_by(ngen) %>%
  #     mutate(maxpos=max(delta_pos)) %>% ungroup() %>%
  #     mutate(lighten=ifelse(ngen< 4, ifelse(delta_pos == maxpos,ordering,NA ),NA))
  #   progu$opacity[progu$ordering %in% na.omit(tmp$lighten)] <- .5
  # }
  
  tp <-
    ggplot() 
  
  if (nshow > 0){
    meanlinetype=c(2,3,4,5)[1:length(datasets)]
    progu <- progu %>%
      mutate(Indicator=ordering <= nshow) %>%
      # filter(ordering %in% c(1:20)) %>%
      mutate(ngen=as.integer(ngen),kval=as.factor(kval))
    print(names(progu))
    tp <- tp +
      geom_line(data=summary, aes(y=mean, x=ngen, linetype=legend.print)) +
    #geom_point(data=filter(progu,Indicator),aes(y=delta, x=ngen, shape=kval), size=.5) +
    geom_line(data=filter(progu,Indicator),
#              aes(y=delta, x=ngen,# linetype="Individual Ordering",
#                  group=ordering, color=as.factor(ordering), alpha=opacity), size=.5) +
#      scale_alpha_identity()+
              aes(y=delta, x=ngen,linetype="Individual Ordering",
                  group=ordering, color=as.factor(ordering)), size=.5) +
      #scale_shape_discrete(name="argmax(k)") +

      guides(color = 'none')
  }
  else{
    tp <- tp +
      geom_line(data=ungroup(summary), aes(y=mean, x=ngen, color=legend.print),linetype="mean") +
      scale_color_discrete(name = NULL) 
  }
  # if (length(datasets) > 1){
  #   tp <- tp + 
  #     scale_linetype_manual(name=paste0("Mean (",norder," Orderings)")) 
  # }
  
  tp <- tp + 

    #scale_linetype_manual(name=paste0("Mean (",norder," Orderings)")) +

    #guides(linetype=guide_legend(title=paste0("Mean (",norder," Orderings)"))) +
    # scale_linetype_manual(name=paste0("Mean (",norder," Orderings)")) +

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
  print(str(progu))
  output <- sapply(unique(progu$dataset),
                   function(x){as.numeric(alpha(filter(progu,dataset==x)))}, USE.NAMES = TRUE)
return(output)
}


alpha <- function(progu.df){
  avg.progu <- ungroup(progu.df) 
  if (! "mean_delta" %in% colnames(progu.df)){
    avg.progu <- progu.df
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
