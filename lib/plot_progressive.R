
.libPaths(c(.libPaths(),"~/R/rstudio/4.0"))
.libPaths(c(.libPaths(),"~/rlibs/4.0.2/","~/rlibs/4.0.2/gcc/9.3.0","~/R/4.0.4","/usr/local/lib/R/site-library","/usr/local/lib/R/library"))

require(tidyr)
require(ggplot2)
require(data.table)
require(openssl)
require(dplyr)

# tag="HVSVC2_2"
# outdir="/home/jbonnie1/scr16_blangme2/jessica/dandd/progressive_union"
# gcount=10
# load(resultsfile)
library("optparse")

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

plotProgressiveUnion <- function(species, out, delta, nshow=Inf){
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
  geom_smooth(aes(y=delta_pos, x=ngenomes, linetype="Loess"), method='loess',formula=y~x) +
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
  ggtitle(label=paste0("Values of \u03b4 over orderings of ",gcount," ",stringr::str_to_title(tag)," Genomes")) +
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
