library(MEGENA)
library(Mergeomics)

################ data load
load("/HENA/hena/results/integration/integrated_int.RData")
load("/HENA/hena/results/integration/node_attributes.RData")

net <- integrated_int[abs(integrated_int$score) > 0.5, 1:3]
net2 <- integrated_int[abs(integrated_int$score) > 0.5,]
net$score <- abs(net$score)
net2[which(net2$score == 10),]$score <- 1
hist(net2$score)

# input parameters
n.cores <- 4; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "pearson" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 

#### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par){
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))}

##### calculate PFN
el <- calculate.PFN(net,doPar = doPar,num.cores = n.cores,keep.track = FALSE)
el_2 <- el
el_2$weight[1:29] <- 1
g <- graph.data.frame(el_2,directed = FALSE)

##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 10,max.size = vcount(g)/2,
                           doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                           save.output = FALSE)

###### unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# annot.table=NULL
annot.table <- na.omit(node_attributes[,1:2])
symbol.col <- 2
id.col <- 1
summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                       min.size = 10,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)


if (!is.null(annot.table)){
  # update annotation to map to gene symbols
  # V(g)$name <- s(V(g)$name,sep = "|")
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- summary.output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"}

print(head(summary.output$modules,2))
print(summary.output$module.table)

pnet.obj <- plot_module(output.summary = summary.output,PFN = g,subset.module = "c1_65",
                        layout = "kamada.kawai",label.hubs.only = TRUE,
                        gene.set = NULL,color.code =  "grey",
                        output.plot = FALSE,out.dir = "modulePlot",col.names = c("lavenderblush","lightsteelblue1","gray88"),label.scaleFactor = 20,
                        hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)


print(pnet.obj[[1]])


module.table <- summary.output$module.table
colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".

hierarchy.obj <- plot_module_hierarchy(module.table = module.table,label.scaleFactor = 0.15,
                                       arrow.size = 0.03,node.label.color = "blue")
#X11();
print(hierarchy.obj[[1]])


################ KDA
job.kda <- list()
job.kda$label<-"Alzh"
## parent folder for results
job.kda$folder<-"Results"
## Input a network
## columns: TAIL HEAD WEIGHT
job.kda$netfile <- "net.txt"

######### for tutorial
# job.kda$netfile<-system.file("extdata","network.mouseliver.mouse.txt",
#                              package="Mergeomics")
## Gene sets derived from ModuleMerge, containing two columns, MODULE, 
## NODE, delimited by tab 
# job.kda$modfile <- something_to_be_replaced
# job.kda$modfile<- system.file("extdata","mergedModules.txt",
#                               package="Mergeomics")
######### for tutorial

## "0" means we do not consider edge weights while 1 is opposite.
job.kda$edgefactor<-0.0
## The searching depth for the KDA
job.kda$depth<-1
## 0 means we do not consider the directions of the regulatory interactions
## while 1 is opposite.
job.kda$direction <- 0
job.kda$nperm <- 2000 # the default value is 2000, use 20 for unit tests

## kda.start() process takes long time while seeking hubs in the given net
## Here, we used a very small subset of the module list (1st 10 mods
## from the original module file):
# moddata <- tool.read(job.kda$modfile)

# modules <- summary.output$modules
# moddata_2 <- list()
# for(idx in 1:length(modules)){
#   moddata_2 <- rbind(moddata_2, cbind(names(modules[idx]), modules[[idx]]))
# }
# 
# mod.names <- unique(moddata$MODULE)[1:min(length(unique(moddata$MODULE)))]
# moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
# ## save this to a temporary file and set its path as new job.kda$modfile:
# tool.save(moddata, "subsetof.supersets.txt")
job.kda$modfile <- "subsetof.supersets.txt"

## Let's run KDA!

job.kda <- kda.configure(job.kda)
job.kda <- kda.start(job.kda)
job.kda <- kda.prepare(job.kda)
job.kda <- kda.analyze(job.kda)
job.kda <- kda.finish(job.kda)


## Remove the temporary files used for the test:
file.remove("subsetof.supersets.txt")
## remove the results folder
unlink("Results", recursive = TRUE)

graph <- job.kda$graph


res <- kda.finish.estimate(job.kda)
res <- kda.finish.trim(res, job.kda)
res <- kda.finish.summarize(res,job.kda)

cyto_job.kda <- job.kda
cyto_job.kda <- kda2cytoscape(cyto_job.kda, node.list = NULL, modules = NULL, ndrivers = 5, depth = 1)


save.image("kda_results.RData")
load("kda_results.RData")
