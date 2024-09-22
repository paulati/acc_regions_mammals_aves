base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/acc_regions_mammals_aves'

library(ape)
nm.base.path <- file.path(base_path, "mammals/conservation/data/neutral_model")
file.name <- "hg38.phastCons100way_label.nwk"
setwd(nm.base.path)
tree <- read.tree(file.name) # read tree, assign to variable tree.


plot(show.tip.label = TRUE, show.node.label = TRUE, tree)
     
plot(tree, type="cladogram") # you can plot this.
plot(tree, type="fan", show.tip.label = TRUE, show.node.label = TRUE) # you can plot this.
plot(tree, type="radial") # you can plot this.


#https://www.rdocumentation.org/packages/ape/versions/5.1/topics/plot.phylo

tree$node.label


