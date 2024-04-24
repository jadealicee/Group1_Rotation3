#frequnecy tree

#making of a neighbour-joining tree
NJtree <- nj(dist(as.matrix(aa.genlight)))

#plot the tree
plot(NJtree, typ = "unrooted", show.tip = TRUE, cex = 0.7)

#add a title 
title(expression("Neighbour-Joining Tree of Arenosa and Lyrata"))

#save the tree
write.tree(NJtree), file = "NJ_Tree"
