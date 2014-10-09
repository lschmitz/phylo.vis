# Introduction to phylogenetic trees in R

# Lars Schmitz, 2014

# 1 Loading libraries
# 2 Trees as "phylo" objects
# 3 Swapping sisterclades, identifying clades/tips, dropping tips
# 4 Fully resolved and polytomous trees
# 5 Modifying tree shape and other plotting options
# 6 Handling multiple trees



################################################################################################

# 1 Loading libraries
    
  # Libraries.
    
    library(ape)
    library(geiger)
    
################################################################################################
   

# 2 Trees as "phylo" objects
  
  # Let's begin by simulating a tree. There are many options for doing this. 
  # We'll use the rtree() function of the 'ape' package.
  # The tree is stored as an object called "phy"

    phy <- rtree(n=30) # n specifies the number of tips we want.

  # How does the tree look like? We can plot it with the following line:

    plot(phy)

  # So, now we wonder how the phylogenetic information is encoded.
  # First let's find out what the class of the object "phy" is.
  
    class("phy") # OK, it's a 'phylo' object.

  # By typing phy in the command line one can retrieve some basic information about the object.

    phy # A tree with n tips and (n-1) nodes, tip labels, the tree is rooted, and there are branch lengths.

  # But how is this information organized within the phylo object?
  # We can find out with the str() function, which displays the structure of an R object.

    str(phy)

  # The output tells us that a phylo object is a list of four components (+ two attributes): 
  # $edge, $tip.label, $edge.length, and $Nnode.

  # We can display these items separately

    phy$edge
    phy$tip.label
    phy$edge.length
    phy$Nnode

  # And one can also store these components as new objects:

    branches <- phy$edge
    species <- phy$tip.label
    brlength <- phy$edge.length
    nodes <- phy$Nnode

  # But what does this all mean? Let's explore this with a hand-written tree.
  # Assume you have a tree with 6 species (A through F).
  # The phylogenetic relationships of species A:F can be described in bracket form ("parenthetic format", or Newick):

    mini.phy <- read.tree(text = "((((A,B), C), (D,E)),F);")
    plot(mini.phy)

  # The structure of the object contains only three components this time because we didn't provide branch length:

    str(mini.phy)

  # But how are the edges defined? When we type...

    mini.phy$edge

  # ... we see a matrix of 10 rows and 2 columns. 
  # This matrix represents unique combinations of node- and tip numbers, defining each branch segment of the tree.

    plot(mini.phy, label.offset=0.2)  # the label.offset argument moves the species names a bit to the right
    nodelabels()                      # add node numbers
    tiplabels()                       # add tip numbers

  # For example, the branch (or edge) leading to species "C" is identified by row 6: 9, 3.

  # Both nodelabel() and tiplabel() function are quite flexible and afford the opportunity to visually anhance your tree plot.
  # Let's explore this further!
  
    ?tiplabels

  # And here's an example for mini.phy
  # We begin by creating an vector of different colors, with the same length as number of species in the tree.
  
    mycol<-c("blue", "blue", "blue", "red", "red", "red")

  # Now we plot the tree, moving the taxon names a bit to the right, and add the tiplabels without text, using symbols instead. 
    
    plot(mini.phy, adj=0, label.offset=0.75, lwd=2)
    tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=2)

################################################################################################


# 3 Swapping sisterclades, identifying clades/tips, dropping tips

  # It sometimes may be useful to rotate the tree about a specific node, i.e. swap sister clades.
  # This can be carried out with the rotate() function. Let's continue to work with mini.phy:

    plot(mini.phy, label.offset=0.2)  # the label.offset argument moves the species names a bit to the right
    nodelabels()                      # add node numbers
    tiplabels()                       # add tip numbers

  # How about we swap clades (D, E) and (A, B, C)? Their most recent common ancestor is found at node 8.
  
    rot.phy <- rotate(mini.phy, node=8)

  # And now let's see what happenend:

    plot(rot.phy, label.offset=0.2)   # the label.offset argument moves the species names a bit to the right
    nodelabels()                      # add node numbers
    tiplabels()                       # add tip numbers

  # It will also be very helpful to select all tips in a given clade.
  # This is implemented in the 'geiger' package; the tips() function finds all descendants of a node.
  
    cladeABC <- tips(rot.phy, node=9) # node 9 defines the clade composed of (A, B, C)
    cladeABC

  # Another helpful command allows for tree pruning, i.e. cutting of tips or entire clades.
  # For example, one can delete the entire cladeABC group:

    pruned.phy <- drop.tip(rot.phy, cladeABC)
    plot(pruned.phy, label.offset=0.2)    # the label.offset argument moves the species names a bit to the right
    nodelabels()                          # add node numbers
    tiplabels()                           # add tip numbers; note that nodes and tips are re-labeled.

  # Or we can drop tips (1 or multiple) randomly. Liam Revell explained how to do this nicely on his blog:
  # http://blog.phytools.org/2011/12/dropping-random-tip-or-set-of-tips-from.html

  # To prune tips, say m=2 random tips, enter:
  
    m=2
    pruned.phy2 <- drop.tip(rot.phy, sample(rot.phy$tip.label)[1:m]) # m=1 drops 1 single tip, of course!
    plot(pruned.phy2, label.offset=0.2)   
  
# It may also be useful to select all branches in a specific group of tips.
  # This is implemented in the 'ape' package; the which.edge() function finds all edges in a specified group.
  # For example, let's identify all branches of the cladeABC group as defined above.

    cladeABCbranches <- which.edge(rot.phy, cladeABC) # cladeABC was defined earlier, using the tips() function
    cladeABCbranches # this should be a numerical vector containing 6, 7, 8, 9
  
  # And as we can see, rows 6-9 of the $edge.length matrix represent the expected branches.
  # Let's first plot the tree, and then look at the $edge matrix for cross-checking.

    plot(rot.phy, label.offset=0.2)   # the label.offset argument moves the species names a bit to the right
    nodelabels()                      # add node numbers
    tiplabels()                       # add tip numbers
    rot.phy$edge
  
  # Here would be a good opportunity to show you how to assign different branch colors.
  # For example, how can one emphasize the branches of the clade formed by A, B, and C?
  # We first create a color vector, only consisting of grey colors.
  # Then we'll assign black to all branches of clade ABC.

    clcolr <- rep("darkgrey", dim(rot.phy$edge)[1]) 
    clcolr[cladeABCbranches] <- "black"
    plot(rot.phy, lwd=3, edge.color=clcolr)

  # Let's conclude this section with one last exercise: combining trees. 
  # Assume you have two different phylogenies, with two different sets of taxa (no overlap).
  # Another assumption is that you have knowledge how the trees may fit together.
  # Then the bind.tree() function of 'ape' package can help. 
  # The function takes in two phylo-objects.
  # The position of where the trees are bound is defined by tip- or node number within the first tree.
  # Note that you can also specify the "root" as binding position.

    tree1 <- rtree(n=10); plot(tree1); nodelabels(); tiplabels()
    tree2 <- rtree(n=10); plot(tree2)
    combined.tree <- bind.tree(tree1, tree2, where=1) 
    plot(combined.tree)


################################################################################################


# 4 Fully resolved and polytomous trees

  # All tree examples so far were fully resolved, i.e. each tree was fully binary. 
  # It's very easy to access visually for small trees, but one can also do this more formally:

    is.binary.tree(mini.phy)

  # Let's make a tree that is not fully resolved.

    poly.phy <- read.tree(text = "(((A,B,C),(D,E)),F);")
    plot(poly.phy)

  # Is this tree binary?

    is.binary.tree(poly.phy) # "FALSE"!

  # OK. Many comparative methods require fully resolved trees. But what to do if that's not the case?
  # The multi2di() function can resolve polytomies, either randomly or in the order in which they appear in the tree.
  # The default setting is to resolve polytomies randomly.

    resolved.phy <- multi2di(poly.phy)
    is.binary.tree(resolved.phy)
    
    plot(resolved.phy) # visual inspection.

  # If you repeat the above lines a few times you will see the effect of randomly resolving the polytomy.


################################################################################################


# 5 Modifying tree shape and other plotting options

  # There are many options for formatting and beautifying trees in R. Here are some basics.
  # Let's begin by simulating a tree once more. 

    phy <- rtree(n=10) # n specifies the number of tips we want.

  # The default plot produces a rightwards tree

    plot(phy)

  # The tree orientation can be changed by modifying the "direction"- argument. Try it out!

    plot(phy, direction="upwards") # other options are "rightwards" (default), "leftwards", and "downwards".

  # The font size of the tip labels (species names) can be changed with the cex argument.
  
    plot(phy, direction="upwards", 
         cex=5) # Try a few different settings!

# If we don't want the species names displayed, you can do the following:

    plot(phy, direction="upwards", 
         show.tip.label=FALSE)
    
  # And if you would like thicker branches, do this:

    plot(phy, direction="upwards", 
         show.tip.label=FALSE, 
         edge.width=2) # Try a few different settings!

  # Color of the branches can be controlled by the edge.color argument:

    plot(phy, direction="upwards", 
         show.tip.label=FALSE, 
         edge.width=2, 
         edge.color="blue") # Try a few different settings!
  
  # There are many other options... type ?plot.phylo in the console to read more.

  # Another useful command I'd like to introduce is the ladderize() function.
  # It reorganizes the tree structure, normally yielding much more readable trees.
  
    ladderized.phy <- ladderize(phy)

  # Let's create a plot with two panels, contrasting 'before" and 'after'.
  # With par() one can include the mfrow option, which specifies number of wors and columns in the plot matrix.

    par(mfrow=(c(1,2))) # 1 row but 2 plot panels

    plot(phy, direction="upwards", 
         show.tip.label=TRUE,
         edge.width=2,
         edge.color="blue")

    plot(ladderized.phy, direction="upwards", 
         show.tip.label=TRUE,
         edge.width=2,
         edge.color="blue")
    
    par(mfrow=c(1,1)) #for subsequent single-panel plots

  # Finally, let's change the type of the tree.
  # This is done by modifying the "type" argument. The defualt is set to "phylogram"
  # Other options include "cladogram", "fan", "unrooted", and "radial". Try all of them!

    plot(ladderized.phy, type="fan", 
         show.tip.label=FALSE,
         edge.width=2,
         edge.color="blue")
    
  # Tons of options!


################################################################################################


# 6 Handling multiple trees

  # In many cases you will be dealing with large distributions of trees, which are normally stored in a multiphylo-object.
  # In this section, we will read-in such a tree distribution, select a single tree, and do some exercises.

    haemulid.phy <- read.nexus("Haemulidae_trsample.nex")

  # What is the class of theis object?

    class(haemulid.phy) #multiPhylo

  # What does it contain?

    str(haemulid.phy)   # A list of 12 different trees, each with edge, edge.length, and node numbers

  # How can a single tree be selected?

    haemulid.phy.2 <- haemulid.phy[[2]] # Selecting the 2nd tree in the list

  # Using this tree, please do the following:

  # Drop Haemulon steindachneri Atlantic
  # Change plotting parameters to your liking, e.g. ladderize the tree, adjust the font size of the tip labels, 
  # optimize thickness of branches, etc. 

    