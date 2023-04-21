using MondrianForests

#function plot_mondrian_tree(tree::MondrianTree)
    #@assert isa(tree, MondrianTree{2})
    #println(tree)
    #return nothing
#end

d = 2
lambda = 2.0
tree = MondrianTree(d, lambda)
#plot_mondrian_tree(tree)
MondrianForests.show(tree)
MondrianForests.show(tree)
