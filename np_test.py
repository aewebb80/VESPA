import numpy as np

test_list = [2,4,6,8]

print np.average(test_list)

import dendropy

test_tree = dendropy.Tree.get_from_string("(A,(B,(C,D)));", "newick")
test_tree.print_plot()
