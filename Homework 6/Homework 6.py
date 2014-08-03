import comp182, FullBiTree, random, copy, collections

def write_newick(t):
    """Returns a Newick string corresponding to rooted,
    binary tree t.

    Arguments:
    t -- binary tree

    Returns:
    A Newick string."""

    
    if t.is_leaf():
        return str(t.get_node_property('taxon')) 
    else:
        nw = "(" #nw = newick string
        nw += write_newick(t.get_left_child())
        nw += ","
        nw += write_newick(t.get_right_child())
        nw += ")"
    return nw


def compute_nni_neighborhood(t):
    """Computes the set of all trees that can be obtained
    from t by a single nearest neighbor interchange.

    Arguments:
    t -- binary tree

    Returns:
    The set of all trees that can be a single NNI."""

    
    return helpercomputenni(t,t)

def helpercomputenni(t,n):
    """Function updates the nearest neighbor interexchange.

    Arguments:
    t -- a binary tree
    n -- a subtree of a binary tree

    Returns:
    An updated tree with the nearest neighbor interexchange."""
    fnni = nni(t,n) #fnni = Full nearest neighbor interchange
    if n.is_leaf():
        return fnni
    fnni.update(nni(t, n.get_left_child())) #checks the left side of tree
    fnni.update(nni(t, n.get_right_child())) # checks the right side of tree
    return fnni

def nni(t, n):
    """Interchanges the nodes on the binary tree.

    Arguments:
    t -- a binary tree
    n - a node on the binary tree

    Returns:
    All the changes at inputted node."""

    nc = set() #nc = node changes
    if tree_height(n) < 2:
        return nc
    nlc = n.get_left_child() #nlc = node left child
    nrc = n.get_right_child() #nrc = node right child

    if not nlc.is_leaf():
        nc.add(swapTrees(t, n, False, True))
        nc.add(swapTrees(t, n, False, False))
    if not nrc.is_leaf():
        nc.add(swapTrees(t, n, True, True))
        nc.add(swapTrees(t, n, True, False))
    return nc

def swapTrees(t, n, c, sc):
    """Swaps trees at the inputted node.

    Arguments:
    t - a binary tree
    n - node in binary tree
    c - child of the node
    sc - child of the child c

    Returns:
    A revised tree with changes in subtrees."""

    tc = copy.deepcopy(t) # tc = tree copy
    ln = locateNode(tc, n.get_name()) # ln = located node
    if c: #checks the child of the parent node
        left = ln.get_left_child()
        right = ln.get_right_child()
        if sc: #subchild = sc 
            ln.set_children(right.get_left_child(), right)
            right.set_children(left, right.get_right_child())
        else:
            ln.set_children(right.get_right_child(), right)
            right.set_children(right.get_left_child(), left)
    else: #if first case is not satisifed, does the reverse by replacing lefts with rights
        right = ln.get_right_child()
        left = ln.get_left_child()
        if sc:
            ln.set_children(left, left.get_left_child())
            left.set_children(right, left.get_right_child())
        else:
            ln.set_children(left, left.get_right_child())
            left.set_children(left.get_left_child(), right)
    return tc

def locateNode(t, node):
    """Locates the inputed node on the tree.

    Arguments:
    t -- a binary tree
    node -- a node on the binary tree

    Returns:
    The subtree located at the node."""

    if t.get_name() == node:
        return t
    if t.is_leaf():
        return None
    
    ls = locateNode(t.get_left_child(), node) #ls = left side of the tree
    if ls:
        return ls
    rs = locateNode(t.get_right_child(), node) #rs = right side of the tree
    if rs:
        return rs
    return None

def tree_height(tree):
    """
    Computes the height of a FullBiTree.

    Arguments:
    tree - a full binary tree in FullBiTree form.

    Returns:
    The height of the tree
    """
    if tree.is_leaf():
        return 0
    else:
        left_height = tree_height(tree.get_left_child())
        right_height = tree_height(tree.get_right_child())
        if left_height > right_height:
            return left_height + 1
        else:
            return right_height + 1



def random_tree(sot):
    """Returns a random evolutionary tree whose every leaf
    is labeled by a taxon name and the DNA sequence associated
    with it.

    Arguments:
    sot -- set of taxa

    Returns:
    A random evolutionary tree."""

    c = 0 #counter to apply a number reference to trees
    t = [] #places the trees into the list
    for x, y in sot.iteritems():

        fbt = FullBiTree.FullBiTree('leaf'+str(c)) #creates a single binary tree
        fbt.set_node_property('taxon', x)
        fbt.set_node_property('sequence', y)
        t.append(fbt)
        c = c + 1
    return connectTrees(t)

def connectTrees(t):
    """Takes a list of single noded trees and creates one tree.

    Arguments:
    t -- list of trees

    Returns:
    A single tree."""

    c = 0

    while len(t) > 1:
        rtf = random.randrange(len(t)) # rtf = random tree first
        ftf = t.pop(rtf) # ftf = first tree found
        rtl = random.randrange(len(t)) # rtl = random tree last
        ftl = t.pop(rtl) # ftl = first tree last
        mt = FullBiTree.FullBiTree('node'+str(c), ftf, ftl) #mt = make tree
        mt.set_node_property('taxon', None)
        mt.set_node_property('sequence', None)
        t.append(mt)
        c = c + 1
    return t[0]



def compute_ps(t, sek, sL):
    """Computes the PS in three phases.

    Arguments:
    t -- binary tree
    sek -- inputted data of sequence
    sL -- the length of the sequence

    Returns:
    The PS score based on the equation."""

    

    if t.is_leaf():
        return "Tree needs to be a full binary tree."
    sk = "sets_key"
    
    helperBottomUp(t, sek, sk, sL)
    
    roots = ""
    rootset = t.get_node_property(sk)

    for x in xrange(sL): # applies the properties to the roots
        rootidx = rootset[x]
        roots += rootidx.pop()
    t.set_node_property(sek, roots)

    

    return helperTopDown(t, t.get_left_child(), sek, sk, sL) + helperTopDown(t, t.get_right_child(), sek, sk, sL)

def helperTopDown(preNode, child, seqKey, setK, sL):
    """Starts at the top to calculate scores.

    Arguments:
    preNode -- the previous node
    child -- child of the previous node
    seqKey -- inputted data of sequence
    sL -- Sequence length
    setK -- property key to apply a taxa

    Returns:
    Computed scores."""

    counter = 0
    preSeq = preNode.get_node_property(seqKey)
    
    if not child.is_leaf():
        counter = 0
        cSeq = "" #child sequence
        childPro = child.get_node_property(setK)
        
        for x in xrange(sL): #looks for similariites in the sequences
            childProIn = childPro[x]
            preSeqIn = preSeq[x]
            if preSeqIn in childProIn:
                cSeq += preSeqIn
            else:

                cSeq += childProIn.pop()
                counter = counter + 1
                

        child.set_node_property(seqKey, cSeq)
        counter = counter + helperTopDown(child, child.get_left_child(), seqKey, setK, sL) + helperTopDown(child, child.get_right_child(), seqKey, setK, sL)
    else:
        cSeq = child.get_node_property(seqKey)
        counter = sum(a != b for a, b in zip(cSeq, preSeq))
    return counter

    
def helperBottomUp(t, seqK, setK, sL):
    """Starts at the down to calculate scores and makes node
    assignments.

    Arguments:
    t -- a binary tree
    seqK -- inputted data of sequence
    sL -- Sequence length
    setK -- property key to apply a taxa


    Returns:
    Computed scores and node assignments."""

    sets = []
    #print t, "t in helperbottomup"
    t.set_node_property(setK, sets)
    #print "First run through of bottomup"

    if t.is_leaf():
        seq = t.get_node_property(seqK)
        #print "loop goes here line 305"
        for x in range(sL):
            mSet = set()
            sets.append(mSet)
            mSet.add(seq[x])
    else:
        #print "loop goes here line 311"
        glc = t.get_left_child()
        #print "left child line 313", glc
        grc = t.get_right_child()
        #print "right child line 315", grc
        helperBottomUp(glc, seqK, setK, sL)
        #print "function does not make it to line 317"
        helperBottomUp(grc, seqK, setK, sL)

        propGlc = glc.get_node_property(setK)
        propGrc = grc.get_node_property(setK)

        #helperBU(sL, sets, propGlc, propGrc, mSet)
        #print "ran through bottomuphelper"
        


def infer_evolutionary_tree(starts, inseq, txtfile):
    """Creates a evolutionary tree.

    Arguments:
    inseq -- a set of DNA sequences and the length of each sequence
    txtfile -- list of taxa to write Newick string
    starts -- integer to specify the number of restarts

    Returns:
    An evolutionary tree on the number of inputted sequences."""

    m, seq = read_phylip(inseq)
    steps = []
    startcompare = {'t': float('inf'), 'score': float('inf')}
    run = 1
    for x in range(starts):
        print "Run Number:", run
        run += 1
        step = 0
        minim = False
        t = random_tree(seq)
        score = compute_ps(t, 'sequence', m)
        while minim == False:
            check = compute_nni_neighborhood(t)
            compare = {'t': t, 'score': score}
            print "Current Found Score:", score
            for x in check:
                z = compute_ps(x, 'sequence', m)
                if z < compare['score']:
                    compare['t'] = x
                    compare['score'] = z
            if compare['t'] == t:
                minim = True
            else:
                t = compare['t']
                score = compare['score']
            step += 1
        steps.append(step)
        print "Steps", step
        print "Best Score:", compare['score']
        if compare['score'] < startcompare['score']:
            startcompare['t'] = compare['t']
            startcompare['score'] = compare['score']
    print "Best Score in all runs:", startcompare['score']
    ne = write_newick(startcompare['t'])
    with open(txtfile, 'w') as f:
        f.write(ne)
    print inseq, "Mean steps taken:", float(sum(steps)) / len(steps)
        
        
    

def read_phylip(filename):
    """
    Read a file in Phylip format and return the length of the
    sequences and the taxa and sequences.

    Arguments:
    filename -- name of file in Phylip format

    Returns:
    A tuple where the first element is the length of the sequences and
    the second argument is a dictionary mapping taxa to sequences.
    """
    # Initialize return values in case file is bogus
    m = 0
    tsmap = {}

    with open(filename) as f:
        # First line contains n and m separated by a space
        nm = f.readline()
        nm = nm.split()
        n = int(nm[0])
        m = int(nm[1])

        # Subsequent lines contain taxon and sequence separated by a space
        for i in range(n):
            l = f.readline()
            l = l.split()
            tsmap[l[0]] = l[1]

    # Return sequence length and mapping of taxa to sequences
    return m, tsmap


### The code below allows you to turn a Newick string into a 
### FullBiTree.  You may find this useful for testing, but it
### is not required for the homework.  To use this, you will
### need to download and install biopython:
###  http://biopython.org/wiki/Download

from Bio import Phylo
from cStringIO import StringIO
import uuid
import FullBiTree

def parse_newick(newickstr, taxon_name_key):
    """
    Creates a FullBiTree representation of a newick string.

    Arguments:
    newickstr - a newick string
    taxon_name_key - the property key that will be used to apply a taxon's name to a leaf's node
                     property when creating FullBiTree leafs.

    Returns:
    A FullBiTree representation of the given newick string.
    """
    seen_names = set()
    tree = Phylo.read(StringIO(newickstr), "newick")
    return process_clade(tree.root, taxon_name_key, seen_names)

def process_clade(clade, taxon_name_key, seen_names):
    """
    Creates a FullBiTree representation of a clade.

    Arguments:
    clade - a clade object as contructed by Phylo.read(...)
    seen_names - a set of string names of each node previously
             processed during construction.
    taxon_name_key - the property key that will be used to apply a taxon's name to a leaf's node
                     property when creating FullBiTree leafs.

    Returns:
    A FullBiTree representation of the given clade.
    """

    if (len(clade.clades) == 0) and (clade.name is None or len(clade.name) < 1):
        # Leaf nodes have to have a name (the taxon)
        raise Exception("Leaf node must have a name.")

    if clade.name is None or len(clade.name) < 1:
        # Give the internal node a unique name
        clade.name = str(uuid.uuid4())

    if len(clade.clades) != 0 and len(clade.clades) != 2:
        raise Exception("Each tree node must have zero or two children.")

    if clade.name in seen_names:
        raise Exception("Every node name in the tree must be unique. " + 
                        clade.name + " is duplicated.")

    seen_names.add(clade.name)

    if len(clade.clades) == 0:
        # This is a leaf node
        tree = FullBiTree.FullBiTree(clade.name)
        tree.set_node_property(taxon_name_key, clade.name)
        return tree
    else:
        # This is an internal node
        left  = process_clade(clade.clades[0], taxon_name_key, seen_names)
        right = process_clade(clade.clades[1], taxon_name_key, seen_names)
        tree = FullBiTree.FullBiTree(clade.name, left, right)
        return tree

def run_files():
    """Runs all the files with sequences. 
    Arguments:
    None

    Returns:
    Results from the run of the sequences.
    """
    infer_evolutionary_tree(50, 'primate_seqs.phylip', 'primate_newick.txt')
    infer_evolutionary_tree(50, 'yeast_gene1_seqs.phylip', 'yeast1_newick.txt')
    infer_evolutionary_tree(50, 'yeast_gene2_seqs.phylip', 'yeast2_newick.txt')
