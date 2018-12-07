/* 
 * PhyloTree.java
 *
 * Defines a phylogenetic tree, which is a strictly binary tree 
 * that represents inferred hierarchical relationships between species
 * 
 * There are weights along each edge; the weight from parent to left child
 * is the same as parent to right child.
 *
 * 
 * Authors: I.C. & I.M.
 *
 */
 
import java.io.File;
import java.util.Stack;
import java.util.Arrays; 
import java.util.Scanner;
import java.util.HashMap;
import java.util.ArrayList;
import java.io.FileNotFoundException;


public class PhyloTree {
    private PhyloTreeNode overallRoot;    // The actual root of the overall tree
    private int printingDepth;            // How many spaces to indent the deepest 
                                          // node when printing    
    // CONSTRUCTOR

    // PhyloTree
    // Pre-conditions:
    //        - speciesFile contains the path of a valid FASTA input file
    //        - printingDepth is a positive number
    // Post-conditions:
    //        - this.printingDepth has been set to printingDepth
    //        - A linked tree structure representing the inferred hierarchical
    //          species relationship has been created, and overallRoot points to
    //          the root of this tree
    // Notes:
    //        - A lot happens in this step!  See assignment description for details
    //          on the input format file and how to construct the tree
    //        - If you encounter a FileNotFoundException, print to standard error
    //          "Error: Unable to open file " + speciesFilename
    //          and exit with status (return code) 1
    //        - Most of this should be accomplished by calls to loadSpeciesFile and buildTree
    public PhyloTree(String speciesFile, int printingDepth) {
                
        this.printingDepth= printingDepth;
        Species[] species = loadSpeciesFile(speciesFile);

        buildTree(species);          
        
        return;        
    }

    // ACCESSORS

    // getOverallRoot
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the overall root
    public PhyloTreeNode getOverallRoot() {
        return this.overallRoot;
    }

    // toString 
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns a string representation of the tree
    // Notes:
    //    - See assignment description for proper format
    //        (it will be a kind of reverse in-order [RNL] traversal)
    //    - Can be a simple wrapper around the following toString
    //    - Hint: StringBuilder is much faster than repeated concatenation
    public String toString() {
        return toString(this.overallRoot, 0.0, getWeightedHeight());
    }

    // toString 
    // Pre-conditions:
    //    - node points to the root of a tree you intend to print
    //    - weightedDepth is the sum of the edge weights from the
    //      overall root to the current root
    //    - maxDepth is the weighted depth of the overall tree
    // Post-conditions:
    //    - Returns a string representation of the tree
    // Notes:
    //    - See assignment description for proper format
    private String toString(PhyloTreeNode node, double weightedDepth, double maxDepth) {        
        
        weightedDepth = getWeightedDepth(node);
        String path = "";
        StringBuilder temp = new StringBuilder();
        int k = (int)(this.printingDepth * (weightedDepth/ maxDepth));
        
        if(node.getRightChild() != null)
            temp.append(toString(node.getRightChild(), weightedDepth, maxDepth));
                           
        for(int s = 0; s < k; s++)
            path+= ".";
        
        temp.append(path);
        temp.append(node.toString() + "\n");
                
        if(node.getLeftChild() != null)
            temp.append(toString(node.getLeftChild(), weightedDepth, maxDepth));                
        
        return temp.toString();
    }
    
    //getWeightedDepth
    //Pre-conditions:
    //  - node is the start reference to start summing weight along 
    //    some path to the overallRoot
    //Post-conditions:
    //  - returns the sum of weights along path
    private double getWeightedDepth(PhyloTreeNode node) {                
        
        double totalSum = 0.0;            
        while(node != this.overallRoot){
            node = node.getParent();
            totalSum += node.getDistanceToChild();
        }
        return totalSum;
    }       
    
    // toTreeString 
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns a string representation in tree format
    // Notes:
    //    - See assignment description for format details
    //    - Can be a simple wrapper around the following toTreeString
    public String toTreeString() {
        return toTreeString(this.overallRoot);
    }
    
    // toTreeString 
    // Pre-conditions:
    //    - node points to the root of a tree you intend to print
    // Post-conditions:
    //    - Returns a string representation in tree format
    // Notes:
    //    - See assignment description for proper format
    private String toTreeString(PhyloTreeNode node) {
        
        String temp = "";
        double tempVal = 0.0;
        StringBuilder sBuilder = new StringBuilder();
        
        if(node.isLeaf() ){            
            tempVal = Double.valueOf(node.getParent().getDistanceToChild());
            sBuilder.append(":" + String.format("%.5f", tempVal) );
       
        } else {            
            sBuilder.append("(" + toTreeString(node.getRightChild() ) + ", " + toTreeString(node.getLeftChild() ) + ")");            
            if(node != this.overallRoot){            
                tempVal = Double.valueOf(node.getParent().getDistanceToChild());
                sBuilder.append(":" + String.format("%.5f", tempVal) );
            }
        }
        return sBuilder.toString();
    }
        

    // getHeight
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the tree height as defined in class
    // Notes:
    //    - Can be a simple wrapper on nodeHeight
    public int getHeight() {
        return nodeHeight(this.overallRoot) ;
    }
 
    // getWeightedHeight
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the sum of the edge weights along the
    //      "longest" (highest weight) path from the root
    //      to any leaf node.
    // Notes:
    //   - Can be a simple wrapper for weightedNodeHeight
    public double getWeightedHeight() {
        return weightedNodeHeight(this.overallRoot);         
    }
    
    // countAllSpecies
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the number of species in the tree
    // Notes:
    //    - Non-terminals do not represent species
    //    - This functionality is provided for you elsewhere
    //      just call the appropriate method
    public int countAllSpecies() {    
        return countSpecies();
    }

    //countSpecies
    //Pre-conditions:
    //  - None
    //Post-conditions:
    //  - Return the number of terminal nodes, or leafs, 
    //  - of the the tree rooted at this.overallRoot.    
    private int countSpecies(){                
        
        int counter = 0;
        PhyloTreeNode node;
        Stack<PhyloTreeNode> stack = new Stack<PhyloTreeNode>();        
        stack.add(this.overallRoot);
        
        while( !stack.isEmpty()){
            node = stack.pop();
                   
            if( !node.getLeftChild().isLeaf())            
                stack.push(node.getLeftChild());                            
            else 
                counter++;            
            if( !node.getRightChild().isLeaf())
                stack.push(node.getRightChild());
            else
                counter++;                                                                  
        }                        
        return counter;
    }        
    
    // getAllSpecies
    // Pre-conditions:a tree you intend to print
    // Post-conditions:
    //    - None
    // Post-conditions:
    //    - Returns an ArrayList containing all species in the tree
    // Notes:
    //    - Non-terminals do not represent species
    public java.util.ArrayList<Species> getAllSpecies() {  
    
        ArrayList<Species> species = new ArrayList<Species>();                
        getAllDescendantSpecies(this.overallRoot, species);                
        return species;
    }    

    // findTreeNodeByLabel
    // Pre-conditions:
    //    - label is the label of a tree node you intend to find
    //    - Assumes labels are unique in the tree
    // Post-conditions:
    //    - If found: returns the PhyloTreeNode with the specified label
    //    - If not found: returns null
    public PhyloTreeNode findTreeNodeByLabel(String label) {          
        return findTreeNodeByLabel(this.overallRoot, label);
    }
            
    // findLeastCommonAncestor
    // Pre-conditions:
    //    - label1 and label2 are the labels of two species in the tree
    // Post-conditions:
    //    - If either node cannot be found: returns null
    //    - If both nodes can be found: returns the PhyloTreeNode of their
    //      common ancestor with the largest depth
    //      Put another way, the least common ancestor of nodes A and B
    //      is the only node in the tree where A is in the left tree
    //      and B is in the right tree (or vice-versa)
    // Notes:
    //    - Can be a wrapper around the static findLeastCommonAncestor
     public PhyloTreeNode findLeastCommonAncestor(String label1, String label2) {
    
        PhyloTreeNode node;    
        PhyloTreeNode node1 = findTreeNodeByLabel(this.overallRoot, label1);
        PhyloTreeNode node2 = findTreeNodeByLabel(this.overallRoot, label2);
        int depth1, depth2;
        if(node1.getParent() == node2.getParent())
            return node1.getParent();
        
        Stack<PhyloTreeNode> stack = new Stack<PhyloTreeNode>();
        
        if(node1 == null || node2 == null)            
            return null;
        
        while( node1 != node2 ){
            depth1 = nodeDepth(node1);
            depth2 = nodeDepth(node2);
            if(depth2 > depth1)
                node2 = node2.getParent();
            else if(depth1 > depth2)
                node1 = node1.getParent();
            else{
                node1 = node1.getParent();
                node2 = node2.getParent();
            }
        }                               
        return node2;
    }
    
    // findEvolutionaryDistance
    // Pre-conditions:
    //    - label1 and label2 are the labels of two species in the tree
    // Post-conditions:
    //    - If either node cannot be found: returns POSITIVE_INFINITY    
    //    - If both nodes can be found: returns the sum of the weights 
    //      along the paths from their least common ancestor to each of
    //      the two nodes
    public double findEvolutionaryDistance(String label1, String label2) {
        
        PhyloTreeNode node1 = findTreeNodeByLabel(label1);
        PhyloTreeNode node2 = findTreeNodeByLabel(label2);
        
        if( node1 == null || node2 == null) 
            return java.lang.Double.POSITIVE_INFINITY;
        else{
            PhyloTreeNode ancestor = findLeastCommonAncestor(node1, node2);
            
            double distance = 0.0;
            distance += getWeightedDistance(ancestor, node1);
            distance += getWeightedDistance(ancestor, node2);        
            return distance;
        }
    }

    // MODIFIER

    // buildTree
    // Pre-conditions:
    //    - species contains the set of species for which you want to infer
    //      a phylogenetic tree
    // Post-conditions:
    //    - A linked tree structure representing the inferred hierarchical
    //      species relationship has been created, and overallRoot points to
    //      the root of said tree
    // Notes:
    //    - A lot happens in this step!  See assignment description for details
    //      on how to construct the tree.
    //    - Be sure to use the tie-breaking conventions described in the pdf
    //    - Important hint: although the distances are defined recursively, you
    //      do NOT want to implement them recursively, as that would be very inefficient
    private void buildTree(Species[] species) {
        
        int forestSize = 0;
        double minDistance;
        double tempMin = 0.0;
        int lPointer, rPointer;        
        String otherTree="";
        double constant1, constant2, saveConstant1, computedDistance;
        ArrayList<PhyloTreeNode> forest = new ArrayList<PhyloTreeNode>();
        MultiKeyMap<Double> distances = new MultiKeyMap<Double>();
        int numSpecies = species.length;        
        
        //step 1: creating forest of terminal nodes, or trees
        for(int m =0; m< numSpecies; m++)
            forest.add(new PhyloTreeNode(null, species[m]));
                
        //step 2: computing pairwise distance
        String species1, species2, tree1, tree2;
        for(int n =0; n< numSpecies; n++){
            species1= species[n].getName();
            for(int p=n+1; p< numSpecies    ; p++){
                species2 = species[p].getName();
                if(distances.get(species1, species2) == null)
                    distances.put(species1, species2, Species.distance(species[n], species[p]));                
            }
        }        
        
        //step 3: merging trees
        while(forest.size()>1){            
            lPointer=0;
            rPointer=0;
            
            //step 3.a
            minDistance = java.lang.Double.POSITIVE_INFINITY;  
            forestSize= forest.size();
            
            for(int i=0; i<forestSize; i++){                
                tree1= forest.get(i).getLabel();
                for(int j=0; j< forestSize; j++){
                
                    tree2 = forest.get(j).getLabel();
                    if(distances.get(tree1, tree2) != null){
                        tempMin= distances.get(tree1, tree2);
                        if(tempMin< minDistance && !tree1.equals(tree2)){
                            minDistance= tempMin;
                            lPointer = i;
                            rPointer=j;
                            
                            if(tree1.compareTo(tree2) < 0){                                
                                lPointer =i;               
                                rPointer =j;
                            }
                            else{
                                lPointer = j;
                                rPointer = i;                            
                            }                                                        
                        }
                        else if(tempMin==minDistance && tree1.compareTo(tree2)>0){
                            lPointer =j;
                            rPointer =i;
                        }
                    }
                }
            }
            
            //step 3.b

            String lName = forest.get(lPointer).getLabel();
            String rName = forest.get(rPointer).getLabel();            
            
            String mergerName = lName + "+" + rName;            
            PhyloTreeNode mergerTree = new PhyloTreeNode( mergerName, null, forest.get(lPointer), forest.get(rPointer), distances.get(lName, rName)/2.0);
            
            PhyloTreeNode left;
            PhyloTreeNode right;
            
            if(rPointer> lPointer){
                right = forest.remove(rPointer);   
                left = forest.remove(lPointer);
            }else{
                left = forest.remove(lPointer);                        
                right = forest.remove(rPointer);                           
            }            
            
            mergerTree.getLeftChild().setParent(mergerTree);
            mergerTree.getRightChild().setParent(mergerTree);

            forest.trimToSize();
            forest.add(mergerTree);
            forestSize--;
            if(forest.size() ==1) 
                this.overallRoot = mergerTree;
                        
            constant1 = (double) mergerTree.getLeftChild().getNumLeafs();
            constant2 = (double) mergerTree.getRightChild().getNumLeafs();        
            saveConstant1 = constant1;
            constant1 = constant1/ (constant1 + constant2);
            constant2 = constant2/ (saveConstant1 + constant2);

            for(int q =0; q< forestSize; q++){
                otherTree = forest.get(q).getLabel();
                if(distances.get(otherTree, lName) != null && distances.get(otherTree, rName) != null){
                    computedDistance = constant1*distances.get(otherTree, lName) + constant2*distances.get(otherTree, rName);
                    distances.put(mergerName, otherTree, computedDistance);
                }
            }                      
        }//end while loop            
        
    }
    // STATIC

    //getWeightedDistance  HELPER
    //Pre-conditions:
    //  - rootPnter is the root of the tree or subtree with terminal 
    //  - node leafPnter.
    //Post-conditions:
    //  - Returns the weight of the path from rootPnter to leafPnter.
    private static double getWeightedDistance(PhyloTreeNode rootPnter, PhyloTreeNode leafPnter){
        double distance=0.0;
        
        while(leafPnter != rootPnter){            
            leafPnter = leafPnter.getParent();
            distance += leafPnter.getDistanceToChild();
        }
        return distance;
    }    
    
    //getHeight  HELPER
    //Pre-conditions:
    //  - node is the node to calculate height of.
    //Post-conditions:
    //  - Returns the the height of node plus the 
    //  - maximum height of either the left or right
    //  - subtree of node.
    private static int getHeight(PhyloTreeNode node){
        
        if(node == null) 
            return -1;
        if(node.isLeaf())
            return 0;
        else{
            int left, right;
            left = getHeight(node.getLeftChild());
            right = getHeight(node.getRightChild());
            if( left>= right)
                return 1 + left;
            else 
                return 1 + right;    
        }
    }    
    
    // nodeDepth
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns -1
    //    - Else: returns the depth of the node within the overall tree
    public static int nodeDepth(PhyloTreeNode node) {
        if(node == null) 
            return -1;
        if(node.getParent()== null) 
            return 0;
        return 1 + nodeDepth(node.getParent());
    }
    
    // nodeHeight
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns -1
    //    - Else: returns the height subtree rooted at node
    
    public static int nodeHeight(PhyloTreeNode node) {
        if(node == null)            
            return -1;        
        else 
            return getHeight(node);      
    }          
        
    // weightedNodeHeight 
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns NEGATIVE_INFINITY
    //    - Else: returns the weighted height subtree rooted at node
    //     (i.e. the sum of the largest weight path from node
    //     to a leaf; this might NOT be the same as the sum of the weights
    //     along the longest path from the node to a leaf)
    public static double weightedNodeHeight(PhyloTreeNode node) {
        
        if(node == null)
            return java.lang.Double.NEGATIVE_INFINITY;
        else{
            double maxDistance = 0.0;
            double tempDistance = 0.0;
            PhyloTreeNode leafNode;
            ArrayList<Species> species = new ArrayList<Species>();
            getAllDescendantSpecies(node, species);
            
            for(Species leaf: species){
                leafNode = findTreeNodeByLabel(node, leaf.getName());                
                tempDistance = getWeightedDistance(node, leafNode);
                if(tempDistance > maxDistance)
                    maxDistance = tempDistance;
            }                        
            return maxDistance;
        }
    }
    
    // loadSpeciesFile
    // Pre-conditions:
    //    - filename contains the path of a valid FASTA input file
    // Post-conditions:
    //    - Creates and returns an array of species objects representing
    //      all valid species in the input file
    // Notes:
    //    - Species without names are skipped
    //    - See assignment description for details on the FASTA format
    // Hints:
    //    - Because the bar character ("|") denotes OR, you need to escape it
    //      if you want to use it to split a string, i.e. you can use "\\|"       
   public static Species[] loadSpeciesFile(String filename) {
    
        Scanner input;
        int listSize = 0;
        int tempSize = 0;
        String token = "";
        Species[] species;        
        String[] tempString;
        String[] tempSequence;        
        ArrayList<Species> speciesList = new ArrayList<Species>();   
        boolean firstline = true;        
        StringBuilder temp = new StringBuilder();        
        StringBuilder sequence = new StringBuilder();
        
        try{
           input = new Scanner (new File (filename));
           while(input.hasNext()){            
                token = input.nextLine();                        
                
                if(token.charAt(0) == '>' && firstline==false){
                
                    tempSize = temp.length();
                    tempString = sequence.deleteCharAt(0).toString().split("\\|");
                    tempSequence = new String[tempSize];
                    
                    for(int i = 0; i< tempSize; i++)
                        tempSequence[i] = temp.charAt(i) + "";                                    
                                            
                    speciesList.add(new Species( tempString[tempString.length-1],tempSequence) ); 
                    temp.delete(0,tempSize);
                    sequence.delete(0,sequence.length());
                    sequence.append(token);
                }
                
                if(token.charAt(0) != '>')                
                    temp.append(token);                 
                else {
                    firstline = false;
                    sequence.append(token);
                }
            }
            tempSize = temp.length();
            tempString = sequence.deleteCharAt(0).toString().split("\\|");                        
            tempSequence = new String[tempSize];
                    
            for(int i = 0; i< tempSize; i++)
                tempSequence[i] = temp.charAt(i) + "";
            
            speciesList.add(new Species( tempString[tempString.length-1],tempSequence) ); 
            species = new Species[speciesList.size()];
            
            listSize = speciesList.size();
            for(int i = 0; i< listSize; i++)
                species[i] = speciesList.get(i);
                        
            return species;         
        } catch(FileNotFoundException e){
            System.out.println("Error: 1 " + filename + " not found.");
            System.exit(1);        
        }
        return null;
    }

    
    // getAllDescendantSpecies
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - descendants is a non-null reference variable to an empty arraylist object
    // Post-conditions:
    //    - descendants is populated with all species in the subtree rooted at node
    //      in in-/pre-/post-order (they are equivalent here)
    private static void getAllDescendantSpecies(PhyloTreeNode node, java.util.ArrayList<Species> descendants) {
        
        if(node == null)
            return;     
        
        PhyloTreeNode nextNode;        
        Stack<PhyloTreeNode> stack = new Stack<PhyloTreeNode>();
        
        stack.push(node);
        
        while( !stack.isEmpty()){
            nextNode = stack.pop();
            
            if( !nextNode.getLeftChild().isLeaf())
                stack.push(nextNode.getLeftChild());
            else
                descendants.add(nextNode.getLeftChild().getSpecies());
            if(!nextNode.getRightChild().isLeaf())
                stack.push(nextNode.getRightChild());
            else
                descendants.add(nextNode.getRightChild().getSpecies());            
        }
        return;
    }

    // findTreeNodeByLabel
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - label is the label of a tree node that you intend to locate
    // Post-conditions:
    //    - If no node with the label exists in the subtree, return null
    //    - Else: return the PhyloTreeNode with the specified label 
    // Notes:
    //    - Assumes labels are unique in the tree
    private static PhyloTreeNode findTreeNodeByLabel(PhyloTreeNode node, String label) {    
    
        PhyloTreeNode searchNode;        
        Stack<PhyloTreeNode> stack = new Stack<PhyloTreeNode>();        
        stack.add(node);        
       
        while( !stack.isEmpty()){
            searchNode = stack.pop();
            
            if( label.equals(searchNode.getLabel()))
                return searchNode;            
            else{
                if( searchNode.getLeftChild() !=null)
                    stack.push(searchNode.getLeftChild());           
                if (searchNode.getRightChild() != null)
                    stack.push(searchNode.getRightChild());
            }
        }                             
        return null;
    }                

    // findLeastCommonAncestor
    // Pre-conditions:
    //    - node1 and node2 point to nodes in the phylogenetic tree
    // Post-conditions:
    //    - If node1 or node2 are null, return null
    //    - Else: returns the PhyloTreeNode of their common ancestor 
    //      with the largest depth
     private static PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode node1, PhyloTreeNode node2) {
                
        if(node1 == null || node2 == null)
            return null;
        else{
            int depth1, depth2;
            
            while(node1.getParent() != node2.getParent()){
                depth1 = nodeDepth(node1);
                depth2 = nodeDepth(node2);
                
                if(depth1> depth2)
                    node1 = node1.getParent();
                if(depth1 < depth2)
                    node2 =node2.getParent();
                if(depth1 == depth2){
                    node1 = node1.getParent();
                    node2 = node2.getParent();
                }                    
            }
        }
        return node1.getParent();
    }    
    
}//End PhyloTree class




