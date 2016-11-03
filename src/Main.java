import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class Main {
	//Genes_relation.data attributes.
	static int GENES_ID = 0,  GENES_ESS = 1, GENES_CLSS = 2, GENES_CPLX = 3, GENES_PHENO = 4,
			GENES_MOT = 5, GENES_CHROM = 6, GENES_FUNC = 7, GENES_CLASS_LOC = 8;
	//Interactions_relation.data attributes.
	static int INTE_ID1 = 0, INTE_ID2 = 1, INTE_TYPE = 2, INTE_CORR = 3;
	
	static double GENES_TUPLES = 4346;

	public static void main(String[] args) throws FileNotFoundException{

		String[][] genes_relation = genesRelation("Genes_relation.data",4346);
		String[][] genes_relation_test = genesRelation("Genes_relation.test", 1929);
		//String[][] interactions_relation = interactionsRelation("Interactions_relation.data");

		/*
		 * Attribute list to be used for generating the decision tree.
		 */
		HashMap<Integer, String> attribute_list = new HashMap<>();
		attribute_list.put(GENES_ESS, "essential");
		attribute_list.put(GENES_CLSS, "class");
		attribute_list.put(GENES_CPLX, "complex");
		attribute_list.put(GENES_PHENO, "phenotype");
		attribute_list.put(GENES_MOT, "motif");
		attribute_list.put(GENES_CHROM, "chromosome number");
		//attribute_list.put(GENES_FUNC, "function"); //not using function attribute!
		

		//Generating decision tree from training data set.
		Node<String[][]> tree = decisionTree(genes_relation, attribute_list);
		
		/*	Printing String representation of tree
		 *  Print format:
		 *  	if not a root
		 *  		parentLabel (if leaf prints "leaf") : attributeValue -> nodeLabel
		 *  	
		 *  	if root
		 *  		nodeLabel 
		 */
		
		System.out.println("START TREE PRINT");
		printTree(tree,0);
		System.out.println("END TREE PRINT");
		
		
		/*
		 * Determining the accuracy of the decision tree.
		 */
		double accuracy = 0;
		
		
		Scanner in = new Scanner(new File("keys.txt"));
		
		String keyid = "";
		String keyClass = "";
		
		HashMap<String, String> keys = new HashMap<>();
		while(in.hasNextLine()){
			String[] commaSplit = in.nextLine().split(",");
			
			keyid = commaSplit[0];
			keyClass = commaSplit[1];
			
			keys.put(keyid, keyClass);
		}

		in.close();
		
		int hits = 0;
		PrintWriter out = new PrintWriter("out.txt");
		for(String[] tuple : genes_relation_test){
			String id = tuple[GENES_ID];
			String clss = predictLocalization(tuple, tree, attribute_list);
			
			out.println(id + ", " + clss);
			
			
			/*
			 * Comparing the values.
			 * 
			 */
			if(keys.containsKey(id)){
				//System.out.println("key: " + keys.get(id) + " test: " +  clss);
				String keyValue = "";
				String testValue = "";
				
				keyValue = keys.get(id) + ".";//adding period to match file reading format
				testValue = clss.replaceAll("\\s+", "");//removing white space to match key reading format.
				
				if(keyValue.equals(testValue)){
					hits++;
				}
			}
		}
		out.close();
		accuracy = (double)hits/(double)keys.size();
		
		System.out.println("Accuracy: " + accuracy);



	}
	
	/*
	 * Prints the given tree with root node @param N.
	 */
	public static void printTree(Node<String[][]> N, int index){
		
		
		
		if(!N.isRoot() && !N.isLeaf())
			System.out.println(N.getParentLabel() + ":" + N.getLabelValue() + "->" + N.getLabel());
		else if(N.isLeaf())
			System.out.println(N.getParentLabel() + " leaf :" + N.getLabelValue() + "->" + N.getLabel());
		else
			System.out.println(N.getLabel());
		
		
		if(!N.isLeaf()){
			for(Node<String[][]> n: N.getChildren()){
				printTree(n, index++);
			}
			
		}
			
		
	}
	
	/*
	 * Predicts the localization attribute of the given
	 * test data set.
	 * 
	 * @param decisionTree: The tree root node that will be used
	 * to predict the localization value of the given
	 * test data set.
	 * 
	 */
	public static String predictLocalization(String[] tuple, Node<String[][]> N, HashMap<Integer, String> attrList){
		
		if(N.isLeaf()){
			return N.getLabel();
		}
		
		int curAttribute = -1;
		if(N.hasLabelIndex())
			curAttribute = N.getLabelIndex();
		else
			throw new IllegalStateException();
		
		
		Node<String[][]> curNode = null;
		
		
		for(Node<String[][]> child : N.getChildren()){
			if(child.getLabelValue().equals(tuple[curAttribute])){
				curNode = child;
				break;
			}
		}
		if(curNode == null){
			curNode = N.getChildren().get(0);
		}
		
		
		return predictLocalization(tuple, curNode, attrList);
		
		
	}

	//Stepping through the Inducing a decision tree algorithm as shown in 8.2 from book.
	public static Node<String[][]> decisionTree(String[][] D, HashMap<Integer, String> attribute_list){
		
		//Creating a new node N.
		Node<String[][]> N = new Node<String[][]>(D);
		
		
		/*
		 *  FIRST TERMINATING CONDITION
		 *  If tuples in root are all of the same class.
		 */
		String classLabel = sameClass(D);
		if(classLabel != null){
			
			N.setLabel(classLabel);
			return N;
		}

		/*
		 * SECOND TERMINATING CONDITION
		 * If attribute_list is empty
		 */
		
		if(attribute_list.isEmpty()){
			N.setLabel(getMajorityClass(D));
			//System.out.println("attr list empty class label:" + getMajorityClass(D));
			return N;
		}

		//Applying attribute selection method: information gain.
		HashMap.Entry<Integer, String> splitCriterion = info(N.getData(), attribute_list);
		N.setLabel(splitCriterion.getValue());
		N.setLabelIndex(splitCriterion.getKey());
		//System.out.println("Split label:" + splitCriterion.getValue());
		
		//for each outcome of the splitting_criterion
		//partition the tuples and grow subtrees for each partition.
		HashMap<String, Integer> splitValues = getValues(N.getData(), splitCriterion.getKey());
		HashMap<Integer, String> newAttList = attribute_list;
	
		newAttList.remove(splitCriterion.getKey());

		for (HashMap.Entry<String, Integer> e : splitValues.entrySet())
		{	
			String[][] subset = getSubset(D, splitCriterion.getKey(), e.getKey());
		
			/*
			 * THIRD TERMINATING CONDITION
			 * If subset is empty then attach a leaf node labeled with the majority class of D.
			 */
			if(subset.length == 0){
				
				//String[][] majorityClass = getSubset(D, GENES_CLASS_LOC, getMajorityClass(N.getData()));
				Node<String[][]> leaf = new Node<String[][]>(new String[0][0]);
				leaf.setLabel(getMajorityClass(D));
				//System.out.println("subset length 0:" + getMajorityClass(D));
				N.addChild(leaf);
				
				
			}
			else{
				
				N.addChild(decisionTree(subset, newAttList), e.getKey());
				
			}
			
			
		}
		
		return N;
		
	}
	
	/*
	 * Returns a subset of the data set.
	 * @param attr: The global static index that points to the attribute
	 * from where all values will be checked.
	 * @param value: The String value that will be checked for every
	 * tuple. If a tuple has this value at the given attribute then that
	 * tuple will be added to the sub set.
	 */
	public static String[][] getSubset(String[][] D, int attr, String value){
		
		ArrayList<String[]> subsetList = new ArrayList<>();
		
		for(int i = 0; i < D.length; i++){
			
			if(D[i][attr].equals(value)){
				subsetList.add(D[i]);
			}
		}
		
		String[][] subset = new String[subsetList.size()][];
		
		subset = subsetList.toArray(subset);
		
		
		return subset;
	}

	/*
	 * Computes the informational gain of the given data set.
	 * @param attribute_list: The input attribute list of the
	 * data set.
	 * Returns a HasMap.Entry that contains the attribute's static
	 * global index that points to that attribute in the data set.
	 * The global index is the key in the Hashmap.Entry return
	 * object, and it also returns a String representation of the
	 * name of the attribute as the value of the Hashmap.Entry.
	 */
	public static HashMap.Entry<Integer, String> info(String[][] D, HashMap<Integer, String> attribute_list){

		double info_D = 0;
		double P_i = 0;
		

		HashMap<String, Integer> classValues = getValues(D, GENES_CLASS_LOC);
		
		
		for (HashMap.Entry<String, Integer> e : classValues.entrySet())
		{	    
			P_i = e.getValue()/GENES_TUPLES;
			info_D += P_i * Math.log10(P_i)/Math.log10(2);
			
		}		
		info_D = info_D * -1;
		
		HashMap.Entry<Integer, String> splitAttribute = null;
		double maxInfo = Integer.MAX_VALUE * -1;
		
		for (HashMap.Entry<Integer, String> e : attribute_list.entrySet())
		{	
			splitAttribute = e;
			HashMap<String, Integer> curAttribute = getValues(D, e.getKey());
			double D_i = 0;
			double d = curAttribute.size();
			double info_Da = 0;
			
			for (HashMap.Entry<String, Integer> a : curAttribute.entrySet())
			{	    				 				
				D_i = a.getValue();
				info_Da += (D_i/d) * info_D;
				
			}
			
			if((info_D - info_Da) > maxInfo){
				maxInfo = info_D - info_Da;
				splitAttribute = e;
			}
			
		}
		
		return splitAttribute;
	}
	
	/*
	 * Returns the values in the data set of the attribute given.
	 * Returns a Hashmap with the value at the given attribute as
	 * the key, and the amount of times that value is present as the
	 * Hashmap value.
	 */
	public static HashMap<String, Integer> getValues(String[][] D, int v){
		
		HashMap<String, Integer> values = new HashMap<>();

		for(int tuple = 0; tuple < D.length; tuple++){

			String value = D[tuple][v];

			if(values.containsKey(value)){
				values.put(value, values.get(value) + 1);
			}
			else{
				values.put(value, 1);
			}



		}
		
		return values;
	}

	/*
	 * Returns the majority class in the data set, as a String.
	 */
	public static String getMajorityClass(String[][] D){

		HashMap<String, Integer> classValues = new HashMap<>();

		for(int tuple = 0; tuple < D.length; tuple++){

			String classValue = D[tuple][GENES_CLASS_LOC];

			if(classValues.containsKey(classValue)){
				classValues.put(classValue, classValues.get(classValue) + 1);
			}
			else{
				classValues.put(classValue, 1);
			}



		}

		String majorityClass = "";
		int max = 0;

		for (HashMap.Entry<String, Integer> e : classValues.entrySet())
		{	    
			int curCount = e.getValue();

			if(curCount > max){
				majorityClass = e.getKey();
			}

		}

		return majorityClass;

	}

	/*
	 * Searches the input dataset's class attribute and
	 * determines if all tuples have the same class.
	 * If all tuples have the same class it returns the class name.
	 * Otherwise it returns null. 
	 *
	 */
	public static String sameClass(String[][] D){

		String curClassValue = D[0][GENES_CLASS_LOC];

		for(int tuple = 1; tuple < D.length; tuple++){

			String nextClassValue = D[tuple][GENES_CLASS_LOC];
			if(!curClassValue.equals(nextClassValue))
				return null;

			curClassValue = nextClassValue;



		}


		return curClassValue;
	}


	public static String[][] genesRelation(String fileName, int numTuples) throws FileNotFoundException{

		Scanner in = new Scanner(new File(fileName));
		String[][] returnSet = new String[numTuples][20];


		int index = 0;
		while(in.hasNextLine()){
			String line = in.nextLine();
			String[] quotationSplit = line.split("\"");
			for(int i = 0; i < quotationSplit.length ; i++){
				if(i != 0 && i%2 != 0){
					quotationSplit[i] = quotationSplit[i].replaceAll(",", " ");
					
				}
			}
			line = "";
			for(String s : quotationSplit){
				line += s;
			}
			
			String[] commaSplit = line.split(",");
			
			
			
			returnSet[index][0] = commaSplit[0];
			returnSet[index][1] = commaSplit[1];
			returnSet[index][2] = commaSplit[2];
			returnSet[index][3] = commaSplit[3];
			returnSet[index][4] = commaSplit[4];
			returnSet[index][5] = commaSplit[5];
			returnSet[index][6] = commaSplit[6];
			returnSet[index][8] = commaSplit[commaSplit.length - 1];


			++index;



		}
		in.close();


		return returnSet;



	}

	/*
	 * Must edit to combine this function with the other file read function.
	 */
	public static String[][] interactionsRelation(String fileName) throws FileNotFoundException{

		Scanner in = new Scanner(new File(fileName));
		String[][] returnSet = new String[910][4];


		int index = 0;
		while(in.hasNextLine()){
			String line = in.nextLine();
			String[] commaSplit = line.split(",");

			returnSet[index][0] = commaSplit[0];
			returnSet[index][1] = commaSplit[1];
			returnSet[index][2] = commaSplit[2];
			returnSet[index][3] = commaSplit[3];


			++index;



		}
		in.close();


		return returnSet;



	}

	/*
	 * Prints a 2D array
	 * @param a: input 2D array to print.
	 */
	public static void printArray(String[][] a){
		for(int i = 0; i < a.length; i++){
			for(int j = 0; j < a[i].length; j++){
				System.out.print(a[i][j] + ",");
			}
			System.out.println();
		}
	}

}
