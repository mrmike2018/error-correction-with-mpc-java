import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Random;
import java.util.Vector;

public class ErrorCorrectionMPC_Ver02 {

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//p: prime number for Z_p
	//static BigInteger p = new BigInteger("290600942603312283586666508900447746303"); // prime number for Z_p
	//static BigInteger p = new BigInteger("4225480274875879477"); // 64 bit
	//static BigInteger p = new BigInteger("3384451423"); // 32 bit
	//static BigInteger p = new BigInteger("37997"); // 16 bit
	static BigInteger p = new BigInteger("19"); // 8 bit
	
	static int threshold = 6;
	static int numOfErrors;// = 3;
	static int numOfPlayers = 2 * numOfErrors + threshold;
	static int maximumNumOfPlayers = 20;
	static int publicSubMatrixDim = numOfErrors + threshold;
	static int primePsizeInBits = 8;
	
	// array for keeping shares of players:
	static BigInteger[] evaluations = new BigInteger[maximumNumOfPlayers];
	static BigInteger[] evaluationsWithoutErrors = new BigInteger[maximumNumOfPlayers];
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args) throws Exception {
		
		//runMPCErrorCorrection();
		runMPCErrorCorrectionWithDifferentErrors();
		

		System.out.println("Done!");
	}	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void runMPCErrorCorrectionWithDifferentErrors() throws Exception{
		
		// creating random prime for p:
		Random rnd = new Random();
		p = new BigInteger(primePsizeInBits, 10, rnd);
		
		System.out.println("prime p: " + p + "\nprimePsizeInBits: " + primePsizeInBits + " bits");
		System.out.println("is p prime: " + p.isProbablePrime(10));
		//System.out.println("numOfErrors:" + numOfErrors +"      threshold:" + threshold + "       numOfPlayers:" + numOfPlayers);
		System.out.println("threshold:" + threshold + "       numOfPlayers:" + numOfPlayers);
		
		
		//System.out.println("publicSubMatrixDim: " + publicSubMatrixDim);
		
		generateShares();
		//generateSharesNonRandom();
		//createSomeErrors();
		
		Boolean useAllPlayers = false;
		// the following part is for getting execution time for different errors etc.:
		for(int i = 1; i <= 4; i++) {
			numOfErrors = i;
			numOfPlayers = 2 * numOfErrors + threshold;
			publicSubMatrixDim = numOfErrors + threshold;
			createSomeErrorsAtRandomLocations(i);
			System.out.println("------------------------------measuring time for finding the location of " + i + " error(s), start:");
			int numOfPlayersRequired = 2 * i + threshold;
			
			// start counting execution time:
			long lStartTime = System.nanoTime();
			findErrorLocations(i, numOfPlayersRequired, useAllPlayers);
			// get time after finishing the task:
			long lEndTime = System.nanoTime();
			long output = lEndTime - lStartTime;
	        System.out.println("Elapsed time in milliseconds: " + output / 1000000);
	        System.out.println("Elapsed time in seconds: " + output / 1000000 / 1000);
	        System.out.println("Elapsed time in minutes: " + output / 1000000 / 1000 / 60);
			
			System.out.println("");
			System.out.println("------------------------------measuring time for finding the location of  " + i + " error(s), ended.");
			
	        for(int j = 0; j < maximumNumOfPlayers; j++) {
	        	//evaluationBeforeError[i] = evaluations[i];
	        	evaluations[j] = evaluationsWithoutErrors[j];        	
	        }
	        System.out.println("\n\n\n----------------------------------------------------------------------");
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void runMPCErrorCorrection() throws Exception{
		System.out.println("prime p: " + p + "      numOfErrors:" + numOfErrors + "      threshold:" + threshold + "       numOfPlayers:" + numOfPlayers);
		System.out.println("publicSubMatrixDim: " + publicSubMatrixDim);
		
		generateShares();
		//generateSharesNonRandom();
		createSomeErrors();
		
		Boolean useAllPlayers = false;
		for(int err = numOfErrors; err <= numOfErrors; err++) {
			System.out.println("------------------------------catching " + err + " errors:");
			int numOfPlayersRequired = 2 * err + threshold;
			findErrorLocations(err, numOfPlayersRequired, useAllPlayers);
			System.out.println("");
			System.out.println("------------------------------end catching " + err + " errors:");
		}
	}	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void createSomeErrorsAtRandomLocations(int numOfErrsToBeCreated) {
		
		BigInteger[] evaluationBeforeError = new BigInteger[maximumNumOfPlayers];
		
        for(int i = 0; i < maximumNumOfPlayers; i++) {
        	//evaluationBeforeError[i] = evaluations[i];
        	evaluationBeforeError[i] = evaluationsWithoutErrors[i];        	
        }
        numOfPlayers = 2 * numOfErrsToBeCreated + threshold;
        // the following array is used to for saving the location of errors to be created
        final int[] ints = new Random().ints(1, numOfPlayers).distinct().limit(numOfErrsToBeCreated).toArray();
        //System.out.println("ints array: ");
        //System.out.println(Arrays.toString(ints));
        
        // creating some errors in the shares:
        for(int i = 0; i < ints.length; i++) {
        	int locationOfError = ints[i];
        	evaluations[locationOfError - 1] = evaluations[locationOfError - 1].add(BigInteger.valueOf(1)).mod(p);
        }
		
        // print shares of players after injecting errors:
        System.out.println("shares of players after injecting errors:");
        for(int i = 0; i < numOfPlayers; i++) {
        	System.out.print(evaluations[i] + "  ");
        }
        System.out.println("");
        
        System.out.println("------------------------------error locations:");
        for(int i = 0; i < maximumNumOfPlayers; i++) {
        	int res = evaluationBeforeError[i].compareTo(evaluations[i]);
        	int playerID = i + 1;
        	if (res == 1 || res == -1)
        		System.out.println("location " + i + "(i.e. player " + playerID + "):   " + evaluationBeforeError[i] + " -->  " + evaluations[i]);
        }
        System.out.println("");
        
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void createSomeErrors() {
		
		BigInteger[] evaluationBeforeError = new BigInteger[maximumNumOfPlayers];
		
        for(int i = 0; i < maximumNumOfPlayers; i++) {
        		evaluationBeforeError[i] = evaluations[i];
        }
		
		// creating some errors in the shares:
        //evaluations[0] = evaluations[0].add(BigInteger.valueOf(1)).mod(p);
        //evaluations[1] = evaluations[1].add(BigInteger.valueOf(1)).mod(p);
        evaluations[2] = evaluations[2].add(BigInteger.valueOf(1)).mod(p);
        evaluations[3] = evaluations[3].add(BigInteger.valueOf(1)).mod(p);
        //evaluations[4] = evaluations[4].add(BigInteger.valueOf(1)).mod(p);
        evaluations[5] = evaluations[5].add(BigInteger.valueOf(1)).mod(p);
		
        // print shares of players after injecting errors:
        System.out.println("shares of players after injecting errors:");
        for(int i = 0; i < numOfPlayers; i++) {
        	System.out.print(evaluations[i] + "  ");
        }
        System.out.println("");
        
        System.out.println("------------------------------error locations:");
        for(int i = 0; i < maximumNumOfPlayers; i++) {
        	int res = evaluationBeforeError[i].compareTo(evaluations[i]);
        	int playerID = i + 1;
        	if (res == 1 || res == -1)
        		System.out.println("location " + i + "(i.e. player " + playerID + "):   " + evaluationBeforeError[i] + " -->  " + evaluations[i]);
        }
        System.out.println("");
        
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void generateShares(){
		
		Vector <BigInteger> randomCoefficients = new Vector <BigInteger>(threshold);
		
		for(int i = 0; i < threshold; i++) {
			randomCoefficients.add(i, randomBigInteger(p));
			//System.out.println("randomCoefficients[" + i+ "]: " + randomCoefficients.get(i) + "  ");
		}
		
		BigIntegerPoly Px = new BigIntegerPoly(randomCoefficients);
		
		System.out.println("Px = " + Px.toPString());
        
        // creating players' shares:
        for(int i = 0; i < maximumNumOfPlayers; i++) {
        	evaluations[i] = Px.valueOf(i+1).mod(p);
        	evaluationsWithoutErrors[i] = Px.valueOf(i+1).mod(p);;
        	//System.out.print("evaluations[" + i + "]: " + evaluations[i] + "    ");
        }
        
        // print shares of players:
        System.out.println("shares of players:");
        for(int i = 0; i < numOfPlayers; i++) {
        	System.out.print(evaluations[i] + "  ");
        }
        System.out.println("");
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static void generateSharesNonRandom(){
//		//141 +14*x +18*x^2
//    	Polynomial p0 = new Polynomial(141 % p.intValue(), 0);
//        Polynomial p1 = new Polynomial(14 % p.intValue(), 1);
//        Polynomial p2 = new Polynomial(18 % p.intValue(), 2);
		
		//1 + 3*x + 2*x^2
    	Polynomial p0 = new Polynomial(1 % p.intValue(), 0);
        Polynomial p1 = new Polynomial(3 % p.intValue(), 1);
        Polynomial p2 = new Polynomial(2 % p.intValue(), 2);
//        Polynomial p3 = new Polynomial(3, 3);
//        Polynomial p4 = new Polynomial(3, 4);
//        Polynomial Px = p0.plus(p1).plus(p2).plus(p3).plus(p4);
        Polynomial Px = p0.plus(p1).plus(p2);
        
        System.out.println("Px = " + Px);
        
        // creating players' shares:
        int PxEval;
        for(int i = 0; i < maximumNumOfPlayers; i++) {
        	PxEval = Px.evaluate(i+1) % p.intValue();
        	evaluations[i] = BigInteger.valueOf(PxEval);
        	//System.out.print("evaluations[" + i + "]: " + evaluations[i] + "    ");
        }
        
        // print shares of players:
        System.out.println("shares of players:");
        for(int i = 0; i < numOfPlayers; i++) {
        	System.out.print(evaluations[i] + "  ");
        }
        System.out.println("");
	}
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public static Vector<BigInteger> findErrorLocations(int possibleNumOfErrors, int numOfPlayersRequired, Boolean allPlayersOrNot) throws Exception {
    	System.out.println("possibleNumOfErrors: " + possibleNumOfErrors + "\nthreshold: " + threshold);
    	System.out.println("numOfPlayersRequired: " + numOfPlayersRequired + "\n");

    	// create the system of equations with maximum possible dimensions
    	if(allPlayersOrNot == true) {
    		numOfPlayersRequired = numOfPlayers;
    		possibleNumOfErrors = numOfErrors;
    	}
    	
    	// matrix of coefficient in the system of equations (Ax = c):
        BigInteger[][] matrixA = new BigInteger[numOfPlayersRequired][numOfPlayersRequired];
        BigInteger[][] matrixACopy = new BigInteger[numOfPlayersRequired][numOfPlayersRequired];
        
        // creating matrix A dynamically:
        // first part of matrix A:
        for (int i = 0; i < numOfPlayersRequired; i++) {
        	BigInteger bigintI = BigInteger.valueOf(i + 1);
           	for (int j = 0; j < possibleNumOfErrors + threshold; j++) {
            	BigInteger bigintJ = BigInteger.valueOf(j);
            	matrixA[i][j] = bigintI.modPow(bigintJ, p);
            	matrixACopy[i][j] = bigintI.modPow(bigintJ, p);
            }
        }
        // second part of matrix A:
        BigInteger negativeOne = BigInteger.valueOf(-1).mod(p);
        for (int k = 0; k < possibleNumOfErrors; k++) {
            for (int i = 0; i < numOfPlayersRequired; i++) {
            	BigInteger bigintI = BigInteger.valueOf(i + 1);
            	BigInteger iToThePowerOfEi = bigintI.modPow(BigInteger.valueOf(k), p);
            	BigInteger ijElementOfMatrixA = negativeOne.multiply(iToThePowerOfEi).multiply(evaluations[i]).mod(p);
            	
            	matrixA[i][possibleNumOfErrors + threshold + k] = ijElementOfMatrixA;
            	matrixACopy[i][possibleNumOfErrors + threshold + k] = ijElementOfMatrixA;
            }
        }
        
        // vector of solutions, c:      
        //System.out.println("\ncVec:");
        BigInteger cVec[] = new BigInteger[numOfPlayersRequired];
        for (int i = 0; i < numOfPlayersRequired; i++) {
        	BigInteger bigintI = BigInteger.valueOf(i + 1);
        	BigInteger iToThePowerOfE = bigintI.modPow(BigInteger.valueOf(possibleNumOfErrors), p);
        	cVec[i] = iToThePowerOfE.multiply(evaluations[i]).mod(p);
        	//System.out.print(cVec[i] + "  ");
        }
        //System.out.println("");
        
        ModMatrix matrixAObject = new ModMatrix(matrixA);
        //ModMatrix matrixAObject2 = new ModMatrix(matrixA);

        // printing matrix A:
        //System.out.println("\nmatrix A:");
        //ModMatrix.printModMatrix(matrixAObject);
        //ModMatrix.saveToFile("matrix A:");
        //ModMatrix.saveToFileModMatrix(matrixAObject);
        
        System.out.println("-------------calculating matrixAdet, please wait...");
        BigInteger matrixAdet = matrixAObject.determinant(matrixAObject);

        //System.out.println("***************matrixAdet: " + matrixAdet + "    mod " + p.toString() + "************");
        
        
        //ModMatrix.saveToFile("matrixAObject2 data array:");
        //ModMatrix.saveToFileBigIntMatrix(matrixAObject2.getData());
        
        //ModMatrix.saveToFile("\n\n\n\n\n-------------------------------------finding bi'1-------------------------------------");
        //BigInteger b0 = solutionsVec[numOfPlayers -1];
        BigInteger[] bi = new BigInteger[possibleNumOfErrors];
        Vector<BigInteger> ExCoefficients = new Vector<BigInteger>(possibleNumOfErrors);
        //for (int i = 0; i < possibleNumOfErrors; i++) {
        //for (int i = possibleNumOfErrors - 1; i > -1; i--) {
        for (int i = 0; i < possibleNumOfErrors; i++) {
        	
        	//Cramers' rule:
        	BigInteger[][] matrixAiThColumnReplacedWithSolutions = new BigInteger[numOfPlayersRequired][numOfPlayersRequired];
        	
        	for (int j = 0; j < numOfPlayersRequired; j++) {
        		for (int k = 0; k < numOfPlayersRequired; k++) {
        			matrixAiThColumnReplacedWithSolutions[j][k] = matrixA[j][k];
        		}
        	}
        	
        	int temp = threshold + possibleNumOfErrors + i;
        	//System.out.println("threshold + possibleNumOfErrors + i:" + temp);
        	for (int j = 0; j < numOfPlayersRequired; j++)
        		matrixAiThColumnReplacedWithSolutions[j][temp] = cVec[j];
        	
        	ModMatrix matrixAiThColumnReplacedWithSolutionsObject = new ModMatrix(matrixAiThColumnReplacedWithSolutions);
        	
        	BigInteger matrixAiThColumnReplacedWithSolutionsDet = matrixAiThColumnReplacedWithSolutionsObject.determinant(matrixAiThColumnReplacedWithSolutionsObject);
        	
            //System.out.println("-------------matrixAiThColumnReplacedWithSolutionsDet: " + matrixAiThColumnReplacedWithSolutionsDet + "    mod " + p.toString() + "************");
        	//ModMatrix.saveToFile("-------------matrixAiThColumnReplacedWithSolutionsDet: " + matrixAiThColumnReplacedWithSolutionsDet + "    mod " + p.toString() + "************");
        	
        	//bi[i] = matrixAiThColumnReplacedWithSolutionsDet.multiply(matrixAdet.modInverse(p)).mod(p);
        	bi[i] = matrixAiThColumnReplacedWithSolutionsDet.multiply(matrixAdet.modInverse(p)).mod(p);
        	//ExCoefficients.add(i, bi[i]);
        	//int biArrayIndex = possibleNumOfErrors - 1 - i;
        	//System.out.println("-------------i: " + i + "     possibleNumOfErrors: " + possibleNumOfErrors);
        	//System.out.println(">>>>>>>>>>>>>>>>>>>>>bi[" + biArrayIndex + "]: " + bi[i]);
        	
        	ExCoefficients.add(i, bi[i]);
        }
        
        
        ExCoefficients.add(possibleNumOfErrors, BigInteger.ONE);
        
//        System.out.println("");
//        for (int i = 0; i < possibleNumOfErrors; i++) {
//        	System.out.print("b" + i + ": " + bi[i] + "    ");
//        }
        
        // constructing the E(x) polynomial:
        BigIntegerPoly Ex = new BigIntegerPoly(ExCoefficients);
        
        
        
        System.out.println("\nEx = " + Ex.toPString());
        
        // finding roots of Ex, errors locations:
        Vector<BigInteger> rootsOfEx = new Vector<BigInteger>(possibleNumOfErrors);
        
        //rootsOfEx = Ex.iroots();
        rootsOfEx = findBigIntegerPolyRoots(Ex, BigInteger.ZERO, p);
        
        if(rootsOfEx.size() > 0) {
	        for (int i = 0; i < possibleNumOfErrors; i++) {
	        	BigInteger temp = rootsOfEx.get(i);
	        	rootsOfEx.set(i, temp.mod(p));
	        }
	        
	        System.out.println("roots of Ex (location of errors):");
	        for (int i = 0; i < possibleNumOfErrors; i++) {
	        	System.out.print(rootsOfEx.get(i) + "\n");
	        	
	        }
        }
        else
        	System.out.println("Ex has no roots in Z" + p);
		return rootsOfEx;
    }
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public static Vector<BigInteger> findErrorLocations00(int possibleNumOfErrors, int numOfPlayersRequired, Boolean allPlayersOrNot) throws Exception {
    	System.out.print("possibleNumOfErrors: " + possibleNumOfErrors + "         numOfPlayersRequired: " + numOfPlayersRequired);

    	// create the system of equations with maximum possible dimensions
    	if(allPlayersOrNot == true) {
    		numOfPlayersRequired = numOfPlayers;
    		possibleNumOfErrors = numOfErrors;
    	}
    	
    	// matrix of coefficient in the system of equations (Ax = c):
        BigInteger[][] matrixA = new BigInteger[numOfPlayersRequired][numOfPlayersRequired];
        
        // creating matrix A dynamically:
        // first part of matrix A:
        for (int i = 0; i < numOfPlayersRequired; i++) {
        	BigInteger bigintI = BigInteger.valueOf(i + 1);
           	for (int j = 0; j < possibleNumOfErrors + threshold; j++) {
            	BigInteger bigintJ = BigInteger.valueOf(j);
            	matrixA[i][j] = bigintI.modPow(bigintJ, p);
            }
        }
        // second part of matrix A:
        BigInteger negativeOne = BigInteger.valueOf(-1).mod(p);
        for (int k = 0; k < possibleNumOfErrors; k++) {
            for (int i = 0; i < numOfPlayersRequired; i++) {
            	BigInteger bigintI = BigInteger.valueOf(i + 1);
            	BigInteger iToThePowerOfEi = bigintI.modPow(BigInteger.valueOf(k), p);
            	BigInteger ijElementOfMatrixA = negativeOne.multiply(iToThePowerOfEi).multiply(evaluations[i]).mod(p);
            	
            	matrixA[i][possibleNumOfErrors + threshold + k] = ijElementOfMatrixA;
            }
        }
        
        // vector of solutions, c:      
        System.out.println("\ncVec:");
        BigInteger cVec[] = new BigInteger[numOfPlayersRequired];
        for (int i = 0; i < numOfPlayersRequired; i++) {
        	BigInteger bigintI = BigInteger.valueOf(i + 1);
        	BigInteger iToThePowerOfE = bigintI.modPow(BigInteger.valueOf(possibleNumOfErrors), p);
        	cVec[i] = iToThePowerOfE.multiply(evaluations[i]).mod(p);
        	System.out.print(cVec[i] + "  ");
        }
        
        ModMatrix matrixAObject= new ModMatrix(matrixA);
        
        // printing matrix A:
        System.out.println("\nmatrix A:");
        matrixAObject.printModMatrix(matrixAObject);
        
        ModMatrix myModMatrix= new ModMatrix(matrixA);
        
        // checking if matrix is invertible, if determinant == 0:
        //System.out.println("matrixAdet: ");
        
        //matrixAObject.saveToFile("----------------matrix:");
        //matrixAObject.saveToFileModMatrix(matrixAObject);
        //matrixAObject.saveToFile("----------------");
        
        
        BigInteger matrixAdet = myModMatrix.determinant(matrixAObject).mod(p);
        System.out.println("***************matrixAdet: " + matrixAdet + "    mod " + p.toString() + "************");
//        if (matrixAdet.equals(BigInteger.valueOf(0))) {
//        	System.out.println("matrix A is not invertible.");
//        	Vector<BigInteger> rootsOfEx = new Vector<BigInteger>();
//        	rootsOfEx = null;
//        	return rootsOfEx;
//        }
        
        
        System.out.println("------------------------------calculating matrixA inverse:");
        // calculating matrix A inverse:
        ModMatrix matrixAinv = myModMatrix.inverse(myModMatrix);
        
        // printing matrix A inverse:
        System.out.println("------------------------------matrix A inv:");
        matrixAinv.printModMatrix(matrixAinv);
        
        BigInteger[][] matrixAinvBigInt = matrixAinv.getData();
        
//        System.out.println("\nmatrixAinvBigInt:");
//        matrixAObject.printMatrixBigInt(matrixAinvBigInt);
        
        // calculating solutions of the system of equations, that is x in Ax = c:
        BigInteger[] solutionsVec = matrixAObject.multiplyMatrixByVectorBigInt(matrixAinvBigInt, cVec);
        
        // printing solutions x:
        System.out.println("solutions of the sys of equation:");
        for (int i = 0; i < solutionsVec.length; i++) {
        	System.out.print(solutionsVec[i]+" ");
        }
        
        //BigInteger b0 = solutionsVec[numOfPlayers -1];
        BigInteger[] bi = new BigInteger[possibleNumOfErrors];
        Vector<BigInteger> ExCoefficients = new Vector<BigInteger>(possibleNumOfErrors);
        for (int i = 0; i < possibleNumOfErrors; i++) {
        	bi[i] = solutionsVec[threshold + possibleNumOfErrors + i];
        	ExCoefficients.add(i, solutionsVec[threshold + possibleNumOfErrors + i]);
        }
        
        
        ExCoefficients.add(possibleNumOfErrors, BigInteger.ONE);
        
        System.out.println("");
        for (int i = 0; i < possibleNumOfErrors; i++) {
        	System.out.print("b" + i + ": " + bi[i] + "    ");
        }
        
        // constructing the E(x) polynomial:
        BigIntegerPoly Ex = new BigIntegerPoly(ExCoefficients);
        
        
        
        System.out.println("\nEx = " + Ex.toPString());
        
        // finding roots of Ex, errors locations:
        Vector<BigInteger> rootsOfEx = new Vector<BigInteger>(possibleNumOfErrors);
        
        //rootsOfEx = Ex.iroots();
        rootsOfEx = findBigIntegerPolyRoots(Ex, BigInteger.ZERO, p);
        
        if(rootsOfEx.size() > 0) {
	        for (int i = 0; i < possibleNumOfErrors; i++) {
	        	BigInteger temp = rootsOfEx.get(i);
	        	rootsOfEx.set(i, temp.mod(p));
	        }
	        
	        System.out.println("roots of Ex (location of errors):");
	        for (int i = 0; i < possibleNumOfErrors; i++) {
	        	System.out.print(rootsOfEx.get(i) + "\n");
	        	
	        }
        }
        else
        	System.out.println("Ex has no roots in Z" + p);
		return rootsOfEx;
    }
	
    public static Vector<BigInteger> findBigIntegerPolyRoots(BigIntegerPoly poly, BigInteger from, BigInteger end) {
    	
    	Vector<BigInteger> rootsOfPoly = new Vector<BigInteger>();
    	
    	if(poly.degree() == 1) {
    		BigInteger root = poly.at(0).negate().mod(p);
    		rootsOfPoly.add(0, root);
    		return rootsOfPoly;
    	}
    	
    	int k = 0;
    	for (BigInteger i = from; i.compareTo(end) < 0; i = i.add(BigInteger.ONE)){
        	if(poly.valueOf(i).mod(p).equals(BigInteger.ZERO)) {
        		rootsOfPoly.add(k, i);
        		//System.out.println("root" + k + ": " + i);
        		k++;
        	}
        	
        	if(poly.degree() == rootsOfPoly.size())
        		return rootsOfPoly;
        }
        
        return rootsOfPoly;
    }
    
    public static ArrayList findPolyRoots(Polynomial poly, int from, int end) {
    	
    	ArrayList<Object> rootsOfPoly = new ArrayList<>();
        for(int i = from; i <= end; i++) {
        	if(poly.evaluate(i) % p.intValue() == 0) {
        		rootsOfPoly.add(i);
        		//System.out.println(i + " ");
        	}
        }
        
        return rootsOfPoly;
    }
    
    public static void printArrayList(ArrayList arrayList) {
    	for(int i = 0; i < arrayList.size(); i++) {
    		System.out.print(arrayList.get(i) + "   ");
    	}
    	System.out.println("");
    }
    
	public static BigInteger randomBigInteger(BigInteger n) {
		// Generate random integers in range 0 to p for coefficient:
        Random rnd = new Random();
        int maxNumBitLength = n.bitLength();
        BigInteger aRandomBigInt;
        do {
            aRandomBigInt = new BigInteger(maxNumBitLength, rnd);
            // compare random number lessthan ginven number
        } while (aRandomBigInt.compareTo(n) > 0); 
        return aRandomBigInt;
    }
	
	
	private static BigInteger[][] intMatrixToBigIntMatrix(int[][] intA) {
		int n = intA.length;
		BigInteger[][] bigIntA = new BigInteger[n][n];
        for (int i = 0; i < n; i++) {
        	for (int j = 0; j < n; j++) {
        		bigIntA[i][j] = BigInteger.valueOf(intA[i][j]);
        	}
        	System.out.println();
        }
		
		return bigIntA;
	}
}
