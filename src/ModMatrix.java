//package determinantUsingLUDecomposition;

import static java.util.stream.IntStream.range;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.*;
import java.util.Vector;

public class ModMatrix {

	//static BigInteger p = new BigInteger("11"); // 8 bit
	static BigInteger p = ErrorCorrectionMPC_Ver02.p;
	
	//static int thresholdInModMat = ErrorCorrectionMPC_Ver02.threshold;
	//static int numOfErrorsInModMat = ErrorCorrectionMPC_Ver02.numOfErrors;
	//static int numOfPlayersInModMat = ErrorCorrectionMPC_Ver02.numOfPlayers;
	//static int publicSubMatrixDim = numOfErrorsInModMat + thresholdInModMat;
	//static int publicSubMatrixDim = 2;
	static int PmatrixDet = 1;
	//BigInteger det = BigInteger.ZERO;
	
    private int nrows;
    private int ncols;
    private BigInteger[][] data;
    private final BigInteger mod= p;

    public ModMatrix(BigInteger[][] dat) {
        this.data = dat;
        this.nrows = dat.length;
        this.ncols = dat[0].length;
    }

    public ModMatrix(int nrow, int ncol) {
        this.nrows = nrow;
        this.ncols = ncol;
        data = new BigInteger[nrow][ncol];
    }

    public int getNrows() {
        return nrows;
    }

    public void setNrows(int nrows) {
        this.nrows = nrows;
    }

    public int getNcols() {
        return ncols;
    }

    public void setNcols(int ncols) {
        this.ncols = ncols;
    }

    public BigInteger[][] getData() {
        return data;
    }

    public void setData(BigInteger[][] data) {
        this.data = data;
    }

    public BigInteger getValueAt(int i, int j) {
        return data[i][j];
    }

    public void setValueAt(int i, int j, BigInteger value) {
        data[i][j] = value;
    }

    public int size() {
        return ncols;
    }

    // Take the transpose of the Matrix..
    public static ModMatrix transpose(ModMatrix matrix) {
        ModMatrix transposedMatrix = new ModMatrix(matrix.getNcols(), matrix.getNrows());
        for (int i = 0; i < matrix.getNrows(); i++) {
            for (int j = 0; j < matrix.getNcols(); j++) {
                transposedMatrix.setValueAt(j, i, matrix.getValueAt(i, j));
            }
        }
        return transposedMatrix;
    }

    public BigInteger determinant(ModMatrix matrix) throws Exception {
    	//det = BigInteger.ZERO;
    	int lastColumnIndex = matrix.getNcols() - 1;
    	
    	//saveToFile("--------------------------------------------------------------------calculating determinant of: ");
    	//String tempStr = "";
    	//for (int i1 = 0; i1 < matrix.nrows; i1++)
    		//tempStr += "----------";
    	//saveToFile(tempStr + "calculating determinant of: ");
    	//saveToFileModMatrix(matrix);
    	
    	//saveToFile("-------------------------------------lastColumnIndex: " + lastColumnIndex);
    	//saveToFile("-------------------------------------determinant func started------------------------------------det: " + det);
//        if (matrix.size() == 1) {
//            return matrix.getValueAt(0, 0).mod(p);
//        }
//        if (matrix.size() == 2) {        	
//        	return matrix.getValueAt(0, 0).multiply(matrix.getValueAt(1, 1)).subtract((matrix.getValueAt(0, 1).multiply(matrix.getValueAt(1, 0)))).mod(p);
//        }
    	int publicSubMatrixDim = ErrorCorrectionMPC_Ver02.publicSubMatrixDim;
    	saveToFile("-------------------------------------publicSubMatrixDim:" + publicSubMatrixDim);
        if (matrix.size() == publicSubMatrixDim) {
        	//saveToFile("-------------------------------------matrix " + publicSubMatrixDim + " by " + publicSubMatrixDim + ": ");
        	//saveToFileModMatrix(matrix);
        	BigInteger temp = calculateDetUsingLUdecom(matrix.getData());
        	//saveToFile("-------------------------------------its det:" + temp + "    and det: " + det);
        	//saveToFile("-------------------------------------temp det:" + temp);
        	return temp;
        }
        
        //int lastColumnIndex = matrix.getNcols() - 1;
        //det = BigInteger.ZERO;
        BigInteger sum = new BigInteger("0");
        //saveToFile("**************************************det before for loop: " + sum);
        for (int i = 0; i < matrix.getNcols(); i++) {
//        	BigInteger subMatrixDet = determinant(createSubMatrix(matrix, i, lastColumnIndex));
////        	det = det.add(changeSign(i).multiply(matrix.getValueAt(i, lastColumnIndex).multiply(subMatrixDet))).mod(p);        	
//        	det = det.add(changeSign(i + lastColumnIndex).multiply(matrix.getValueAt(i, lastColumnIndex).multiply(subMatrixDet))).mod(p);
//        	saveToFile("********************matrix.getValueAt(" + i + ","  + lastColumnIndex + "): " + matrix.getValueAt(i, lastColumnIndex));
//        	saveToFile("********************changeSign(" + i + "+"  + lastColumnIndex + "): " + changeSign(i + lastColumnIndex));
//        	saveToFile("********************det inside for loop: " + det);
        	sum = sum.add(changeSign(i + lastColumnIndex).multiply(matrix.getValueAt(i, lastColumnIndex).multiply(determinant(createSubMatrix(matrix, i, lastColumnIndex)))));
        	sum = sum.mod(p);
        }
        //saveToFile("**************************************det after for loop: " + sum);
        //this.det = det;
        //saveToFile("-------------------------------------determinant func ended---------------det: " + det);
        //saveToFile(tempStr + "determinant func ended---------------det: " + sum);
        return sum;
    }

    private static BigInteger changeSign(int i) {
        if (i % 2 == 0) {
            return new BigInteger("1"); 
        } else {
            return new BigInteger("-1");
        }
    }

    public static ModMatrix createSubMatrix(ModMatrix matrix, int excluding_row, int excluding_col) {
        ModMatrix mat = new ModMatrix(matrix.getNrows() - 1, matrix.getNcols() - 1);
        int r = -1;
        for (int i = 0; i < matrix.getNrows(); i++) {
            if (i == excluding_row) {
                continue;
            }
            r++;
            int c = -1;
            for (int j = 0; j < matrix.getNcols(); j++) {
                if (j == excluding_col) {
                    continue;
                }
                mat.setValueAt(r, ++c, matrix.getValueAt(i, j));
            }
        }
        return mat;
    }

    public ModMatrix cofactor(ModMatrix matrix) throws Exception {
        ModMatrix mat = new ModMatrix(matrix.getNrows(), matrix.getNcols());
        for (int i = 0; i < matrix.getNrows(); i++) {
            for (int j = 0; j < matrix.getNcols(); j++) {
                mat.setValueAt(i, j, (changeSign(i).multiply(changeSign(j)).multiply(determinant(createSubMatrix(matrix, i, j)))).mod(mod));
            }
        }

        return mat;
    }
    
    public static BigInteger calculateDeterminantUsingLUdecomposition(ModMatrix mymatrix) {
    	
    	BigInteger detResult = BigInteger.ONE;
    	System.out.println("---------------------matrix for LU deccom");
    	printModMatrix(mymatrix);
    	
    	BigInteger[][] zeroMatrix = new BigInteger[mymatrix.nrows][mymatrix.ncols];
    	for (int i = 0; i < zeroMatrix.length; i++) {
    		for (int j = 0; j < zeroMatrix.length; j++) {
    			zeroMatrix[i][j] = BigInteger.ZERO;
    		}
    	}
    	
    	BigInteger[][] zeroMatrix2 = new BigInteger[mymatrix.nrows][mymatrix.ncols];
    	for (int i = 0; i < zeroMatrix.length; i++) {
    		for (int j = 0; j < zeroMatrix.length; j++) {
    			zeroMatrix2[i][j] = BigInteger.ZERO;
    		}
    	}
    	
//    	System.out.println("----------mmm-----------zeroMatrix: ");
//    	for (int i = 0; i < zeroMatrix.length; i++) {
//    		for (int j = 0; j < zeroMatrix.length; j++) {
//    			System.out.print(zeroMatrix[i][j] + " ");
//    		}
//    		System.out.println(" ");
//    	}
    	
    	ModMatrix lower = new ModMatrix(zeroMatrix);
    	ModMatrix upper = new ModMatrix(zeroMatrix2);
//    	for (int i = 0; i < lower.nrows; i++) {
//    		for (int j = 0; j < lower.ncols; j++) {
//    			//System.out.println("--------mmmm-------------i: " + i + "   j: " + j);
//    			lower.setValueAt(i, j, BigInteger.ZERO);
//    		}
//    	}
    	
//    	for (int i = 0; i < upper.nrows; i++) {
//    		for (int j = 0; j < upper.ncols; j++) {
//    			//System.out.println("--------mmmm-------------i: " + i + "   j: " + j);
//    			upper.setValueAt(i, j, BigInteger.ZERO);
//    		}
//    	}
    	
//    	System.out.println("----------mmm-----------upper: ");
//    	printModMatrix(upper);
//    	System.out.println("-----------mmm----------lower: ");
//    	printModMatrix(lower);
    	
    	
        // Decomposing matrix into Upper and Lower
        // triangular matrix
        for (int i = 0; i < mymatrix.nrows; i++) {
     
            // Upper Triangular
            for (int k = i; k < mymatrix.nrows; k++) {
     
                // Summation of L(i, j) * U(j, k)
                BigInteger sum = BigInteger.ZERO;
                for (int j = 0; j < i; j++) {
                    //sum = sum.add(lower.getValueAt(i, j).multiply(upper.getValueAt(j, k)).mod(p)).mod(p);
                	BigInteger lowerIJ =  lower.getValueAt(i, j);
                	BigInteger upperJK =  upper.getValueAt(j, k);
                	BigInteger lowerIJTimesupperJK = lowerIJ.multiply(upperJK).mod(p);
                	sum = sum.add(lowerIJTimesupperJK).mod(p);
                }
                
                //System.out.println("---------------------sum: " + sum);
                // Evaluating U(i, k)
                //upper[i][k] = mat[i][k] - sum;
                //upper.setValueAt(i, k, mymatrix.getValueAt(i, k).subtract(sum).mod(p));
                BigInteger temp2 = mymatrix.getValueAt(i, k);
                temp2 = temp2.subtract(sum).mod(p);
                System.out.println("---------------------temp2: " + temp2);
                upper.setValueAt(i, k, temp2);
            }
            
//        	System.out.println("---------------------upper: ");
//        	printModMatrix(upper);
//        	System.out.println("---------------------lower: ");
//        	printModMatrix(lower);
            
     
            // Lower Triangular
            for (int k = i; k < mymatrix.nrows; k++) {
                if (i == k)
                    lower.setValueAt(i, i, BigInteger.ONE); // Diagonal as 1
                else {
     
                    // Summation of L(k, j) * U(j, i)
                    //int sum = 0;
                	BigInteger sum = BigInteger.ZERO;
                    for (int j = 0; j < i; j++) {
                        //sum += (lower[k][j] * upper[j][i]);
                    	//sum = sum.add(lower.getValueAt(k, j).multiply(upper.getValueAt(j, i)).mod(p)).mod(p);
                    	BigInteger lowerKJ = lower.getValueAt(k, j);
                    	BigInteger upperJI = upper.getValueAt(j, i);
                    	BigInteger lowerKJTimesupperJI = lowerKJ.multiply(upperJI).mod(p);
                    	sum = sum.add(lowerKJTimesupperJI).mod(p);
                    	
                    	System.out.println("lower[" + k + "][" + j + "]: " + lower.getValueAt(k, j));
                    }
                    
                    System.out.println("");
                    // Evaluating L(k, i)
                    //lower[k][i] = (mat[k][i] - sum) / upper[i][i];
                    BigInteger upperii = upper.getValueAt(i, i);
                    System.out.println("---------------------upperii: " + upperii);
                    //BigInteger upperiiInverse = upperii.modInverse(p);

                    BigInteger temp = mymatrix.getValueAt(k, i);
                    temp = temp.subtract(sum);
                    //temp = temp.multiply(upperiiInverse).mod(p);
                    temp = temp.divide(upperii).mod(p);
                    System.out.println("---------------------temp: " + temp);
                    lower.setValueAt(k, i, temp);
                    //lower.setValueAt(k, i, mymatrix.getValueAt(k, i).subtract(sum).mod(p).multiply(upper.getValueAt(i, i).modInverse(p)));
                }
            }
            
        	System.out.println("----------nnn-----------upper: ");
        	printModMatrix(upper);
        	System.out.println("-----------nnn----------lower: ");
        	printModMatrix(lower);
        }
        
//    	System.out.println("---------------------upper: ");
//    	printModMatrix(upper);
//    	System.out.println("---------------------lower: ");
//    	printModMatrix(lower);
        
        BigInteger lowerDet = BigInteger.ONE;
        BigInteger upperDet = BigInteger.ONE;
    	for (int i = 0; i < lower.getNrows(); i++) {
    		lowerDet = detResult.multiply(lower.getValueAt(i, i)).mod(p);
    		upperDet = detResult.multiply(upper.getValueAt(i, i)).mod(p);
    	}
        detResult = lowerDet.multiply(upperDet).mod(p);
        
    	return detResult;
    }
    
    
    public static BigInteger calculateDeterminantUsingLUdecomposition00(ModMatrix mymatrix) {
    	
//    	System.out.println("-------------------mymatrix.nrows: " + mymatrix.nrows);
//    	System.out.println("---------------------matrix.nrows: " + mymatrix.nrows + "   matrix.ncols: " + mymatrix.ncols);
    	
    	BigInteger detResult = BigInteger.ONE;
    	System.out.println("---------------------matrix for LU deccom");
    	printModMatrix(mymatrix);
    	
    	//BigInteger[][] zeroMatrix = new BigInteger[mymatrix.nrows][mymatrix.ncols];
    	BigInteger[][] zeroMatrix = new BigInteger[mymatrix.nrows][mymatrix.ncols];
    	for (int i = 0; i < zeroMatrix.length; i++) {
    		for (int j = 0; j < zeroMatrix.length; j++) {
    			//System.out.println("--------mmmm-------------i: " + i + "   j: " + j);
    			zeroMatrix[i][j] = BigInteger.ZERO;
    		}
    	}
    	
    	BigInteger[][] zeroMatrix2 = new BigInteger[mymatrix.nrows][mymatrix.ncols];
    	for (int i = 0; i < zeroMatrix.length; i++) {
    		for (int j = 0; j < zeroMatrix.length; j++) {
    			//System.out.println("--------mmmm-------------i: " + i + "   j: " + j);
    			zeroMatrix2[i][j] = BigInteger.ZERO;
    		}
    	}
    	
    	System.out.println("----------mmm-----------zeroMatrix: ");
    	for (int i = 0; i < zeroMatrix.length; i++) {
    		for (int j = 0; j < zeroMatrix.length; j++) {
    			System.out.print(zeroMatrix[i][j] + " ");
    		}
    		System.out.println(" ");
    	}
    	
    	ModMatrix lower = new ModMatrix(zeroMatrix);
    	ModMatrix upper = new ModMatrix(zeroMatrix2);
//    	for (int i = 0; i < lower.nrows; i++) {
//    		for (int j = 0; j < lower.ncols; j++) {
//    			//System.out.println("--------mmmm-------------i: " + i + "   j: " + j);
//    			lower.setValueAt(i, j, BigInteger.ZERO);
//    		}
//    	}
    	
    	for (int i = 0; i < upper.nrows; i++) {
    		for (int j = 0; j < upper.ncols; j++) {
    			//System.out.println("--------mmmm-------------i: " + i + "   j: " + j);
    			upper.setValueAt(i, j, BigInteger.ZERO);
    		}
    	}
    	
    	System.out.println("----------mmm-----------upper: ");
    	printModMatrix(upper);
    	System.out.println("-----------mmm----------lower: ");
    	printModMatrix(lower);
    	
//    	System.out.println("---------------------lower: ");
//    	printModMatrix(lower);
    	
    	BigInteger temp;
        // Decomposing matrix into Upper and Lower
        // triangular matrix
        for (int i = 0; i < mymatrix.nrows; i++) {
     
            // Upper Triangular
            for (int k = i; k < mymatrix.nrows; k++) {
     
                // Summation of L(i, j) * U(j, k)
                BigInteger sum = BigInteger.ZERO;
                for (int j = 0; j < i; j++) {
                    sum = sum.add(lower.getValueAt(i, j).multiply(upper.getValueAt(j, k)).mod(p)).mod(p);
                }
                
                //System.out.println("---------------------sum: " + sum);
                // Evaluating U(i, k)
                //upper[i][k] = mat[i][k] - sum;
                //upper.setValueAt(i, k, mymatrix.getValueAt(i, k).subtract(sum).mod(p));
                BigInteger temp2 = mymatrix.getValueAt(i, k);
                temp2 = temp2.subtract(sum).mod(p);
                upper.setValueAt(i, k, temp2);
            }
            
//        	System.out.println("---------------------upper: ");
//        	printModMatrix(upper);
//        	System.out.println("---------------------lower: ");
//        	printModMatrix(lower);
            
     
            // Lower Triangular
            for (int k = i; k < mymatrix.nrows; k++) {
                if (i == k)
                    lower.setValueAt(i, i, BigInteger.ONE); // Diagonal as 1
                else {
     
                    // Summation of L(k, j) * U(j, i)
                    //int sum = 0;
                	BigInteger sum = BigInteger.ZERO;
                    for (int j = 0; j < i; j++) {
                        //sum += (lower[k][j] * upper[j][i]);
                    	sum = sum.add(lower.getValueAt(k, j).multiply(upper.getValueAt(j, i)).mod(p)).mod(p);
                    }
     
                    // Evaluating L(k, i)
                    //lower[k][i] = (mat[k][i] - sum) / upper[i][i];
                    BigInteger upperii = upper.getValueAt(i, i).modInverse(p);
                    System.out.println("**************************upperii: " + upperii);
                    BigInteger upperiiInverse = upperii.modInverse(p);
                    temp = mymatrix.getValueAt(k, i).subtract(sum).mod(p).multiply(upperiiInverse);
                    lower.setValueAt(k, i, temp);
                    //lower.setValueAt(k, i, mymatrix.getValueAt(k, i).subtract(sum).mod(p).multiply(upper.getValueAt(i, i).modInverse(p)));
                }
            }
            
//        	System.out.println("----------nnn-----------upper: ");
//        	printModMatrix(upper);
//        	System.out.println("-----------nnn----------lower: ");
//        	printModMatrix(lower);
        }
        
    	System.out.println("---------------------upper: ");
    	printModMatrix(upper);
    	System.out.println("---------------------lower: ");
    	printModMatrix(lower);        
    	
        BigInteger lowerDet = BigInteger.ONE;
        BigInteger upperDet = BigInteger.ONE;
    	for (int i = 0; i < lower.getNrows(); i++) {
    		lowerDet = detResult.multiply(lower.getValueAt(i, i)).mod(p);
    		upperDet = detResult.multiply(upper.getValueAt(i, i)).mod(p);
    	}
        detResult = lowerDet.multiply(upperDet).mod(p);
        
    	return detResult;
    }
    
    public static BigInteger calculateDetFastUsingReducedRowEchelonForm(ModMatrix matrix) throws Exception {
    	BigInteger detResult = BigInteger.ONE;
    	
    	ModMatrix matrixReducedRowForm = gaussianSolve(matrix);
    	System.out.println("matrixReducedRowForm: ");
    	printModMatrix(matrixReducedRowForm);
    	
    	for (int i = 0; i < matrixReducedRowForm.getNrows(); i++) {
    		detResult = detResult.multiply(matrixReducedRowForm.getValueAt(i, i)).mod(p);
    		System.out.println("detResult: " + detResult + "      matrixReducedRowForm.getValueAt(i, i): " + matrixReducedRowForm.getValueAt(i, i));
    	}
    	
    	return detResult;
    }
    
    // this method produces the row reduced echelon form of a matrix, which get as input:
    public static ModMatrix gaussianSolve(ModMatrix matrix) throws Exception{   
        //This method only works when the modulus is prime   
        if (!p.isProbablePrime(16)) throw new IllegalArgumentException("Gaussian elimination method currently requires modulus to be prime!");
        
        //Work the rows, starting with the first row   
        int currentRow=1;   
        while (currentRow < matrix.nrows) {   
           int i=currentRow;   
           //Make sure diagonal element is nonzero, if possible, by swapping   
           //while (i<=matrix.ncols && matrix.array[i][currentRow].equals(BigInteger.ZERO)) i++;
           while (i<matrix.ncols && matrix.getValueAt(i, currentRow).equals(BigInteger.ZERO)) i++;
           
           //System.out.println("*******************i ? mat.numRows:" + i + "  and  " + mat.numRows);
           if (i > matrix.nrows) throw new Exception("Linear dependence exists here, i.e. determinant is zero!");
           
           if (i >= matrix.nrows) {
               System.out.println("*******************Matrix mat (GaussianSolve Result):");
               printModMatrix(matrix);
               System.out.println("****************************************************");
           }
           //Swap with a row not having zero in diagonal position   
           if (currentRow != i) swapRows(matrix, currentRow, i);   
           //Now, you must produce all zeros below and above the diagonal element   
           i = 0;   
           //Multiply each row by the proper scalar
           System.out.println("*******************i and matrix.nrows:" + i + "   " + matrix.nrows);
           while (i < matrix.nrows) {
        	   System.out.println("*******************matrix.nrows, i and currentRow:" + matrix.nrows + "    " + i + "   " + currentRow);
        	   printModMatrix(matrix);
              if (i!=currentRow) {
                 BigInteger scalar=matrix.getValueAt(i, currentRow);  
                 if (!scalar.equals(BigInteger.ZERO)) {   
                    multiplyRow(matrix,i,matrix.getValueAt(currentRow, currentRow));   
                    multiplyRow(matrix,currentRow,scalar);   
                    //Replace row i with row i minus diagonal row   
                    subtractRow(matrix,i,currentRow);   
                 }   
              }   
              i++;   
           }   
           currentRow++;   
        }   
//        //Now, produce 1's along main diagonal by multiplying by an inverse   
//        for (int index=0;index< matrix.nrows;index++) {   
//           multiplyRow(matrix,index,matrix.getValueAt(index, index).modInverse(p));   
//        }   
        //Remember, b may be a square matrix-polymorphism takes care of this here
        return matrix;   
     }
   
    //Used by gaussianSolve to swap two rows   
    private static void swapRows(ModMatrix mat, int r1,int r2) {   
       BigInteger temp;
       System.out.println("*******************r1:" + r1 + "        " + "r2: " + r2);
       for (int j=0 ;j< mat.ncols;j++) {   
          temp=mat.getValueAt(r1, j);
          //mat.array[r1][j] = mat.array[r2][j];
          mat.setValueAt(r1, j, mat.getValueAt(r2, j));
          //mat.array[r2][j]=temp;
          mat.setValueAt(r2, j, temp);
       }  
    } 
    
    //Used by gaussianSolve to multiply a row by some scalar   
    private static void multiplyRow(ModMatrix mat,int i,BigInteger scalar) {   
       //Multiplies row i by scalar-answer replaces i-th row
    	//for (int k=1;k<=mat.numCols;k++)
    		//mat.array[i][k]=BigIntegerMathMatrixOperations.lnr(mat.array[i][k].multiply(scalar),mat.modulus);
       for (int k=0; k<mat.ncols; k++)
    	   mat.setValueAt(i, k, mat.getValueAt(i, k).multiply(scalar).mod(p));      
    }
    
    //Used by gaussianSolve to subtract one row from another   
    private static void subtractRow(ModMatrix mat,int i,int j) {   
       //Subtracts row j from row i; answer replaces row i
    	//mat.array[i][k]=BigIntegerMathMatrixOperations.lnr(mat.array[i][k].subtract(mat.array[j][k]),mat.modulus);
       for (int k=0; k<mat.ncols; k++)
    	   mat.setValueAt(i, k, mat.getValueAt(i, k).subtract(mat.getValueAt(j, k)).mod(p));
    }
    
    public ModMatrix inverse(ModMatrix matrix) throws Exception {
//    	System.out.println(">>>>>>>>>>>>>>>>>>>>>>matrix in inverse func: ");
//    	printModMatrix(matrix);
        //return (transpose(cofactor(matrix)).dc(determinant(matrix)));
    	
    	//saveToFile(">>>>>>>>>>>>>>>>>>>>>>matrix in inverse func: ");
    	saveToFileModMatrix(matrix);
    	
    	BigInteger det = determinant(matrix);
    	//System.out.println(">>>>>>>>>>>>>>>>>>>>>>matrix det in inverse func: " + det);
    	//saveToFile(">>>>>>>>>>>>>>>>>>>>>>matrix det in inverse func: " + det);
    	return (transpose(cofactor(matrix)).dc(det));
    }

    private ModMatrix dc(BigInteger d) {
        BigInteger inv = d.modInverse(mod);
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                data[i][j] = (data[i][j].multiply(inv)).mod(mod); 
            }
        }
        return this;
    }
    
    ////////////////////////////////////////////////////////////////////////my functions
    public BigInteger[] multiplyMatrixByVectorBigInt(BigInteger[][] matrix, BigInteger[] vector){
    	
    	int rows = matrix.length;
        int columns = matrix[0].length;
        
        BigInteger[] result = new BigInteger[rows];
        
        for (int row = 0; row < rows; row++) {
            BigInteger sum = BigInteger.valueOf(0);
            for (int column = 0; column < columns; column++) {
//            	System.out.println("matrix[row][column]: " + matrix[row][column] +
//            					   "    vector[column]: " + vector[column] + 
//            					   "    mult result: " + matrix[row][column].multiply(vector[column]));
                sum = sum.add(matrix[row][column].multiply(vector[column]).mod(p)).mod(p);
            }
            //System.out.println("sum: " + sum);
            result[row] = sum;
        }
    	
    	return result;
    }
    
    
    // print ModMatrix
    public static void printModMatrix(ModMatrix matrix) {
	    for (int i = 0; i < matrix.getNrows(); i++) {
	    	System.out.print("{");
	        for (int j = 0; j < matrix.getNcols() - 1; j++) {
	        	System.out.print(matrix.getData()[i][j]+", ");
	        }
	        System.out.print(matrix.getData()[i][matrix.getNcols() - 1] + "},");
	        System.out.println("");
	    }
    }
    
    // print ModMatrix
    public static void printModMatrix00(ModMatrix matrix) {
	    for (int i = 0; i < matrix.getNrows(); i++) {
	        for (int j = 0; j < matrix.getNcols(); j++) {
	        	System.out.print(matrix.getData()[i][j]+" ");
	        }
	        System.out.println("");
	    }
    }
    
    // save to file ModMatrix
    public static void saveToFileModMatrix(ModMatrix matrix) {
    	
    	try {
            PrintStream out = new PrintStream(new FileOutputStream("OutFile.txt", true));
            
    	    for (int i = 0; i < matrix.getNrows(); i++) {
    	        for (int j = 0; j < matrix.getNcols(); j++) {
    	        	//out.print(matrix.getData()[i][j]+" ");
    	        	out.printf("%4d", matrix.getData()[i][j]);
    	        }
    	        out.println("");
    	    }

            out.close();
          } catch (FileNotFoundException e) {
            e.printStackTrace();
          }
    }
    
    // save to file ModMatrix
    public static void saveToFileBigIntMatrix(BigInteger[][] matrix) {
    	
    	try {
            PrintStream out = new PrintStream(new FileOutputStream("OutFile.txt", true));
            
    	    for (int i = 0; i < matrix.length; i++) {
    	        for (int j = 0; j < matrix.length; j++) {
    	        	//out.print(matrix.getData()[i][j]+" ");
    	        	out.printf("%4d", matrix[i][j]);
    	        }
    	        out.println("");
    	    }

            out.close();
          } catch (FileNotFoundException e) {
            e.printStackTrace();
          }
    }
    
    // print biginteger matrix
    public void printMatrixBigInt(BigInteger[][] matrix) {
	    for (int i = 0; i < matrix.length; i++) {
	        for (int j = 0; j < matrix.length; j++) {
	                System.out.print(matrix[i][j]+" ");
	        }
	        System.out.println("");
	    }
    }
    
    public static void saveToFile(String str) {
        try {
          PrintStream out = new PrintStream(new FileOutputStream("OutFile.txt", true));
          out.println(str);

          out.close();
        } catch (FileNotFoundException e) {
          e.printStackTrace();
        }
      }

    static BigInteger dotProduct(BigInteger[] a, BigInteger[] b) {
        //return range(0, a.length).mapToDouble(i -> a[i] * b[i]).sum();
    	BigInteger sum = BigInteger.ZERO;
    	for (int i = 0; i < a.length; i++)
    		sum = sum.add(a[i].multiply(b[i]));
    	
    	return sum;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    static BigInteger[][] matrixMul(BigInteger[][] A, BigInteger[][] B) {
    	BigInteger[][] result = new BigInteger[A.length][B[0].length];
    	BigInteger[] aux = new BigInteger[B.length];
 
        for (int j = 0; j < B[0].length; j++) {
 
            for (int k = 0; k < B.length; k++)
                aux[k] = B[k][j];
 
            for (int i = 0; i < A.length; i++)
                result[i][j] = dotProduct(A[i], aux);
        }
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    static BigInteger[][] pivotize(BigInteger[][] m) {
        int n = m.length;
        PmatrixDet = 1;
        double[][] id = range(0, n).mapToObj(j -> range(0, n)
                .mapToDouble(i -> i == j ? 1 : 0).toArray())
                .toArray(double[][]::new);
        
        BigInteger[][] id2 = new BigInteger[n][n];;
        
        for (int i = 0; i < n; i++) {
        	BigInteger maxm = m[i][i];
            int row = i;
            for (int j = i; j < n; j++)
                if(m[j][i].compareTo(maxm) > 0) {//if (m[j][i] > maxm) {
                    maxm = m[j][i];
                    row = j;
                }
 
            if (i != row) {
            	double[] tmp = id[i];
                id[i] = id[row];
                id[row] = tmp;
                PmatrixDet *= -1; // each row (or column) permutation changes the sign of determinant 
            }
        }
        for (int i = 0; i < n; i++) {
        	for (int j = 0; j < n; j++) {
        		if(id[i][j] == 1.0)
        			id2[i][j] = BigInteger.ONE;
        		else if(id[i][j] == 0.0)
        			id2[i][j] = BigInteger.ZERO;
        		else
        			System.out.println("element not 1.0 nor 0.0");
        	}
        }
        return id2;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    static BigInteger[][][] lu(BigInteger[][] A) {
        int n = A.length;
        BigInteger[][] L = new BigInteger[n][n];
        BigInteger[][] U = new BigInteger[n][n];
        BigInteger[][] P = pivotize(A);
        BigInteger[][] A2 = matrixMul(P, A);
        
        for (int i = 0; i < n; i++) {
        	for (int j = 0; j < n; j++) {
        		L[i][j] = BigInteger.ZERO;
        		U[i][j] = BigInteger.ZERO;
        	}
        }
        
//        System.out.println("L: ");
//        print(L);
        for (int j = 0; j < n; j++) {
            L[j][j] = BigInteger.ONE;
            for (int i = 0; i < j + 1; i++) {
            	BigInteger s1 = BigInteger.ZERO;
                for (int k = 0; k < i; k++)
                    s1 = s1.add(U[k][j].multiply(L[i][k])); //s1 += U[k][j] * L[i][k];
                U[i][j] = A2[i][j].subtract(s1).mod(p);//U[i][j] = A2[i][j] - s1;
                //System.out.println("U[" + i + "][" + j + "]: " + U[i][j]);
            }
            for (int i = j; i < n; i++) {
            	BigInteger s2 = BigInteger.ZERO;
                for (int k = 0; k < j; k++)
                    s2 = s2.add(U[k][j].multiply(L[i][k]));//s2 += U[k][j] * L[i][k];
                //L[i][j] = (A2[i][j].subtract(s2)).divide(U[j][j]);//L[i][j] = (A2[i][j] - s2) / U[j][j];
                
                // when U[j][j] = 0, this does not work, probably det(A) zero ...
                //L[i][j] = (A2[i][j].subtract(s2).mod(primeP)).multiply(U[j][j].modInverse(primeP)).mod(primeP);//L[i][j] = (A2[i][j] - s2) / U[j][j];
                
                if (U[j][j].equals(BigInteger.ZERO) == true)
                	L[i][j] = BigInteger.ZERO;
                else
                	L[i][j] = (A2[i][j].subtract(s2).mod(p)).multiply(U[j][j].modInverse(p)).mod(p);
                
                
                //System.out.println("L[" + i + "][" + j + "]: " + L[i][j]);
            }
        }
        return new BigInteger[][][]{L, U, P};
    }
 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    static void print(BigInteger[][] A) {
    	int n = A.length;
//        stream(m).forEach(a -> {
//            stream(a).forEach(n -> System.out.printf(Locale.US, "%5.1f ", n));
//            System.out.println();
//        });
        for (int i = 0; i < n; i++) {
        	for (int j = 0; j < n; j++) {
        			System.out.print(A[i][j] + " ");
        	}
        	System.out.println();
        }
        System.out.println();
    }
 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    public BigInteger calculateDetUsingLUdecom(BigInteger[][] mat) {
    	BigInteger detU = BigInteger.ONE;
    	//saveToFile("----------detU --- in calculateDetUsingLUdecom func at beiginning: " + detU);
		// getting L, U and P
		BigInteger[][][] LUP = null;
		LUP = lu(mat);
		//BigInteger[][] L = LUP[0];
		BigInteger[][] U = LUP[1];
		//BigInteger[][] P = LUP[2];
		
//		saveToFile("L, U, P:");
//		saveToFile("L:");
//		saveToFileBigIntMatrix(L);
//		saveToFile("U:");
//		saveToFileBigIntMatrix(U);
//		saveToFile("P:");
//		saveToFileBigIntMatrix(P);
		
		
		// calculating the determinant of A (original matrix) using LU and det of P:
		//BigInteger detU = BigInteger.ONE;
		for (int i = 0; i < U.length; i++) {
			//saveToFile("U[" + i + "][" + i + "]:" + U[i][i] + "      detU: " + detU);
			detU = detU.multiply(U[i][i]);
			detU = detU.mod(p);
		}
		//saveToFile("----------detU --- in calculateDetUsingLUdecom: " + detU);
		if (PmatrixDet == -1) {
			//saveToFile("----------PmatrixDet --- in calculateDetUsingLUdecom: " + PmatrixDet);
			detU = p.subtract(detU).mod(p);
		}
		
		//saveToFile("----------detU --- in calculateDetUsingLUdecom func at end: " + detU);
		//saveToFile("-------det with LU: " + detU);
		
		return detU;

    }

}
