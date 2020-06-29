// this seems works

import static java.util.Arrays.stream;
import java.util.Locale;
import static java.util.stream.IntStream.range;

import java.math.BigInteger;
 
public class luDecompositionBigIntModP {
 
	static BigInteger primeP = new BigInteger("13");
	static int PmatrixDet = 1;
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
	public luDecompositionBigIntModP() {
		
	}
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
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
                U[i][j] = A2[i][j].subtract(s1).mod(primeP);//U[i][j] = A2[i][j] - s1;
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
                	L[i][j] = (A2[i][j].subtract(s2).mod(primeP)).multiply(U[j][j].modInverse(primeP)).mod(primeP);
                
                
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
    public static BigInteger calculateDetUsingLUdecom(BigInteger[][] mat) {
    	BigInteger detU = BigInteger.ONE;
		// getting L, U and P
		BigInteger[][][] LUP = null;
		LUP = lu(mat);
		BigInteger[][] L = LUP[0];
		BigInteger[][] U = LUP[1];
		BigInteger[][] P = LUP[2];
		
//		System.out.println("L, U, P: ");
//		System.out.println("L:");
//		print(L);
//		System.out.println("U:");
//		print(U);
//		System.out.println("P:");
//		print(P);
		
		// calculating the determinant of A (original matrix) using LU and det of P:
		//BigInteger detU = BigInteger.ONE;
		for (int i = 0; i < U.length; i++) {
			detU = detU.multiply(U[i][i]).mod(primeP);
		}
		if (PmatrixDet == -1)
			detU = primeP.subtract(detU);
		
//		System.out.println("detU:" + detU);
		
		return detU;

    }
}