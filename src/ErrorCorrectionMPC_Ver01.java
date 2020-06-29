import java.lang.*;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.jlinalg.Matrix;
import org.jlinalg.Vector;
import org.jlinalg.field_p.*;
import org.jlinalg.field_p.FieldP;
import org.jlinalg.field_p.FieldPAbstractFactory;
import org.jlinalg.field_p.FieldPFactoryMap;
//import org.jlinalg.polynomial.Polynomial;
//import org.jlinalg.polynomial.PolynomialFactory;
import org.jlinalg.rational.Rational;
import org.jlinalg.rational.Rational.RationalFactory;
import org.jscience.mathematics.function.Polynomial;
import org.jscience.mathematics.function.Variable;
import org.jscience.mathematics.number.Complex;
//import org.jscience.mathematics.function.Polynomial;

public class ErrorCorrectionMPC_Ver01 {
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	static int threshold = 3;
	static int numOfErrors = 1;
	static int numOfPlayers = 2 * numOfErrors + threshold;
	
	public static void main(String[] args) {
		findErrorLocations();
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	@SuppressWarnings("unchecked")
	public static void findErrorLocations() {
		
		// The PField require the creation of a factory. Here, a factory for G(7) is used.
		FieldPAbstractFactory<?> factory7 = FieldPFactoryMap.getFactory(new Long(7));

		// By the means of the factory, instances are created.
		FieldP<?> p1, p2, p3, p4, p5, p6;
		p1 = factory7.get(1);
		p2 = factory7.get(2);
		p3 = factory7.get(3);
		p4 = factory7.get(4);
		p5 = factory7.get(5);
		p6 = factory7.get(6);
		
		int alpha1 = 3;
		int alpha2 = 0 + 1;
		int alpha3 = 6;
		int alpha4 = 0;
		int alpha5 = 3;
		
		Matrix<?> matrixA = new Matrix(new FieldP[][]
				{ {factory7.get(1), factory7.get(1),    factory7.get(1),    factory7.get(1),    factory7.get(-1*alpha1)},
				  {factory7.get(1), factory7.get(Math.pow(2,1)),    factory7.get(Math.pow(2,2)),    factory7.get(Math.pow(2,3)),    factory7.get(-1*alpha2)},
				  {factory7.get(1), factory7.get(Math.pow(3,1)),    factory7.get(Math.pow(3,2)),    factory7.get(Math.pow(3,3)),    factory7.get(-1*alpha3)},
				  {factory7.get(1), factory7.get(Math.pow(4,1)),    factory7.get(Math.pow(4,2)),    factory7.get(Math.pow(4,3)),    factory7.get(-1*alpha4)},
				  {factory7.get(1), factory7.get(Math.pow(5,1)),    factory7.get(Math.pow(5,2)),    factory7.get(Math.pow(5,3)),    factory7.get(-1*alpha5)},
				});
		
		System.out.println("matrixA:\n" + matrixA);
		
		Vector<?> c = new Vector(new FieldP[] {
				factory7.get(alpha1), factory7.get(2*alpha2), factory7.get(3*alpha3), factory7.get(4*alpha4), factory7.get(5*alpha5)});
		
		System.out.println("\nc:\n" + c);
		
		System.out.println("\nmatrixA_inv:\n" + matrixA.inverse());
		
		Vector<?> sysOfEqSolutions = (Vector) ((Matrix) matrixA.inverse()).multiply(c);
		
		System.out.println("\nmatrixA_inv . c = " + sysOfEqSolutions);
		
		// constructing the error locator polynomial, E(x), we need to extract b0 from sysOfEqSolutions:
		FieldP<?> b0 = factory7.get(0);
		b0 = (FieldP<?>) sysOfEqSolutions.getEntry(numOfPlayers);
		
		System.out.println("\nb0 = " + b0);
		
//		//Then, the polynomial "1 + 2*x^2 + 2*x^3" can be created as follows:
//		PolynomialFactory<Rational> rationalPolyFactory = PolynomialFactory.getFactory(Rational.FACTORY);
//		Map<Integer, Rational> coeff1 = new HashMap<Integer, Rational>();
//		coeff1.put(0, Rational.FACTORY.get(1));
//		coeff1.put(1, Rational.FACTORY.get(2));
//		coeff1.put(2, Rational.FACTORY.get(3));
//		Polynomial<Rational> polynomial = rationalPolyFactory.get(coeff1);
//		
//		System.out.println("\npoly: " + polynomial);
		
		
//	    // Defines two local variables (x, y).
//	    Variable<Complex> varX = new Variable.Local<Complex>("x");
//	    Variable<Complex> varY = new Variable.Local<Complex>("y");
//		Polynomial<Complex> x = Polynomial.valueOf(Complex.ONE, varX);
//	    Polynomial<Complex> fx = x.pow(2).times(Complex.I).plus(x.times(Complex.valueOf(2, 0)).plus(Complex.ONE));
		
	    // Defines two local variables (x, y).
	    Variable<Complex> varX = new Variable.Local<Complex>("x");
	    Variable<Complex> varY = new Variable.Local<Complex>("y");
		Polynomial<Complex> x = Polynomial.valueOf(Complex.ONE, varX);
	    Polynomial<Complex> fx = x.pow(2).times(Complex.I).plus(x.times(Complex.valueOf(2, 0)).plus(Complex.ONE));
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// finite field example using JLinAlg: http://jlinalg.sourceforge.net/
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// solving linear system of equation using Jama: https://math.nist.gov/javanumerics/jama/
	/**
	public static void JamaTest(){
		Matrix A = new Matrix(new double[][] { { 0,  1,  1 },
											   { 2,  4, -2 },
											   { 0,  3, 15 } }); 

		Matrix b = new Matrix(new double[][] { {  4 },
		              {  2 },
		              { 36 } });
		Matrix x = A.solve(b);
		Matrix residual = A.times(x).minus(b);
		double rnorm = residual.normInf();
		
		StdOut.println("A");
		A.print(9, 6);                // printf("%9.6f");
		
		StdOut.println("b");
		b.print(9, 6);
		
		StdOut.println("x");
		x.print(9, 6);
		StdOut.println(x);
		
		StdOut.println("residual");
		residual.print(9, 6);
	}*/
	
}
