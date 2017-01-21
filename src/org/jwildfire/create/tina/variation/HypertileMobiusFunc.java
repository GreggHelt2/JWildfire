package org.jwildfire.create.tina.variation;

import static java.lang.Math.random;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.base.Layer;

// import org.jwildfire.base.mathlib.Complex;
import org.nfunk.jep.type.Complex;

/* using org.nfunk.jep.type.Complex for complex numbers 
//     instead of org.jwildfire.base.mathlib.Complex because
//     its methods return new complex numbers (instead of modifying in place like ...mathlib.Complex)
//     which makes chaining methods together for long calcs possible 
// there is a performance hit relative to ...mathlib.Complex, since creating lots of 
//     intermediary complex objects. But compared performance within transform() loop to 
//     using ...mathlib.Complex instead, and only get ~20% slower. If really becomes a problem 
//     could later convert to using ...mathlib.Complex
*/    

/**
 * Given positive integers K, L, M such that (1/K + 1/L + 1/M) < 1
 * creates a hyperbolic tiling where tiled triangles have angles of 
 * k = pi/k, l = pi/l, m = pi/m and 
 * triangle veritices have point symmetry of K, L, M
 * math based on Norberto Gavioli's triangular_generators.py script
 * internally uses two Mobius transformations 
 *   (and calculates a third that could be used for the "untiled" part)
 * Author: Gregg Helt
 */
public class HypertileMobiusFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_K = "K";
  private static final String PARAM_L = "L";
  private static final String PARAM_M = "M";

  private static final String[] paramNames = { PARAM_K, PARAM_L, PARAM_M };


  // flattened matrices for mapping
  //    [0, 1, 2, 3] = [a, b, c, d]  ==> f(z)= (az+b)/(cz+d)
  protected Complex[] mat_a = new Complex[4];
  protected Complex[] mat_b = new Complex[4];
  protected Complex[] mat_c = new Complex[4];
  // matrix inversion (not yet needed)
  // protected Complex[] mat_a_inv = new Complex[4];
  // protected Complex[] mat_b_inv = new Complex[4];
  // protected Complex[] mat_c_inv = new Complex[4];

  protected Object[] mtransforms;
  
  // K, L, M should be integers, and 1/K + 1/L + 1/M should be < 1
  //    but for now leaving as double to see what happens if set them non-integer
  protected double K = 3;
  protected double L = 4;
  protected double M = 5;
  protected Complex kp, lm, mu;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    // randomly pick one of the two calculated Mobius transform matrices, a, b
    int mindex = pContext.random(2);
    Complex[] mat = (Complex[]) mtransforms[mindex];
    // then use selected matrix for Mobius transformation:
    //  f(z) = (az + b) / (cz + d);

    Complex a = mat[0];
    Complex b = mat[1];
    Complex c = mat[2];
    Complex d = mat[3];
    
    Complex zin = new Complex(pAffineTP.x, pAffineTP.y);
    Complex zout = zin.mul(a).add(b).div(zin.mul(c).add(d));
    
    pVarTP.x += pAmount * zout.re();
    pVarTP.y += pAmount * zout.im();
    
    if (pContext.isPreserveZCoordinate()) {
      pVarTP.z += pAmount * pAffineTP.z;
    }
  }
  
  static Complex re1 = new Complex(1, 0);
  static Complex im1 = new Complex(0, 1);
  static Complex zero = new Complex(0, 0);
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    /* From Norberto Gavioli's triangular_generators.py script
    kp=exp(-1j*pi/k)
    lm=exp(1j*pi/l)
    mu=exp(1j*pi/m)

    rho=(( (mu+1/mu+lm*kp+1/(lm*kp))/(mu+1/mu+lm/kp+kp/lm).real   )**.5).real
    r=1/rho-rho
    
    a=[[kp,0],[0,kp.conjugate()]]
    b11=rho*lm.conjugate()/r-lm/(r*rho)
    b12=(lm-lm.conjugate())/r
    b=[[b11,b12],[b12.conjugate(), b11.conjugate()]]
    c=matprod2(a,b)
    */
    // trying to set kp = e^(-i*pi/k)    
    // kp = new Complex(Math.E, 0);
    // kp.power(new Complex(0, -1 * M_PI / K));
    // not sure why above doesn't work, but testing it get kp = e instead
    // so trying different way, based on Euler's formulas
    //      e^( i*t) = cos(t) + i*sin(t)
    //      e^(-i*t) = cos(t) - i*sin(t)
    // which works...
    kp = new Complex(Math.cos(M_PI/K), -1 * Math.sin(M_PI/K));
    
    // lm = new Complex(Math.E, 0);
    // lm.power(new Complex(0, M_PI / L));
    lm = new Complex(Math.cos(M_PI/L), Math.sin(M_PI/L));
    
    // mu = new Complex(Math.E, 0);
    // mu.power(new Complex(0, M_PI / M));
    mu = new Complex(Math.cos(M_PI/M), Math.sin(M_PI/M));
    
    // rho=(( (mu + 1/mu + lm*kp + 1/(lm*kp)) / (mu+ 1/mu + lm/kp + kp/lm).real   )**.5).real
    // r=1/rho-rho
    Complex rhonum =   mu.add(re1.div(mu)).add(lm.mul(kp)).add(re1.div(lm.mul(kp)));
    Complex rhodenom = (mu.add(re1.div(mu)).add(lm.div(kp)).add(kp.div(lm)));
    double realdenom = rhodenom.re();
    double rho = rhonum.mul(1/realdenom).power(0.5).re();
    double r = 1/rho - rho;
    
    // a=[[kp,0],[0,kp.conjugate()]]
    mat_a[0] = kp;
    mat_a[1] = zero;
    mat_a[2] = zero;
    mat_a[3] = kp.conj();
    
    // b11 = rho*lm.conjugate()/r - lm/(r*rho)
    //  b12=(lm-lm.conjugate())/r 
    // b=[[b11,b12],[b12.conjugate(), b11.conjugate()]]
    mat_b[0] = lm.conj().mul(rho).mul(1/r).sub(lm.mul(1/(r*rho)));            
    mat_b[1] = lm.sub(lm.conj()).mul(1/r);
    mat_b[2] = mat_b[1].conj();
    mat_b[3] = mat_b[0].conj();

    // matrix C is composition (matrix product) of A and B
    mat_c[0] = mat_a[0].mul(mat_b[0]).add(mat_a[1].mul(mat_b[2]));
    mat_c[1] = mat_a[0].mul(mat_b[1]).add(mat_a[1].mul(mat_b[3]));
    mat_c[2] = mat_a[2].mul(mat_b[0]).add(mat_a[3].mul(mat_b[2]));
    mat_c[3] = mat_a[2].mul(mat_b[1]).add(mat_a[3].mul(mat_b[3]));
    
    // since c is 0 in mat_a, get Infinity (or NaN) returned form getMobiusFixedPoints()
    Complex[] fixed_pointsA = getMobiusFixedPoints(mat_a);
    Complex[] fixed_pointsB = getMobiusFixedPoints(mat_b);
    Complex[] fixed_pointsC = getMobiusFixedPoints(mat_c);
    
    System.out.println("a0: " + mat_a[0]);
    System.out.println("a1: " + mat_a[1]);
    System.out.println("a2: " + mat_a[2]);
    System.out.println("a3: " + mat_a[3]);
    
    System.out.println("b0: " + mat_b[0]);
    System.out.println("b1: " + mat_b[1]);
    System.out.println("b2: " + mat_b[2]);
    System.out.println("b3: " + mat_b[3]);
    
    System.out.println("c0: " + mat_c[0]);
    System.out.println("c1: " + mat_c[1]);
    System.out.println("c2: " + mat_c[2]);
    System.out.println("c3: " + mat_c[3]);
    
    System.out.println("B fixed points 0: " + fixed_pointsB[0]);
    System.out.println("B fixed points 1: " + fixed_pointsB[1]);
    
    System.out.println("C fixed points 0: " + fixed_pointsC[0]);
    System.out.println("C fixed points 1: " + fixed_pointsC[1]);
    
    System.out.println("===================");
    mtransforms = new Object[] {mat_a, mat_b, mat_c};
  }
  
  /*
  *  Given an array m of 4 complex numbers representing parameters of a Mobius transform 
  *          f(z) = (az+b)/(cz+d)
  *     where m[0] = a, m[1] = b, m[2] = c, m[3] = d
  *     find and return the two fixed points for the Mobius transform
  */
  public Complex[] getMobiusFixedPoints(Complex[] m) {
    Complex[] fixed_points = new Complex[2];
    // solving for fixed points is solving for point where: 
    //     z = (az+b)/(cz+d), 
    //     which works out to solving quadratic equation: 
    //     cz2 + (d-a)z - b = 0, so:
    //     z = (a - d +- (sqrt((d-a)^2 + 4bc))) / 2c
    //     note that when c = 0, this yield Infinity (or NaN using Complex class)
    Complex a = m[0];
    Complex b = m[1];
    Complex c = m[2];
    Complex d = m[3];
    
    // t = sqrt((d-a)^2 + 4bc)
    Complex t = d.sub(a).power(2).add(b.mul(c).mul(4)).sqrt();
    // z+ = (a - d + t)/2c
    fixed_points[0] = a.sub(d).add(t).div(c.mul(2));
    // z- = (a - d - t)/2c
    fixed_points[1] = a.sub(d).sub(t).div(c.mul(2));
    return fixed_points;
  }
   
  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { K, L, M };
    
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_K.equalsIgnoreCase(pName)) {
      K = pValue;
    }
    else if (PARAM_L.equalsIgnoreCase(pName)) {
      L = pValue;
    }
    else if (PARAM_M.equalsIgnoreCase(pName)) {
      M = pValue;
    }
    else {
      throw new IllegalArgumentException(pName);
    }
  }

  @Override
  public String getName() {
    return "hyperptile_mobius";
  }

}
