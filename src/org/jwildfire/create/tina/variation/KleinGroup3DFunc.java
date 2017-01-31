package org.jwildfire.create.tina.variation;

import static org.jwildfire.base.mathlib.MathLib.sqrt;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.base.Layer;

/**
 * 
 * KleinGroup3dFunc is meant to make it easier to create "interesting" 3D Kleinian group limit sets 
 *     internally uses three 3D Mobius transformations as generators, plus their inverses
 *         (3D Mobius transformations use quaternions, instead of complex numbers like 2D Mobius transformations)
 * 
 *  Author: Gregg Helt
 */
public class KleinGroup3DFunc extends VariationFunc {
  
  /** 
   * Inner class for quaternion math
   * Q(n0, n1, n2, n3) = n0 + n1*i + n2*j + n3*k
   */
  class Quaternion {
    double n0, n1, n2, n3;
    
    public Quaternion(double n0, double n1, double n2, double n3) {
      this.n0 = n0;
      this.n1 = n1;
      this.n2 = n2;
      this.n3 = n3;
    }
    
    public Quaternion add(Quaternion b) {
      return new Quaternion(n0+b.n0, n1+b.n1, n2+b.n2, n3+b.n3);
    }

    public Quaternion mult(Quaternion b) {
       // m0 = a0b0−a1b1−a2b2−a3b3
      double new0 = n0*b.n0 - n1*b.n1 - n2*b.n2 - n3*b.n3;
      //  m1 = a0b1+a1b0+a2b3−a3b2
      double new1 = n0*b.n1 + n1*b.n0 + n2*b.n3 - n3*b.n2;
      //  m2 = a0b2−a1b3+a2b0+a3b1
      double new2 = n0*b.n2 - n1*b.n3 + n2*b.n0 + n3*b.n1;
      //  m3 = a0b3+a1b2−a2b1+a3b0
      double new3 = n0*b.n3 + n1*b.n2 - n2*b.n1 + n3*b.n0;
      return new Quaternion(new0, new1, new2, new3);
    }
    
    public Quaternion reciprocal() {
      double norm = sqrt(n0*n0 + n1*n1 + n2*n2 + n3*n3);
      double normsq = norm * norm;
      // reciprocal = conjugate / (norm^2)
      return new Quaternion(n0/normsq, -n1/normsq, -n2/normsq, -n3/normsq);
    }
    
    /** 
     *  Two kinds of quaternion division, since quaternion multiplication is not commutative
     *  div1 = reciprocal(b) * this
     */
    public Quaternion div1(Quaternion b) {
      return b.reciprocal().mult(this);
    }

    /**
     *  Two kinds of quaternion division, since quaternion multiplication is not commutative
     *  div2 = this * reciprocal(b)
     */
    public Quaternion div2(Quaternion b) {
      return this.mult(b.reciprocal());
    }

  }

  private static final long serialVersionUID = 1L;

  private static final String PARAM_TR = "tr";
  private static final String PARAM_TI = "ti";
  private static final String PARAM_TJ = "tj";
  private static final String PARAM_TK = "tk";
  private static final String PARAM_BTRANS = "btrans";
  private static final String PARAM_CTRANS = "ctrans";

  private static final String[] paramNames = { PARAM_TR, PARAM_TI, PARAM_TJ, PARAM_TK, PARAM_BTRANS, PARAM_CTRANS };

  protected double tr = 2;
  protected double ti = 0;
  protected double tj = 0;
  protected double tk = 0;
  protected double btrans = 2;
  protected double ctrans = 2;
  
  protected Quaternion[] 
          mat_a, mat_a_inv, 
          mat_b, mat_b_inv, 
          mat_c, mat_c_inv;

  protected Object[] mtransforms;

  @Override
/**
 * Using Mobius transformations, but with quaternions instead of complex numbers
 */
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    
    // randomly pick one of the six calculated Mobius transformaton matrices  
    //    a, b, c, A, B, C where A = inverse(a), B = inverse(b), C = inverse(c)
    int mindex = pContext.random(6);
    Quaternion[] quat = (Quaternion[])mtransforms[mindex];
    
    // then use selected matrix for Mobius transformation:
    //  f(z) = (az + b) / (cz + d)  where a,b,c,d are quaternions
    Quaternion a = quat[0];
    Quaternion b = quat[1];
    Quaternion c = quat[2];
    Quaternion d = quat[3];
    Quaternion qin = new Quaternion(pAffineTP.x, pAffineTP.y, pAffineTP.z, 0);
 
    // use div1
    Quaternion qout = qin.mult(a).add(b).div1(qin.mult(c).add(d));
    // use div2
    // Quaternion qout2 = qin.mult(a).add(b).div2(qin.mult(c).add(d));

    pVarTP.x += pAmount * qout.n0;
    pVarTP.y += pAmount * qout.n1;
    pVarTP.z += pAmount * qout.n2;

  }
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    mat_a = new Quaternion[4];
    mat_b = new Quaternion[4];
    mat_c = new Quaternion[4];

    mat_a[0] = new Quaternion(tr, ti, tj, tk);
    //  mat_a[0] = new Quaternion(tr, ti, 0, tk);
    mat_a[1] = new Quaternion(0, -1, 0, 0);
    mat_a[2] = new Quaternion(0, -1, 0, 0);
    mat_a[3] = new Quaternion(0, 0, 0, 0);
    
    mat_b[0] = new Quaternion(1, 0, 0, 0);
    // mat_b[1] = new Quaternion(2, 0, 0, 0);
    mat_b[1] = new Quaternion(btrans, 0, 0, 0);
    mat_b[2] = new Quaternion(0, 0, 0, 0);
    mat_b[3] = new Quaternion(1, 0, 0, 0);
    
    mat_c[0] = new Quaternion(1, 0, 0, 0);
    // previous error, but gave interestingly different results 
    // mat_c[1] = new Quaternion(1, 0, 2, 0);
    // mat_c[1] = new Quaternion(0, 0, 2, 0);
    mat_c[1] = new Quaternion(0, 0, ctrans, 0);
    mat_c[2] = new Quaternion(0, 0, 0, 0);
    mat_c[3] = new Quaternion(1, 0, 0, 0);

    mat_a_inv = getInverseMatrix(mat_a);
    mat_b_inv = getInverseMatrix(mat_b);
    mat_c_inv = getInverseMatrix(mat_c);
    
    mtransforms = new Object[] {mat_a, mat_a_inv, mat_b, mat_b_inv, mat_c, mat_c_inv };
  }

  Quaternion[] getInverseMatrix(Quaternion[] mat)  {
    Quaternion[] invmat = new Quaternion[4];
    Quaternion m0 = mat[0];
    Quaternion m1 = mat[1];
    Quaternion m2 = mat[2];
    Quaternion m3 = mat[3];

    invmat[0] = new Quaternion( m3.n0,  m3.n1,  m3.n2,  m3.n3);
    invmat[1] = new Quaternion(-m1.n0, -m1.n1, -m1.n2, -m1.n3);
    invmat[2] = new Quaternion(-m2.n0, -m2.n1, -m2.n2, -m2.n3);
    invmat[3] = new Quaternion( m0.n0,  m0.n1,  m0.n2,  m0.n3);

    return invmat;
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { tr, ti, tj, tk, btrans, ctrans };
    
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_TR.equalsIgnoreCase(pName)) {
      tr = pValue;
    }
    else if (PARAM_TI.equalsIgnoreCase(pName)) {
      ti = pValue;
    }
    else if (PARAM_TJ.equalsIgnoreCase(pName)) {
      tj = pValue;
    }
    else if (PARAM_TK.equalsIgnoreCase(pName)) {
      tk = pValue;
    }
    else if (PARAM_BTRANS.equalsIgnoreCase(pName)) {
      btrans = pValue;
    }
    else if (PARAM_CTRANS.equalsIgnoreCase(pName)) {
      ctrans = pValue;
    }
    else {
      throw new IllegalArgumentException(pName);
    }
  }

  @Override
  public String getName() {
    return "klein_group_3D";
  }

}
