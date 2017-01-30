package org.jwildfire.create.tina.variation;

import static org.jwildfire.base.mathlib.MathLib.sqr;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.base.Layer;

class Quaternion {
  double t;
  double x;
  double y;
  double z;
  
  public Quaternion(double t, double x, double y, double z) {
    this.t = t;
    this.x = x;
    this.y = y;
    this.z = z;
  }

}

/**
 * 
 * KleinGroup3dFunc is meant to make it easier to create "interesting" 3D Kleinian group limit sets 
 *     internally uses three 3D Mobius transformations as generators, plus their inverses
 *         (3D Mobius transformations use quaternions in place of complex numbers )
 * 
 *  Author: Gregg Helt
 */
public class KleinGroup3DFunc extends VariationFunc {
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
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    // using (modified) quaternion calcs from Mobiq by zephyrtronium http://zephyrtronium.deviantart.com/art/Mobiq-Apophysis-Plugin-170449212 converted to work in JWildfire by Brad Stefanov
   /*  
    Mobiq notes:
    Qlib uses the notation T + X i + Y j + Z k, so I used the following while
    simplifying. I took a usual Mobius transform (ax + b) / (cx + d) and made
    a, b, c, and d quaternions instead of just complex numbers, with x being
    a quaternion FTx + FTy i + FTz j + 0 k. Multiplying quaternions happens to
    be a very complex thing, and dividing is no simpler, so I needed to use
    names that are easily and quickly typed while simplifying to prevent
    massive insanity/violence/genocide. I then found switching back from those
    names in my head to be rather difficult and confusing, so I decided to
    use these macros instead.
*/

    // randomly pick one of the six calculated Mobius transformaton matrices  
    //    a, b, c, A, B, C where A = inverse(a), B = inverse(b), C = inverse(c)
    int mindex = pContext.random(6);
    Quaternion[] quat = (Quaternion[])mtransforms[mindex];
    
    // then use selected matrix for Mobius transformation:
    //  f(z) = (az + b) / (cz + d);
    Quaternion a = quat[0];
    Quaternion b = quat[1];
    Quaternion c = quat[2];
    Quaternion d = quat[3];
    Quaternion qin = new Quaternion(pAffineTP.x, pAffineTP.y, pAffineTP.z, 0);
 
    double nt = a.t * qin.t - a.x * qin.x - a.y * qin.y + b.t;
    double nx = a.t * qin.x + a.x * qin.t - a.z * qin.y + b.x;
    double ny = a.t * qin.y + a.y * qin.t + a.z * qin.x + b.y;
    double nz = a.z * qin.t + a.x * qin.y - a.y * qin.x + b.z;
    double dt = c.t * qin.t - c.x * qin.x - c.y * qin.y + d.t;
    double dx = c.t * qin.x + c.x * qin.t - c.z * qin.y + d.x;
    double dy = c.t * qin.y + c.y * qin.t + c.z * qin.x + d.y;
    double dz = c.z * qin.t + c.x * qin.y - c.y * qin.x + d.z;
    
    double normalizer = sqr(dt) + sqr(dx) + sqr(dy) + sqr(dz);

    pVarTP.x += pAmount * (nt * dt + nx * dx + ny * dy + nz * dz) / normalizer;
    pVarTP.y += pAmount * (nx * dt - nt * dx - ny * dz + nz * dy) / normalizer;
    pVarTP.z += pAmount * (ny * dt - nt * dy - nz * dx + nx * dz) / normalizer;

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
    invmat[0] = new Quaternion( m3.t,  m3.x,  m3.y,  m3.z);
    invmat[1] = new Quaternion(-m1.t, -m1.x, -m1.y, -m1.z);
    invmat[2] = new Quaternion(-m2.t, -m2.x, -m2.y, -m2.z);
    invmat[3] = new Quaternion( m0.t,  m0.x,  m0.y,  m0.z);
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
    return "klein3D_group";
  }

}
