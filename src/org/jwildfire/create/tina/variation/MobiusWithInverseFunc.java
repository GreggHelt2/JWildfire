
package org.jwildfire.create.tina.variation;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.nfunk.jep.type.Complex;

/**
 * MobisuWithInverseFunc uses parameters to create a Mobius transformation M:
 *    M(z) = (az + b) / (cz + d) where a,b,c,d are complex numbers
 * and also creates the inverse Mobius transormation M'
 *    such that M'(M(z)) = z
 * when transform() called, randomly chooses between M and M' for transformation 
 * 
 * Author: Gregg Helt
 */
public class MobiusWithInverseFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_A_RE = "a_re";
  private static final String PARAM_A_IM = "a_im";
  private static final String PARAM_B_RE = "b_re";
  private static final String PARAM_B_IM = "b_im";
  private static final String PARAM_C_RE = "c_re";
  private static final String PARAM_C_IM = "c_im";
  private static final String PARAM_D_RE = "d_re";
  private static final String PARAM_D_IM = "d_im";


  private static final String[] paramNames = { 
    PARAM_A_RE, PARAM_A_IM, 
    PARAM_B_RE, PARAM_B_IM, 
    PARAM_C_RE, PARAM_C_IM, 
    PARAM_D_RE, PARAM_D_IM 
  };

  private double a_re = 1;
  private double a_im = 1;
  private double b_re = 0;
  private double b_im = 0;
  private double c_re = 1;
  private double c_im = 0;
  private double d_re = 1;
  private double d_im = 0;
  
  private Complex[] mobius;
  private Complex[] mobius_inverse;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {

    // randomly select either the Mobius transformation or its inverse
    Complex[] mat;
    double rand = pContext.random();
    if (rand < 0.5) {
      mat = mobius;
    }
    else {
      mat = mobius_inverse;
    }
    // then use selected matrix for Mobius transformation:
    //    f(z) = (az + b) / (cz + d);
    // for the generator matrices 
    //    [0, 1, 2, 3] = [a, b, c, d]  ==> f(z)= (az+b)/(cz+d)
    //     
    Complex zin = new Complex(pAffineTP.x, pAffineTP.y);
    Complex a = mat[0];
    Complex b = mat[1];
    Complex c = mat[2];
    Complex d = mat[3];
    Complex zout = zin.mul(a).add(b).div(zin.mul(c).add(d));
    
    pVarTP.x += pAmount * zout.re();
    pVarTP.y += pAmount * zout.im();

    if (pContext.isPreserveZCoordinate()) {
      pVarTP.z += pAmount * pAffineTP.z;
    }
  }
  
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    // normalize mobius matrix
    Complex[] mobius_prenorm = new Complex[4];
    mobius_prenorm[0] = new Complex(a_re, a_im);
    mobius_prenorm[1] = new Complex(b_re, b_im);
    mobius_prenorm[2] = new Complex(c_re, c_im);
    mobius_prenorm[3] = new Complex(d_re, d_im);
    mobius = normalize(mobius_prenorm);
    // create inverse matrix
    mobius_inverse = matrixInverse(mobius);
  }
  
  /* assumes 2x2 matrix represented by 4-element array: [a b c d] */
  public Complex[] normalize(Complex[] mat) {
    Complex[] normalized = new Complex[4];
    Complex a = mat[0];
    Complex b = mat[1];
    Complex c = mat[2];
    Complex d = mat[3];
    // determinant = ad - bc
    Complex det = a.mul(d).sub(b.mul(c));
    Complex denom = det.sqrt();
    normalized[0] = a.div(denom);
    normalized[1] = b.div(denom);
    normalized[2] = c.div(denom);
    normalized[3] = d.div(denom);
    return normalized;
  }
  
  /* assumes 2x2 matrix represented by 4-element array: [a b c d] */
  public Complex[] matrixInverse(Complex[] mat) {
    Complex[] matinv = new Complex[4];
    matinv[0] = mat[3];
    matinv[1] = mat[1].mul(-1);
    matinv[2] = mat[2].mul(-1);
    matinv[3] = mat[0];
    return matinv;
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { a_re, a_im, b_re, b_im, c_re, c_im, d_re, d_im };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_A_RE.equalsIgnoreCase(pName)) {
      a_re = pValue;
    }
    else if (PARAM_A_IM.equalsIgnoreCase(pName)) {
      a_im = pValue;
    }
    else if (PARAM_B_RE.equalsIgnoreCase(pName)) {
      b_re = pValue;
    }
    else if (PARAM_B_IM.equalsIgnoreCase(pName)) {
      b_im = pValue;
    }
    else if (PARAM_C_RE.equalsIgnoreCase(pName)) {
      c_re = pValue;
    }
    else if (PARAM_C_IM.equalsIgnoreCase(pName)) {
      c_im = pValue;
    }
    else if (PARAM_D_RE.equalsIgnoreCase(pName)) {
      d_re = pValue;
    }
    else if (PARAM_D_IM.equalsIgnoreCase(pName)) {
      d_im= pValue;
    }
    else {
      throw new IllegalArgumentException(pName);
    }
  }

  @Override
  public String getName() {
    return "mobius_with_inverse";
  }

}