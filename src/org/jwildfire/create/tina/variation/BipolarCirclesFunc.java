/*
  JWildfire - an image and animation processor written in Java 
  Copyright (C) 1995-2011 Andreas Maschke

  This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser 
  General Public License as published by the Free Software Foundation; either version 2.1 of the 
  License, or (at your option) any later version.
 
  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License along with this software; 
  if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02110-1301 USA, or see the FSF site: http://www.fsf.org.
*/
package org.jwildfire.create.tina.variation;

import org.jwildfire.base.mathlib.Complex;
import org.jwildfire.base.mathlib.MathLib;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.cosh;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.sinh;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
/**
 * CPow revisited
 * straight implementation of complex numbers raised by complex powers
 * input x/y is considered as complex number z = x + i*y
 * params re, im are treated as complex number p = re + i*re
 * calculates z' = z^p
 * x' = real of z'
 * y' = imaginary of z'
 * 
 * 
 * Taken from Field Theory Handbook, Chapter 2, "Bipolar Circles" 
 * 
 */
public class BipolarCirclesFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_A = "a";
  private static final String PARAM_DIFF_MODE = "diff mode";
  
  private static final String[] paramNames = { PARAM_A, PARAM_DIFF_MODE };
  
  protected double a = 1;
  protected boolean diff_mode = false;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    /*   
    Complex z = new Complex(pAffineTP.x, pAffineTP.y);
    Complex p = new Complex(re, im);
    z.CPow(p);  // raise z to the p power

    pVarTP.x += pAmount * z.re;
    pVarTP.y += pAmount * z.im;
    */
    
    /*
    Complex w1 = new Complex(pAffineTP.x, pAffineTP.y);
    w1.Exp();
    Complex w2 = new Complex(w1);
    w1.Add(new Complex(1, 0));
    w2.Sub(new Complex(1,0));
    w1.Div(w2);
    w1.Mul(new Complex(a, 0));
    double x = w1.re;
    double y = w1.im;
    */
    
    double u = pAffineTP.x;
    double v = pAffineTP.y;
    double x = (a * sinh(u)) / (cosh(u) - cos(v));
    double y = (a * sin(v))  / (cosh(u) - cos(v));
  
        // Add final values in to variations totals
    if (diff_mode) {
      pVarTP.x = pAffineTP.x + (pAmount * (x - pAffineTP.x));
      pVarTP.y = pAffineTP.y + (pAmount * (y - pAffineTP.y));
    }
    else {
      pVarTP.x += pAmount * x;
      pVarTP.y += pAmount * y;
    }

    if (pContext.isPreserveZCoordinate()) {
      pVarTP.z += pAmount * pAffineTP.z;
    }
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { a, (diff_mode ? 1 : 0) };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_A.equalsIgnoreCase(pName))
      a = pValue;
    else if (PARAM_DIFF_MODE.equalsIgnoreCase(pName))
      diff_mode = (pValue >= 1);
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "bipolar_circles";
  }

}
