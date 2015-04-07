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

import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.cosh;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.sinh;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class Sin2Func extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_XS_SHIFT = "xs_shift";
  private static final String PARAM_XC_SHIFT = "xc_shift";
  private static final String PARAM_YS_SHIFT = "ys_shift";
  private static final String PARAM_YC_SHIFT = "yc_shift";
  private static final String[] paramNames = { PARAM_XS_SHIFT, PARAM_XC_SHIFT, 
                                               PARAM_YS_SHIFT, PARAM_YC_SHIFT };

  double xs_shift = 0.0;
  double xc_shift = 0.0;
  double ys_shift = 0.0;
  double yc_shift = 0.0;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    /* Sin2Func: extending SinFunc to add user-adjustable parameters, by CozyG */
    /* complex vars by cothe */
    /* exp log sin cos tan sec csc cot sinh cosh tanh sech csch coth */
    double sinsin = sin(pAffineTP.x + xs_shift);
    double sincos = cos(pAffineTP.x + xc_shift);
    double sinsinh = sinh(pAffineTP.y + ys_shift);
    double sincosh = cosh(pAffineTP.y + yc_shift);
    pVarTP.x += pAmount * sinsin * sincosh;
    pVarTP.y += pAmount * sincos * sinsinh;
    if (pContext.isPreserveZCoordinate()) {
      pVarTP.z += pAmount * pAffineTP.z;
    }
  }

  @Override
  public String getName() {
    return "sin2";
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { xs_shift, xc_shift, 
                          ys_shift, yc_shift };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (pName.equalsIgnoreCase(PARAM_XS_SHIFT)) {
      xs_shift = pValue;
    }
    else if (pName.equalsIgnoreCase(PARAM_XC_SHIFT)) {
      xc_shift = pValue;
    }
    else if (pName.equalsIgnoreCase(PARAM_YS_SHIFT)) {
      ys_shift = pValue;
    }
    else if (pName.equalsIgnoreCase(PARAM_YC_SHIFT)) {
      yc_shift = pValue;
    }
    else
      throw new IllegalArgumentException(pName);
  }


}
