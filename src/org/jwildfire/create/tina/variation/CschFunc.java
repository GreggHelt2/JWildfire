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

public class CschFunc extends SimpleVariationFunc {
  private static final long serialVersionUID = 1L;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    /* complex vars by cothe */
    /* exp log sin cos tan sec csc cot sinh cosh tanh sech csch coth */
    //Hyperbolic Cosecant CSCH
    double cschsin = sin(pAffineTP.y);
    double cschcos = cos(pAffineTP.y);
    double cschsinh = sinh(pAffineTP.x);
    double cschcosh = cosh(pAffineTP.x);
    double d = (cosh(2.0 * pAffineTP.x) - cos(2.0 * pAffineTP.y));
    if (d == 0) {
      return;
    }
    double cschden = 2.0 / d;
    pVarTP.x += pAmount * cschden * cschsinh * cschcos;
    pVarTP.y -= pAmount * cschden * cschcosh * cschsin;
    if (pContext.isPreserveZCoordinate()) {
  pVarTP.z += pAmount * pAffineTP.z;
}
  }

  @Override
  public String getName() {
    return "csch";
  }

}
