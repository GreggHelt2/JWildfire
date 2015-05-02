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

import static org.jwildfire.base.mathlib.MathLib.EPSILON;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.exp;
import static org.jwildfire.base.mathlib.MathLib.sin;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class FociFunc extends SimpleVariationFunc {
  private static final long serialVersionUID = 1L;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    /* Foci in the Apophysis Plugin Pack */

    double expx = exp(pAffineTP.x) * 0.5;
    double expnx = 0.25 / expx;
    if (expx <= EPSILON || expnx <= EPSILON) {
      return;
    }
    double siny = sin(pAffineTP.y);
    double cosy = cos(pAffineTP.y);

    double tmp = (expx + expnx - cosy);
    if (tmp == 0)
      return;
    tmp = pAmount / tmp;

    pVarTP.x += (expx - expnx) * tmp;
    pVarTP.y += siny * tmp;
    pVarTP.z = pAmount * pAffineTP.z;
  }

  @Override
  public String getName() {
    return "foci";
  }

}
