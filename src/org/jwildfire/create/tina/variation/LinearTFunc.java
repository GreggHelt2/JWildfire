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

import static org.jwildfire.base.mathlib.MathLib.fabs;
import static org.jwildfire.base.mathlib.MathLib.pow;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class LinearTFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_POWX = "powX";
  private static final String PARAM_POWY = "powY";

  private static final String[] paramNames = { PARAM_POWX, PARAM_POWY };

  private double powX = 1.2;
  private double powY = 1.2;

  private double sgn(double arg) {
    if (arg > 0)
      return 1.0;
    else
      return -1.0;
  }

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    // linearT by FractalDesire, http://fractaldesire.deviantart.com/journal/linearT-plugin-219864320
    pVarTP.x += sgn(pAffineTP.x) * pow(fabs(pAffineTP.x), this.powX) * pAmount;
    pVarTP.y += sgn(pAffineTP.y) * pow(fabs(pAffineTP.y), this.powY) * pAmount;
    pVarTP.z = pAmount * pAffineTP.z;
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { powX, powY };
  }

  @Override
  public String[] getParameterAlternativeNames() {
    return new String[] { "lT_powX", "lT_powY" };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_POWX.equalsIgnoreCase(pName))
      powX = pValue;
    else if (PARAM_POWY.equalsIgnoreCase(pName))
      powY = pValue;
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "linearT";
  }

}
