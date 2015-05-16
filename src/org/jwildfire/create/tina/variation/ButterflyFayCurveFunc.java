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

import static java.lang.Math.abs;
import java.math.BigInteger;
import java.util.ArrayList;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.atan2;
import static org.jwildfire.base.mathlib.MathLib.exp;
import static org.jwildfire.base.mathlib.MathLib.pow;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.floor;
import static org.jwildfire.base.mathlib.MathLib.sqrt;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class ButterflyFayCurveFunc extends AbstractPolarCurveFunc {
  private static final long serialVersionUID = 1L;
    
  private static final String PARAM_RADIAL_OFFSET = "radial_offset";
  private static final String[] additionalParamNames = { PARAM_RADIAL_OFFSET };

  private double radial_offset = 0;
  //  private double radians_to_close = 2 * M_PI * M_PI * M_PI; // 2(PI)^3
  // private double cycles_to_close = radians_to_close / cycle_length; // = PI^2
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    super.init(pContext, pLayer, pXForm, pAmount);  // calls recalcCycles() and recalcCurveIntersects()
  }

  /* called from super.init() */
  @Override
  public void recalcCycles() {
    double radians_to_close = 2 * M_PI * M_PI * M_PI;   // radians to close is 2 * (PI^3)
    cycles_to_close = radians_to_close / (2 * M_PI);     // so cycles to close is PI^2
    if (cycles_param == 0) {
      // cycles = cycles_to_close * metacycles; // moved metacycle calcs to superclass
      cycles = cycles_to_close;
    }
    else {
      cycles = cycles_param;
    }
    super.recalcCycles();
  }

  /* called from super.transform() */
  public void calcCurvePoint(FlameTransformationContext pContext, double theta, XYZPoint pResult) {
    // pResult should be zero'd out before getting here
    double r = 0.5 * (exp(cos(theta)) - (2 * cos(4 * theta)) - pow(sin(theta / 12), 5) + radial_offset);
    pResult.x = r * sin(theta);
    pResult.y = -1 * r * cos(theta); // multiplying by -1 to flip butterfly to point up
    pResult.x *= curve_scale;
    pResult.y *= curve_scale;
    // z unchanged?
    super.calcCurvePoint(pContext, theta, pResult);
  }
  
  @Override
  public String[] getParameterNames() {
    return joinArrays(additionalParamNames, paramNames);
  }

  @Override
  public Object[] getParameterValues() {
    return joinArrays(new Object[] { radial_offset }, super.getParameterValues());
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_RADIAL_OFFSET.equalsIgnoreCase(pName))
      radial_offset = pValue;
    else {
      super.setParameter(pName, pValue);
    }
  }

  @Override
  public String getName() {
    return "butterfly_fay_curve";
  }

}
