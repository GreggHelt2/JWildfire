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

import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sin;
import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class CircleCurveFunc extends AbstractPolarCurveFunc {
  private static final long serialVersionUID = 1L;
  protected boolean randomize = true;

  // private double radial_offset = 0;
  //  private double radians_to_close = 2 * M_PI * M_PI * M_PI; // 2(PI)^3
  // private double cycles_to_close = radians_to_close / cycle_length; // = PI^2
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    super.init(pContext, pLayer, pXForm, pAmount);  // calls recalcCycles() and recalcCurveIntersects()
  }

  /* called from super.init() */
  @Override
  public void recalcCycles() {
    cycles_to_close = 1;  
    if (cycles_param == 0) {
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
    double r = curve_scale;
    if (randomize) {
      // totally radomize theta to force curve:
      theta = (pContext.random() * M_2PI) - M_PI;
    }
    // if (theta < 0) { theta = theta + M_2PI; }
    // double r = (M_PI + theta) / M_2PI;
    pResult.x = r * cos(theta);
    pResult.y = r * sin(theta); 
    // if (pResult.y < 0) { pResult.y += -1; }
   // pResult.x *= curve_scale;
   // pResult.y *= curve_scale;
    // z unchanged?
    super.calcCurvePoint(pContext, theta, pResult);
  }
  
  /*  @Override
  public String[] getParameterNames() {
    return 
    return joinArrays(additionalParamNames, paramNames);
  }
  */

  /*
  @Override
  public Object[] getParameterValues() {
    return joinArrays(new Object[] { radial_offset }, super.getParameterValues());
  }
  */

  /*
  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_RADIAL_OFFSET.equalsIgnoreCase(pName))
      radial_offset = pValue;
    else {
      super.setParameter(pName, pValue);
    }
  }
  */

  @Override
  public String getName() {
    return "circle_curve";
  }

}
