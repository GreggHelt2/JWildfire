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
import static org.jwildfire.base.mathlib.MathLib.sin;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

/**
 *  CircleFlowFunc
 *     by CozyG
 *  this variation opens a circular hole 
 *  by transforming coordinates to "flow" around the hole, 
 * 
 * inner_radius: radius of the hole
 * outer_radius: radius of extent of the flow effect
 *
 */
public class CircleFlowFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_OUTER_RADIUS = "outer_radius";
  private static final String PARAM_INNER_RADIUS = "inner_radius";

  private static final String[] paramNames = { PARAM_INNER_RADIUS, PARAM_OUTER_RADIUS };

  private double outer_radius = 2.0;
  private double inner_radius = 1.0;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double r = pAffineTP.getPrecalcSqrt();
    double t = pAffineTP.getPrecalcAtanYX();
    if (r < outer_radius && outer_radius != 0) {
      r = inner_radius + ((outer_radius - inner_radius)/outer_radius) * r;
    }
    double x = r * cos(t);
    double y = r * sin(t);

    pVarTP.x += pAmount * x;
    pVarTP.y += pAmount * y;

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
    return new Object[] { inner_radius, outer_radius };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_INNER_RADIUS.equalsIgnoreCase(pName))
      inner_radius = pValue;
    else if (PARAM_OUTER_RADIUS.equalsIgnoreCase(pName))
      outer_radius = pValue;
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "circle_flow";
  }


}
