/*
  JWildfire - an image and animation processor written in Java 
  Copyright (C) 1995-2015 Andreas Maschke

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
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sin;

import org.jwildfire.base.Tools;
import static org.jwildfire.base.mathlib.MathLib.atan2;
import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class SymmetricDrosteFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_CENTRE_X = "centre_x";
  private static final String PARAM_CENTRE_Y = "centre_y";
  private static final String PARAM_SCALE = "scale";
  private static final String PARAM_ORDER = "order";
  private static final String PARAM_ROTATION = "rotation_degrees";
  private static final String PARAM_PASSTHROUGH = "pass_through";
  private static final String[] paramNames = { PARAM_CENTRE_X, PARAM_CENTRE_Y, PARAM_ROTATION, PARAM_ORDER, PARAM_PASSTHROUGH };

  private double centre_x = 0;
  private double centre_y = 0;
  private int order = 3;
  double pass_through = 0.5;
  double scale = 0.5;
  double rotation = 0;  // rotation in degrees
  double theta;  // rotation in radians
  
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double rnd = pContext.random();
    double xin, yin, tin;
    xin = pAffineTP.x;
    yin = pAffineTP.y;
    tin = atan2(yin, xin); // atan2 range is [-Pi..+Pi]
    
    double xout, yout;
    if (rnd < pass_through) {
      pVarTP.x += pAffineTP.x;
      pVarTP.y += pAffineTP.y;
      pVarTP.z += pAffineTP.z;
    }
    else {
      // double dx = (pVarTP.x - centre_x) * pAmount;
      // double dy = (pVarTP.y - centre_y) * pAmount;
      // double dx = (pAffineTP.x - centre_x) * pAmount;
      // double dy = (pAffineTP.y - centre_y) * pAmount
      
      // rotate (in degrees)
      double xrot = xin * cos(theta) - yin * sin(theta);
      double yrot = xin * sin(theta) + yin * cos(theta);
      
      // scale
      double xscaled = xrot * pAmount;
      double yscaled = yrot* pAmount;


      int index = (int)(Math.random()*(order+1));
      /*
      double xrot = (xscaled * _cosa[idx]) + (yscaled * _sina[idx]);
      double yrot = (xscaled * _sina[idx]) + (yscaled * _cosa[idx]);
//      pVarTP.x += centre_x + dx * _cosa[idx] + dy * _sina[idx];
//      pVarTP.y += centre_y + dy * _cosa[idx] - dx * _sina[idx];
      */
      // rotate result copies
      double tout = ((double)index / (double)order) * M_2PI;
      double offset = centre_x;
      // translate
      xout = xscaled + (offset * cos(tout));
      yout = yscaled + (offset * sin(tout));
      pVarTP.x += xout;
      pVarTP.y += yout;
      pVarTP.z += pAmount * pAffineTP.z;
    }
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { centre_x, centre_y, rotation, order, pass_through };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_CENTRE_X.equalsIgnoreCase(pName)) {
      centre_x = pValue;
    }
    else if (PARAM_CENTRE_Y.equalsIgnoreCase(pName)) {
      centre_y = pValue;
    }
    else if (PARAM_ROTATION.equalsIgnoreCase(pName)) {
      rotation = pValue;
    }
    else if (PARAM_ORDER.equalsIgnoreCase(pName)) {
      order = limitIntVal(Tools.FTOI(pValue), 1, Integer.MAX_VALUE);
    }

    else if (PARAM_PASSTHROUGH.equalsIgnoreCase(pName)) {
      pass_through = pValue;
    }
    else {
      throw new IllegalArgumentException(pName);
    }
  }

  @Override
  public String getName() {
    return "symmetric_droste";
  }

  @Override
  public int getPriority() {
    return 1;
  }

  private double _sina[], _cosa[];

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    theta = rotation * M_2PI / 360;
    _sina = new double[order];
    _cosa = new double[order];
    double da = M_2PI / (double) order;
    double angle = 0.0;
    for (int i = 0; i < order; i++) {
      _sina[i] = sin(angle);
      _cosa[i] = cos(angle);
      angle += da;
    }
  }

}
