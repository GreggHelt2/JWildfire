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

/**
  Epitrochoid
  Implemented by CozyG, April 2015
  For reference see http://en.wikipedia.org/wiki/Epitrochoid

  Interesting relationships:

  Loops always cross origin (center of circle A?) at c_scale = cusps + 1; 
       so when c_radius = a_radius + (1/b_radius)
  Inner loops touch at c_scale ~= 4.5 fro larger # of cusps, and down to a limit of 3 for cusps = 2, ~3.67 for cusps = 3...
          this is independent of b_radius
          I think there's probably a way to figure this out into an equation, some factor that weights each succeeding cusp less?
          (c_radius = 4.5 * b_radius)

    a_radius = b_radius * cusps;
    c_radius = b_radius * c_scale;

  c_scale = c_radius / b_radius
  cusps = a_radius / b_radius
  c_radius / b_radius = (a_radius / b_radius) + 1
  c_radius = a_radius + (1/b_radius)
*/
public class EpitrochoidFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  //  private static final String PARAM_A_RADIUS = "a_radius";
  private static final String PARAM_CUSPS = "cusps";  // a_radius/b_radius (so cusps and b_radius determines a_radius
  private static final String PARAM_B_RADIUS = "b_radius";
  private static final String PARAM_C_SCALE = "c_scale";   // c_radius (distance from point to center of circle B) is determined by b_radius and c_scale
  private static final String PARAM_UNIFIED_INNER_OUTER = "unified_inner_outer";
  private static final String PARAM_INNER_MODE = "inner_mode";
  private static final String PARAM_OUTER_MODE = "outer_mode";
  private static final String PARAM_INNER_SPREAD = "inner_spread";
  private static final String PARAM_OUTER_SPREAD = "outer_spread";
  private static final String PARAM_INNER_SPREAD_RATIO = "inner_spread_ratio";
  private static final String PARAM_OUTER_SPREAD_RATIO = "outer_spread_ratio";
  private static final String PARAM_SPREAD_SPLIT = "spread_split";
  private static final String PARAM_FILL = "fill";

  // my standard approach if there is a cycles variable is that if cycles is set to 0, 
  //     that means function should decide cycle value automatically
  //     for curves, function will make best effort to calculate minimum number of cycles needed 
  //     to close curve, or a somewhat arbitrary number if cannot 
  private static final String PARAM_CYCLES = "cycles";


  private static final String[] paramNames = { PARAM_CUSPS, PARAM_B_RADIUS, PARAM_C_SCALE,
                                               PARAM_UNIFIED_INNER_OUTER, PARAM_INNER_MODE, PARAM_OUTER_MODE, 
                                               PARAM_INNER_SPREAD, PARAM_OUTER_SPREAD, 
                                               PARAM_INNER_SPREAD_RATIO, PARAM_OUTER_SPREAD_RATIO, PARAM_SPREAD_SPLIT,
                                               PARAM_CYCLES, PARAM_FILL }; 

  private double cyclesParam = 0;  // number of cycles (2*PI radians, circle circumference), if set to 0 then number of cycles is calculated automatically
  // private double a_radius = 1;   // radius of circle "A" (the stationary circle)

  private double b_radius = 0.2;   // radius of circle "B" (circle that is rolling around outside of circle "A")
  private double a_radius; // radius of circle "A", dermined by b_radius and cusps params
  private double c_radius; // "radius" of point "C" (fixed distnace from center of circle "B" to point "C")
  private double cusps = 5;

  private double c_scale = 1.5; 
  private int unified_inner_outer = 1;
  private int inner_mode = 1;
  private int outer_mode = 1;
  private double inner_spread = 0; // deform based on original x/y
  private double outer_spread = 0; // deform based on original x/y
  private double inner_spread_ratio = 1; // how much outer_spread applies to x relative to y
  private double outer_spread_ratio = 1; // how much outer_spread applies to x relative to y
  private double spread_split = 1;
  private double fill = 0;
  private double cycles;  // 1 cycle = 2*PI

  private double cycle_length = 2 * M_PI; // 2(PI)
  private double radians_to_close = 2 * M_PI;  
  private double cycles_to_close = radians_to_close / cycle_length;  // = PI^2

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    // setting cyclesParam to 0 is meant to indicate the function should determine how many cycles 
    //     (actually setting cycles = 0 will just yield a single point)
    a_radius = b_radius * cusps;
    c_radius = b_radius * c_scale;
    if (cyclesParam == 0) {  
      cycles = cycles_to_close;
    }
    else {
      cycles = cyclesParam;
    }
  }

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double tin = atan2(pAffineTP.y, pAffineTP.x);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    double theta = cycles * tin;
    double rin = spread_split * sqrt((pAffineTP.x  * pAffineTP.x) + (pAffineTP.y * pAffineTP.y));
        
    double x = ((a_radius + b_radius) * cos(theta)) - (c_radius * cos(((a_radius + b_radius)/b_radius) * theta));
    double y = ((a_radius + b_radius) * sin(theta)) - (c_radius * sin(((a_radius + b_radius)/b_radius) * theta));
    double r = sqrt(x*x + y*y);
    double t = atan2(y, x);

    if (fill != 0) { 
      r = r + (fill * (pContext.random() - 0.5));
    }

    double xin, yin;
    double rinx, riny;

    if ((abs(rin) > abs(r)) || (unified_inner_outer == 1))  { // incoming point lies "outside" of curve OR ignoring inner_outer distinction
      switch(outer_mode) {
        case 0: // no spread
          pVarTP.x += pAmount * x;
          pVarTP.y += pAmount * y;
          break;
        case 1:
          rinx = (rin * outer_spread * outer_spread_ratio) - (outer_spread * outer_spread_ratio) + 1;
          riny = (rin * outer_spread) - outer_spread + 1;
          pVarTP.x += pAmount * rinx * x;
          pVarTP.y += pAmount * riny * y;
          break;
        case 2:
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (x + (outer_spread * outer_spread_ratio * (xin-x)));
          pVarTP.y += pAmount * (y + (outer_spread * (yin-y)));
          break;
        case 3:
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (x + (outer_spread * outer_spread_ratio * xin));
          pVarTP.y += pAmount * (y + (outer_spread * yin));
          break;
        case 4:
          rinx = (0.5 * rin) + (outer_spread * outer_spread_ratio);
          riny = (0.5 * rin) + outer_spread;
          pVarTP.x += pAmount * rinx * x;
          pVarTP.y += pAmount * riny * y;
          break;
        case 5: // same as outer_mode 3, but without the sign modifications
          pVarTP.x += pAmount * (x + (outer_spread * outer_spread_ratio * pAffineTP.x));
          pVarTP.y += pAmount * (y + (outer_spread * pAffineTP.y));
          break;
        case 6: 
          pVarTP.x += pAmount * pAffineTP.x;
          pVarTP.y += pAmount * pAffineTP.y;
          break;
        case 8:
          pVarTP.x += pAmount * rin * cos(t) * outer_spread * outer_spread_ratio;
          pVarTP.y += pAmount * rin * sin(t) * outer_spread;
          break;
        default:
          pVarTP.x += pAmount * x;
          pVarTP.y += pAmount * y;
          break;
      }
    }
    else  { // incoming point lies "inside" or "on" curve
      switch(inner_mode) {
        case 0: // no spread
          pVarTP.x += pAmount * x;
          pVarTP.y += pAmount * y;
          break;
        case 1:
          rinx = (rin * inner_spread * inner_spread_ratio) - (inner_spread * inner_spread_ratio) + 1;
          riny = (rin * inner_spread) - inner_spread + 1;
          pVarTP.x += pAmount * rinx * x;
          pVarTP.y += pAmount * riny * y;
          break;
        case 2:
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (x - (inner_spread * inner_spread_ratio * (x-xin)));
          pVarTP.y += pAmount * (y - (inner_spread * (y-yin)));
          break;
        case 3:
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (x - (inner_spread * inner_spread_ratio * xin));
          pVarTP.y += pAmount * (y - (inner_spread * yin));
          break;
        case 4:
          rinx = (0.5 * rin) + (inner_spread * inner_spread_ratio);
          riny = (0.5 * rin) + inner_spread;
          pVarTP.x += pAmount * rinx * x;
          pVarTP.y += pAmount * riny * y;
          break;
        case 5: // same as inner_mode 3, but without the sign modifications
          pVarTP.x += pAmount * (x + (inner_spread * inner_spread_ratio * pAffineTP.x));
          pVarTP.y += pAmount * (y + (inner_spread * pAffineTP.y));
          break;
        case 6: 
          pVarTP.x += pAmount * pAffineTP.x;
          pVarTP.y += pAmount * pAffineTP.y;
          break;
        case 8:
          pVarTP.x += pAmount * rin * cos(t) * (inner_spread * inner_spread_ratio);
          pVarTP.y += pAmount * rin * sin(t) * inner_spread;
          break;
        default:
          pVarTP.x += pAmount * x;
          pVarTP.y += pAmount * y;
          break;
      }
    }
    pVarTP.z += pAmount * pAffineTP.z;
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { cusps, b_radius, c_scale,
                          unified_inner_outer, inner_mode, outer_mode, inner_spread, outer_spread,
                          inner_spread_ratio, outer_spread_ratio, spread_split,
                          cyclesParam, fill };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_CUSPS.equalsIgnoreCase(pName))
      cusps = pValue;
    else if (PARAM_B_RADIUS.equalsIgnoreCase(pName))
      b_radius = pValue;
    else if (PARAM_C_SCALE.equalsIgnoreCase(pName))
      c_scale = pValue;
    else if (PARAM_UNIFIED_INNER_OUTER.equalsIgnoreCase(pName)) {
      unified_inner_outer = (pValue == 0 ? 0 : 1);
    }
    else if (PARAM_OUTER_MODE.equalsIgnoreCase(pName)) {
      outer_mode = (int)floor(pValue);
      if (outer_mode > 10 || outer_mode < 0) { outer_mode = 0; }
    }
    else if (PARAM_INNER_MODE.equalsIgnoreCase(pName)) {
      inner_mode = (int)floor(pValue);
      if (inner_mode > 10 || inner_mode < 0) { inner_mode = 0; }
    }
    else if (PARAM_OUTER_SPREAD.equalsIgnoreCase(pName))
      outer_spread = pValue;
    else if (PARAM_INNER_SPREAD.equalsIgnoreCase(pName))
      inner_spread = pValue;
    else if (PARAM_OUTER_SPREAD_RATIO.equalsIgnoreCase(pName))
      outer_spread_ratio = pValue;
    else if (PARAM_INNER_SPREAD_RATIO.equalsIgnoreCase(pName))
      inner_spread_ratio = pValue;    
    else if (PARAM_SPREAD_SPLIT.equalsIgnoreCase(pName))
      spread_split = pValue;
    else if (PARAM_CYCLES.equalsIgnoreCase(pName))
      cyclesParam = pValue;
    else if (PARAM_FILL.equalsIgnoreCase(pName))
      fill = pValue;
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "epitrochoid";
  }

}
