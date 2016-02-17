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

import static java.lang.Math.abs;

/**
 *  Hearts4 variation function
 */
public class Hearts4SimpleFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_A = "a"; 
  private static final String PARAM_B = "b"; 
  private static final String PARAM_C = "c"; 
  private static final String PARAM_D = "d"; 
  private static final String PARAM_INNER_MODE = "inner_mode";
  private static final String PARAM_OUTER_MODE = "outer_mode";
  private static final String PARAM_INNER_SPREAD = "inner_spread";
  private static final String PARAM_OUTER_SPREAD = "outer_spread";
  private static final String PARAM_INNER_SPREAD_RATIO = "inner_spread_ratio";
  private static final String PARAM_OUTER_SPREAD_RATIO = "outer_spread_ratio";
  private static final String PARAM_SPREAD_SPLIT = "spread_split";
  private static final String PARAM_THICKNESS = "thickness";
  private static final String PARAM_DIFF_MODE = "diff_mode";

  private static final String[] paramNames = { PARAM_A, PARAM_B, PARAM_C, PARAM_D, 
                                               PARAM_INNER_MODE, PARAM_OUTER_MODE, 
                                               PARAM_INNER_SPREAD, PARAM_OUTER_SPREAD, 
                                               PARAM_INNER_SPREAD_RATIO, PARAM_OUTER_SPREAD_RATIO, PARAM_SPREAD_SPLIT, PARAM_THICKNESS, PARAM_DIFF_MODE };
  private static int MIN_MODE = 0;
  private static int MAX_MODE = 8;
  private static int DEFAULT_MODE = MIN_MODE;

  private double cyclesParam = 0;  // number of cycles (2*PI radians, circle circumference), if set to 0 then number of cycles is calculated automatically

  private double a = 1;
  private double b = 1;
  private double c = 1;
  private double d = 1;
  
  private int inner_mode = DEFAULT_MODE;
  private int outer_mode = DEFAULT_MODE;
  private double inner_spread = 0.5; // deform based on original x/y
  private double outer_spread = 0.5; // deform based on original x/y
  private double inner_spread_ratio = 1; // how much inner_spread applies to x relative to y
  private double outer_spread_ratio = 1; // how much outer_spread applies to x relative to y
  private double spread_split = 1;
  private double thickness = 0; // amount to thicken curve by randomizing input

  private double cycles;  // 1 cycle = 2*PI
  private double radians; // = 2*PI*cycles

  private double cycle_length = 2 * M_PI; // 2(PI)
  private double cycles_to_close = 0; // 0 indicates unknown, -1 indicates curve will never close

  private boolean DEBUG = false;
  private boolean diff_mode = false;

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    cyclesParam = d;
    if (cyclesParam == 0) {  
        cycles_to_close = 1;
        cycles = cycles_to_close;
    }
    else {
      cycles = cyclesParam;
    }
    radians = cycles * 2 * M_PI;
  }
  
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double xin = pAffineTP.x;
    double yin = pAffineTP.y;
    double tin = atan2(yin, xin);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    double theta = cycles * tin;
    double raw_rin = sqrt((xin * xin) + (yin * yin));
    double rin = spread_split * raw_rin;

// HEARTS2
// x = 16sin^3(t)
// y = 13cost-5cos(2t)-2cos(3t)-cos(4t)
// double sign = (sin(theta) >= 0 ? -1 : 1);
// double x = -0.08 * sign * (16 * Math.pow(sin(theta), radius * 3));
double x = -0.08 * (16 * Math.pow(sin(theta), a * 3));
double y = -0.08 * ((13 * cos(b * theta)) - (5 * cos(2*theta)) - (2 * cos(3*theta)) - cos(4*theta));
double r = sqrt(x*x + y*y);

double t = atan2(y, x);

    double rinx, riny;
    double xout, yout, rout;

    if (abs(rin) > abs(r))  { // incoming point lies "outside" of curve
      switch(outer_mode) {
        case 0: // no spread
          if (thickness == 0) {
            xout = x;
            yout = y;
          }
          else {
            rout = r + (thickness * (pContext.random() - 0.5));
            xout = rout * cos(t);
            yout = rout * sin(t);
          }
          break;
        case 1:
          rinx = (rin * outer_spread * outer_spread_ratio) - (outer_spread * outer_spread_ratio) + 1;
          riny = (rin * outer_spread) - outer_spread + 1;
          xout = rinx * x;
          yout = riny * y;
          break;
        case 2:
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          xout = (x + (outer_spread * outer_spread_ratio * (xin-x)));
          yout = (y + (outer_spread * (yin-y)));
          break;
        case 3:
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          xout = (x + (outer_spread * outer_spread_ratio * xin));
          yout = (y + (outer_spread * yin));
          break;
        case 4:
          rinx = (0.5 * rin) + (outer_spread * outer_spread_ratio);
          riny = (0.5 * rin) + outer_spread;
          xout = rinx * x;
          yout = riny * y;
          break;
        case 5: // same as outer_mode 3, but without the sign modifications
          xout = (x + (outer_spread * outer_spread_ratio * xin));
          yout = (y + (outer_spread * yin));
          break;
        case 6: 
          xout = xin;
          yout = yin;
          break;
        case 7:
          xout = rin * cos(t) * outer_spread * outer_spread_ratio;
          yout = rin * sin(t) * outer_spread;
          break;
        case 8: 
          rout = r + ((rin-r) * outer_spread);
          xout = rout * cos(t);
          yout = rout * sin(t);
          break;
        default:
          xout = x;
          yout = y;
          break;
      }
    }
    else  { // incoming point lies "inside" or "on" curve
      switch(inner_mode) {
        case 0: // no spread
          if (thickness == 0) {
            xout = x;
            yout = y;
          }
          else {
            rout = r + (thickness * (pContext.random() - 0.5));
            xout = rout * cos(t);
            yout = rout * sin(t);
          }
          break;
        case 1:
          rinx = (rin * inner_spread * inner_spread_ratio) - (inner_spread * inner_spread_ratio) + 1;
          riny = (rin * inner_spread) - inner_spread + 1;
          xout = rinx * x;
          yout = riny * y;
          break;
        case 2:
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          xout = (x - (inner_spread * inner_spread_ratio * (x-xin)));
          yout = (y - (inner_spread * (y-yin)));
          break;
        case 3:
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          xout = (x - (inner_spread * inner_spread_ratio * xin));
          yout = (y - (inner_spread * yin));
          break;
        case 4:
          rinx = (0.5 * rin) + (inner_spread * inner_spread_ratio);
          riny = (0.5 * rin) + inner_spread;
          xout = rinx * x;
          yout = riny * y;
          break;
        case 5: // same as inner_mode 3, but without the sign modifications
          xout = (x + (inner_spread * inner_spread_ratio * xin));
          yout = (y + (inner_spread * yin));
          break;
        case 6: 
          xout = xin;
          yout = yin;
          break;
        case 7:
          xout = rin * cos(t) * (inner_spread * inner_spread_ratio);
          yout = rin * sin(t) * inner_spread;
          break;
        case 8: 
          rout = r + ((rin - r) * inner_spread);
          xout = rout * cos(t);
          yout = rout * sin(t);
          break;
        default:
          xout = x;
          yout = y;
          break;
      }
    }
    // Add final values in to variations totals
    if (diff_mode) {
      pVarTP.x = pAffineTP.x + (pAmount * (xout - pAffineTP.x));
      pVarTP.y = pAffineTP.y + (pAmount * (yout - pAffineTP.y));
    }
    else {
      pVarTP.x += pAmount * xout;
      pVarTP.y += pAmount * yout;
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
    return new Object[] { a, b, c, d, 
                          inner_mode, outer_mode, inner_spread, outer_spread,
                          inner_spread_ratio, outer_spread_ratio, spread_split, thickness, (diff_mode ? 1 : 0)  };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_A.equalsIgnoreCase(pName))
      a = pValue;
    else if (PARAM_B.equalsIgnoreCase(pName))
      b = pValue;
    else if (PARAM_C.equalsIgnoreCase(pName))
      c = pValue;
    else if (PARAM_D.equalsIgnoreCase(pName))
      d = pValue;
    else if (PARAM_OUTER_MODE.equalsIgnoreCase(pName)) {
      outer_mode = (int)floor(pValue);
      if (outer_mode < MIN_MODE || outer_mode > MAX_MODE) { outer_mode = DEFAULT_MODE; }
    }
    else if (PARAM_INNER_MODE.equalsIgnoreCase(pName)) {
      inner_mode = (int)floor(pValue);
      if (inner_mode < MIN_MODE || inner_mode > MAX_MODE) { inner_mode = DEFAULT_MODE; }
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
    else if (PARAM_THICKNESS.equalsIgnoreCase(pName) || pName.equalsIgnoreCase("fill"))
      thickness = pValue;
    else if (PARAM_DIFF_MODE.equalsIgnoreCase(pName) || pName.equalsIgnoreCase("diff mode")) {
      diff_mode = (pValue >= 1);
    }
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "hearts4_simple";
  }

}
