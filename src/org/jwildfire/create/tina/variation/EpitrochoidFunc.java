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

/**
  Epitrochoid
  Implemented by CozyG, April 2015
  For reference see http://en.wikipedia.org/wiki/Epitrochoid

  Interesting relationships:
  If c_radius = b_radius, curve is an epicycloid (no retrograde cycles, just point cusps)
  If a_radius = b_radius, curve is a limacon
  If a_radius = b_radius = c_radius, curve is a cardioid

  Loops always cross origin (center of circle A?) at c_scale = cusps + 1; 
       so when c_radius = a_radius + (1/b_radius)
  Inner loops touch at c_scale ~= 4.5 fro larger # of cusps, and down to a limit of 3 for cusps = 2, ~3.67 for cusps = 3...
          this is independent of b_radius
          I think there's probably a way to figure this out into an equation, some factor that weights each succeeding cusp less?
          (c_radius = 4.5 * b_radius)

 Figuring out better parameters that map to obvious features:
  if a = radius of stationary circle
     b = radius of rolling circle
     c = radius of point attached to circle b

  then max radius of the epitrochoid will be a+b+c
       min radius of the epitrochoid will be a+b-c
       "average" radius of the epitrochoid will be a+b
  
  cusps = k = a/b is the number of retrograde cycles in the closed curve
  cscale = s = c/b is the radius of point C (relative to origin of circle B)
  cusp_size = cscale-1 is cscale adjusted, so that at cusp_size == 0 the cusps terminate in a point instead of having a retrograde cycle

  so make ajustable values be:
       radius = r = a+b+c
       cusps  = k = a/b
       cscale = s = c/b
  
   so 
   given r, k, and c:
       r = a + b + c  ==> b = r - a - c
       k = a/b        ==> a = bk
       s = c/b        ==> c = bs
       substituting:
           b = r - bk - bs
           b = r -b(k+s)
           b + b(k+s) = r
           b(1 + k + s) = r
           b = r/(k + s + 1)  
           solve for b, then a = bk, c = bs
*/
public class EpitrochoidFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  // cusps, radius, and cusp_size are used to determine:
  //    a_radius: radius of stationary circle
  //    b_radius: radius of rolling circle
  //    c_radius: radius of rolling point attached to center of rolling circle
  private static final String PARAM_CUSPS = "cusps";  //  = a_radius/b_radius 
  private static final String PARAM_CUSP_SIZE = "cusp_size";  //  = (c_radius/b_radius) - 1
  private static final String PARAM_RADIUS = "radius";  // = a_radius + b_radius + c_radius
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


  private static final String[] paramNames = { PARAM_CUSPS, PARAM_CUSP_SIZE, PARAM_RADIUS, 
                                               PARAM_UNIFIED_INNER_OUTER, PARAM_INNER_MODE, PARAM_OUTER_MODE, 
                                               PARAM_INNER_SPREAD, PARAM_OUTER_SPREAD, 
                                               PARAM_INNER_SPREAD_RATIO, PARAM_OUTER_SPREAD_RATIO, PARAM_SPREAD_SPLIT,
                                               PARAM_CYCLES, PARAM_FILL }; 

  private double cyclesParam = 0;  // number of cycles (2*PI radians, circle circumference), if set to 0 then number of cycles is calculated automatically
  private double cusps = 5;
  private double cusp_size = 0.5;
  private double radius = 1.2;

  private double a_radius; // radius of circle "A" (stationary circle);
  private double b_radius; // radius of circle "B" (rolling circles, rolling around circumference of A)
  private double c_radius; // "radius" of point "C" (fixed distance from center of circle "B" to point "C")
  private double c_scale;  // c_radius/b_radius, = cusp_size + 1
  
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

  // vars for determining inner/outer via even-odd rule
  int sampleCount = 36000;
  int binCount = 720;
  ArrayList<ArrayList<Double>> theta_intersects = null;

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    // see class comments for derivation of a, b, c radius from radius, cusp, cusp_size parameters
    c_scale = cusp_size + 1;
    b_radius = radius/(cusps + c_scale + 1)  ;
    a_radius = b_radius * cusps;
    c_radius = b_radius * c_scale;
    
    if (cyclesParam == 0) {  
      cycles = cycles_to_close;
    }
    else {
      cycles = cyclesParam;
    }

    if (inner_mode == 10 || outer_mode == 10) {
      recalcCurveIntersects();
    }
    else {
      theta_intersects = null;
    }
  }

  public void recalcCurveIntersects() {
    // System.out.println("recalcing curves: " + this);
    theta_intersects = new ArrayList<ArrayList<Double>>(binCount);
    for (int i=0; i<binCount; i++) { 
      theta_intersects.add(new ArrayList<Double>());
    }
    ArrayList<Double> prev_tsects = null;
    for (int i=0; i<sampleCount; i++) {
      // double theta = map(i, 0,sampleCount, 0, 2 * M_PI);
      double theta = ((double)i/(double)sampleCount) * M_2PI;

      double x = ((a_radius + b_radius) * cos(theta)) - (c_radius * cos(((a_radius + b_radius)/b_radius) * theta));
      double y = ((a_radius + b_radius) * sin(theta)) - (c_radius * sin(((a_radius + b_radius)/b_radius) * theta));
      //      x = pAmount * x;
      //      y = pAmount * y;
      double r = sqrt(x*x + y*y);
      double angle = atan2(y, x);
      int anglebin =  (int)Math.floor(((angle + M_PI)/M_2PI) * (binCount-1));
      ArrayList<Double> tsects = theta_intersects.get(anglebin);

      // still rotating through same bin, merge results
      if (prev_tsects == tsects) {  
        // try ignoring for now -- should try averaging later?
      }
      else {
        tsects.add(r);  // autoboxing float r to Float object      
      }
      prev_tsects = tsects;
    }
  }

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double tin = atan2(pAffineTP.y, pAffineTP.x);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    double theta = cycles * tin;
    double raw_rin = sqrt((pAffineTP.x  * pAffineTP.x) + (pAffineTP.y * pAffineTP.y));
    double rin = spread_split * raw_rin;
        
    double x = ((a_radius + b_radius) * cos(theta)) - (c_radius * cos(((a_radius + b_radius)/b_radius) * theta));
    double y = ((a_radius + b_radius) * sin(theta)) - (c_radius * sin(((a_radius + b_radius)/b_radius) * theta));
    double r = sqrt(x*x + y*y);
    double t = atan2(y, x);

    if (fill != 0) { 
      r = r + (fill * (pContext.random() - 0.5));
    }

    double xin, yin;
    double rinx, riny;
    boolean outer_handled = false;
    boolean inner_handled = false;
    if (inner_mode == 10 || outer_mode == 10) {
      int anglebin =  (int)Math.floor(((tin + M_PI)/M_2PI) * (binCount-1));
      ArrayList<Double> tsects = theta_intersects.get(anglebin);
      int shorter = 0;
      int longer = 0;
      for (double ir : tsects) {
        if (ir <= raw_rin) { shorter++; }
        else { longer++; }
      }

      // modified Even-Odd rule: 
      //    cast ray from origin through incoming point to infinity
      //    count how many times curve intersects ray further out than incoming point (longer)
      //    if number is odd then point is inside, if number is even then point is outside
      if (longer % 2 == 0) {
        // point is outside
        // place on curve
        pVarTP.x += pAmount * x;
        pVarTP.y += pAmount * y;
      }
      else {
        // point is inside
        pVarTP.x += pAmount * pAffineTP.x;
        pVarTP.y += pAmount * pAffineTP.y;
      }
      
      /* more heavily modified Even-Odd rule, 
      if (longer % 2 == 0) {
        if (longer == 0) {  } // point is outside 
        else if (longer == 2 && tsects.size() == 2) {  } // point is inside (assumes curve is closed and encloses 0?) 
        else if (longer == 2 && tsects.size() == 3) {  } // point is outside (assumes curve is closed and encloses 0?) 
        else { } // point is weird? but usually outside
      }
      else {   } // point is inside
      */

    }
    else {

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
        case 9: 
          double rout = r + ((rin-r) * outer_spread);
          pVarTP.x += pAmount * rout * cos(t);
          pVarTP.y += pAmount * rout * sin(t);
          break;
        case 10:   // already handled before rin<->r comparison, and outer_mode conditional
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
        case 9: 
          // double rout = (rin * inner_spread) + (r * (1 - inner_spread));
          double rout = r + ((rin - r) * inner_spread);
          pVarTP.x += pAmount * rout * cos(t);
          pVarTP.y += pAmount * rout * sin(t);
          break;
        case 10:   // already handled before rin<->r comparison, and outer_mode conditional
          break;
        default:
          pVarTP.x += pAmount * x;
          pVarTP.y += pAmount * y;
          break;
      }
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
    return new Object[] { cusps, cusp_size, radius, 
                          unified_inner_outer, inner_mode, outer_mode, inner_spread, outer_spread,
                          inner_spread_ratio, outer_spread_ratio, spread_split,
                          cyclesParam, fill };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_CUSPS.equalsIgnoreCase(pName))
      cusps = pValue;
    else if (PARAM_CUSP_SIZE.equalsIgnoreCase(pName))
      cusp_size = pValue;
    else if (PARAM_RADIUS.equalsIgnoreCase(pName))
      radius = pValue;
    else if (PARAM_UNIFIED_INNER_OUTER.equalsIgnoreCase(pName)) {
      unified_inner_outer = (pValue == 0 ? 0 : 1);
    }
    else if (PARAM_OUTER_MODE.equalsIgnoreCase(pName)) {
      outer_mode = (int)floor(pValue);
      if (outer_mode > 15 || outer_mode < 0) { outer_mode = 0; }
    }
    else if (PARAM_INNER_MODE.equalsIgnoreCase(pName)) {
      inner_mode = (int)floor(pValue);
      if (inner_mode > 15 || inner_mode < 0) { inner_mode = 0; }
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
