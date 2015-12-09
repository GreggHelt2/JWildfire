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

/**
  Epitrochoid
  Implemented by CozyG, started April 2015
  For reference see http://en.wikipedia.org/wiki/Epitrochoid 
  The basic equations for an epitrochoid curve are identical to the basic equations in the Spirograph variation 
      (the Spirograph toy is a system of gears to trace epitrochoid and hypotrochoid curves on paper)
  However, how the Epitrochoid and Spirograph variations transform incoming points based on the curve can be very different

 This Epitrochoid variation transforms incoming points differently, 
     depending on whether they fall inside or outside the epitrochoid curve
 In this version for both inside and outside points there are 9 different rendering modes 
     to choose from (inner_mode and outer_mode values of 0-8)

 ---------------------------------------------------------

 Relation of radius, cusps, and cusp_size params to epitrochoid equations
 Translating from equation parameters to input params that better map to obvious features
  if a = radius of stationary circle A
     b = radius of outer rolling circle B
     c = distance from point C attached to circle B to center of B

  then max radius of the epitrochoid will be a+b+c
       min radius of the epitrochoid will be a+b-c
    
  cusps = k = a/b is the number of retrograde cycles in the closed curve
  cscale = s = c/b is the ratio of (distance of point C from origin of B) to (radius of B)
  cusp_size = cscale-1 is cscale adjusted, so that at cusp_size == 0 the cusps terminate in a point instead of having a retrograde cycle

  so make ajustable values be:
       radius = r = a+b+c
       cusps  = k = a/b
       cusp_size = s-1 = (c/b)-1
 ---------------------------------------------------------

 PS: I'm working on a newer, enhanced version of Epitrochoid re-implemented on top of a more general polar curve class 
      for sharing much of the code for a set of different curves based on polar coordinates. But I'm not sure when that 
      work will be ready for a JWildfire release, so I thought I'd put out this earlier version as a custom variation for now
*/
public class EpitrochoidSimpleFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  // cusps, radius, and cusp_size are used to determine:
  //    a_radius: radius of stationary circle
  //    b_radius: radius of rolling circle
  //    c_radius: radius of rolling point attached to center of rolling circle
  private static final String PARAM_RADIUS = "radius";  // = a_radius + b_radius + c_radius
  private static final String PARAM_CUSPS = "cusps";  //  = a_radius/b_radius 
  private static final String PARAM_CUSP_SIZE = "cusp_size";  //  = (c_radius/b_radius) - 1
  private static final String PARAM_CUSP_DIVISOR = "cusp_divisor"; 
  private static final String PARAM_INNER_MODE = "inner_mode";
  private static final String PARAM_OUTER_MODE = "outer_mode";
  private static final String PARAM_INNER_SPREAD = "inner_spread";
  private static final String PARAM_OUTER_SPREAD = "outer_spread";
  private static final String PARAM_INNER_SPREAD_RATIO = "inner_spread_ratio";
  private static final String PARAM_OUTER_SPREAD_RATIO = "outer_spread_ratio";
  private static final String PARAM_SPREAD_SPLIT = "spread_split";
  private static final String PARAM_THICKNESS = "thickness";
  private static final String PARAM_DIFF_MODE = "diff_mode";

  private static final String[] paramNames = { PARAM_RADIUS, PARAM_CUSPS, PARAM_CUSP_SIZE, PARAM_CUSP_DIVISOR, 
                                               PARAM_INNER_MODE, PARAM_OUTER_MODE, 
                                               PARAM_INNER_SPREAD, PARAM_OUTER_SPREAD, 
                                               PARAM_INNER_SPREAD_RATIO, PARAM_OUTER_SPREAD_RATIO, PARAM_SPREAD_SPLIT, PARAM_THICKNESS, PARAM_DIFF_MODE };
  private static int MIN_MODE = 0;
  private static int MAX_MODE = 8;
  private static int DEFAULT_MODE = MIN_MODE;

  private double cyclesParam = 0;  // number of cycles (2*PI radians, circle circumference), if set to 0 then number of cycles is calculated automatically
  private double cusps = 6;
  private double cusp_size = 1;
  private double radius = 1;
  private double cusp_divisor = 1;

  private double a_radius; // radius of circle "A" (stationary circle);
  private double b_radius; // radius of circle "B" (rolling circle, rolling around outside of circumference of A)
  private double c_radius; // "radius" of point "C" (fixed distance from center of circle "B" to point "C")
  private double c_scale;  // c_radius/b_radius, = cusp_size + 1
  private double k;  // cusps/cusp_divisor
  
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

    // see class comments for derivation of a, b, c radius from radius, cusp, cusp_size, cusp_divisor, k parameters
    k = cusps/cusp_divisor;
    c_scale = cusp_size + 1;
    // b_radius = radius/(k + c_scale + 1);
    b_radius = -1 * radius / (1 - k - c_scale);
    a_radius = b_radius * k;
    c_radius = b_radius * c_scale;
    
    if (cyclesParam == 0) {  
      if (k % 1 == 0) { // k is an integer
        cycles_to_close = 1;
        cycles = cycles_to_close;
      }
      else if ((cusps % 1 == 0) && (cusp_divisor % 1 == 0)) {
        // if cusp and cusp_divisor are integers, then k is rational
        // reduce cusps and cups_divisor to simplest terms using greatest common denominator
        //   using builtin gcd() function for BigIntegers in Java
        BigInteger cusp_big = BigInteger.valueOf((long)cusps);
        BigInteger cusp_divisor_big = BigInteger.valueOf((long)cusp_divisor);
        int gcd = cusp_big.gcd(cusp_divisor_big).intValue();
        // simplify terms using gcd
        //     (and if no common denominator then gcd = 1, so values wills be unchanged)
        double cusps_normalized = cusps/gcd;
        double cusp_divisor_normalized = cusp_divisor/gcd;
        cycles_to_close = cusp_divisor_normalized;
        cycles = cycles_to_close;
      }
      else {
        // cannot figure out if k is rational (could be irrational or just not easily determined)
        // pick a high number relative to cusps and cusps_divisor
        cycles_to_close = 0;
        cycles = cusps * cusp_divisor;
      }
    }
    else {
      cycles = cyclesParam;
    }
    radians = cycles * 2 * M_PI;

    if (DEBUG) {
      System.out.println("a: " + a_radius + ", b: " + b_radius + ", c: " + c_radius + 
              ", a/b: " + (a_radius/b_radius) + ", c/b: " + (c_radius/b_radius) + ", cycles: " + cycles) ;
    }

  }
  
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double xin = pAffineTP.x;
    double yin = pAffineTP.y;
    double tin = atan2(yin, xin);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    double theta = cycles * tin;
    double raw_rin = sqrt((xin * xin) + (yin * yin));
    double rin = spread_split * raw_rin;

    double x = ((a_radius + b_radius) * cos(theta)) - (c_radius * cos(((a_radius + b_radius)/b_radius) * theta));
    double y = ((a_radius + b_radius) * sin(theta)) - (c_radius * sin(((a_radius + b_radius)/b_radius) * theta));
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
    return new Object[] { radius, cusps, cusp_size, cusp_divisor, 
                          inner_mode, outer_mode, inner_spread, outer_spread,
                          inner_spread_ratio, outer_spread_ratio, spread_split, thickness, (diff_mode ? 1 : 0)  };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_RADIUS.equalsIgnoreCase(pName))
      radius = pValue;
    else if (PARAM_CUSPS.equalsIgnoreCase(pName))
      cusps = pValue;
    else if (PARAM_CUSP_SIZE.equalsIgnoreCase(pName))
      cusp_size = pValue;
    else if (PARAM_CUSP_DIVISOR.equalsIgnoreCase(pName))
      cusp_divisor = pValue;
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
    return "epitrochoid_simple";
  }

}
