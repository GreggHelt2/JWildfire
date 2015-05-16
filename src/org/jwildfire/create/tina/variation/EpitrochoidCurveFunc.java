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
public class EpitrochoidCurveFunc extends AbstractPolarCurveFunc {
  private static final long serialVersionUID = 1L;
  //  public enum IntersectMode { SIMPLE_RADIUS, MIN_RADIUS, MAX_RADIUS, EVEN_ODD, MODIFIED_EVEN_ODD }
  
  // cusps, radius, and cusp_size are used to determine:
  //    a_radius: radius of stationary circle
  //    b_radius: radius of rolling circle
  //    c_radius: radius of rolling point attached to center of rolling circle
  private static final String PARAM_RADIUS = "radius";  // = a_radius + b_radius + c_radius
  private static final String PARAM_CUSPS = "cusps";  //  = a_radius/b_radius 
  private static final String PARAM_CUSP_SIZE = "cusp_size";  //  = (c_radius/b_radius) - 1
  private static final String PARAM_CUSP_DIVISOR = "cusp_divisor"; 
  private static final String[] additionalParamNames = { PARAM_CUSPS, PARAM_CUSP_SIZE, PARAM_CUSP_DIVISOR, PARAM_RADIUS }; 

  private double cusps = 5;
  private double cusp_size = 0.5;
  private double radius = 1.2;
  private double cusp_divisor = 1.0;

  private double a_radius; // radius of circle "A" (stationary circle);
  private double b_radius; // radius of circle "B" (rolling circles, rolling around circumference of A)
  private double c_radius; // "radius" of point "C" (fixed distance from center of circle "B" to point "C")
  private double c_scale;  // c_radius/b_radius, = cusp_size + 1
  private double k;  // cusps/cusp_divisor
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {

    // see class comments for derivation of a, b, c radius from radius, cusp, cusp_size, cusp_divisor, k parameters
    k = cusps/cusp_divisor;
    c_scale = cusp_size + 1;
    b_radius = radius/(k + c_scale + 1)  ;
    a_radius = b_radius * k;
    c_radius = b_radius * c_scale;
    if (DEBUG) {
      System.out.println("a: " + a_radius + ", b: " + b_radius + ", c: " + c_radius + 
              ", a/b: " + (a_radius/b_radius) + ", c/b: " + (c_radius/b_radius) + ", cycles: " + cycles) ;
    }
    super.init(pContext, pLayer, pXForm, pAmount);  // calls recalcCycles() and recalcCurveIntersects()
  }

  @Override
  public void recalcCycles() {
    if (cycles_param == 0) {  
      if (k % 1 == 0) { // k is an integer
        cycles_to_close = 1;
        //  cycles = cycles_to_close;
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
        // cycles = cycles_to_close * metacycles;
      }
      else {
        // cannot figure out if k is rational (could be irrational or just not easily determined)
        // should really manually set cycles in this case
        cycles_to_close = 0;
      }
      if (cycles_to_close > 0) {
        cycles = cycles_to_close;
      }
      else {
        // cycles to close is unknown, so just pick a high number relative to cusps and cusps_divisor
        cycles = cusps * cusp_divisor;
        if (cycles < cusps || cycles < cusp_divisor) { cycles = Math.max(cusps, cusp_divisor); }
        if (cycles < 3) { cycles = 3; }
      }
    }
    else {
      cycles = cycles_param;
    }
    // radians = cycles * 2 * M_PI;
    super.recalcCycles();
  }

  /* called from super.transform() */
  public void calcCurvePoint(FlameTransformationContext pContext, double theta, XYZPoint pResult) {
    // pResult should be zero'd out before getting here
    pResult.clear();
    pResult.x = ((a_radius + b_radius) * cos(theta)) - (c_radius * cos(((a_radius + b_radius)/b_radius) * theta));
    pResult.y = ((a_radius + b_radius) * sin(theta)) - (c_radius * sin(((a_radius + b_radius)/b_radius) * theta));
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
    return joinArrays(new Object[] { cusps, cusp_size, cusp_divisor, radius }, super.getParameterValues());
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_CUSPS.equalsIgnoreCase(pName))
      cusps = pValue;
    else if (PARAM_CUSP_SIZE.equalsIgnoreCase(pName))
      cusp_size = pValue;
    else if (PARAM_CUSP_DIVISOR.equalsIgnoreCase(pName))
      cusp_divisor = pValue;
    else if (PARAM_RADIUS.equalsIgnoreCase(pName))
      radius = pValue;
    else {
      super.setParameter(pName, pValue);
    }
  }

  @Override
  public String getName() {
    return "epitrochoid_curve";
  }

}
