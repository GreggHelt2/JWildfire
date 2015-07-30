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

import java.math.BigInteger;

import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sin;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class RhodoneaCurveFunc extends AbstractPolarCurveFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_KNUMER = "knumer";
  private static final String PARAM_KDENOM = "kdenom";
  private static final String PARAM_RADIAL_OFFSET = "radial_offset";

  private static final String[] additionalParamNames = { PARAM_KNUMER, PARAM_KDENOM, PARAM_RADIAL_OFFSET };

  private double kn = 3; // numerator of k in rose curve equations,   k = kn/kd
  private double kd = 4; // denominator of k in rose curve equations, k = kn/kd
  private double k; // k = kn/kd
  private double radial_offset = 0; // often called "c" in rose curve modifier equations
  private double petal_count = 0; // petal_count = 0 means petal count is unknown
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    k = kn / kd;
    super.init(pContext, pLayer, pXForm, pAmount);  // calls recalcCycles() and recalcCurveIntersects()
  }
  
  @Override
  public void recalcCycles() {
    // attempt to calculate minimum cycles manually, or reasonable upper bound if unsure
    if ((k % 1) == 0) { // k is integer 
      if ((k % 2) == 0) { // k is even integer, will have 2k petals
        cycles_to_close = 1; // (2PI)
        petal_count = 2 * k;
      }
      else { // k is odd integer, will have k petals (or sometimes 2k with offset) and cycles closes in 1PI
        if (radial_offset != 0) { // unless adding radial_offset, then will need 2PI cycles to close
          cycles_to_close = 1;
        } 
        else {
          cycles_to_close = 0.5;
        } // (1PI)
        petal_count = k;
      }
    }
    else if ((kn % 1 == 0) && (kd % 1 == 0)) {
      // if kn and kd are integers,
      //   determine if kn and kd are relatively prime (their greatest common denominator is 1)
      //   using builtin gcd() function for BigIntegers in Java
      // and if they're not, make them
      BigInteger bigkn = BigInteger.valueOf((long) kn);
      BigInteger bigkd = BigInteger.valueOf((int) kd);
      int gcd = bigkn.gcd(bigkd).intValue();
      if (gcd != 1) {
        kn = kn / gcd;
        kd = kd / gcd;
      }

      // paraphrased from http://www.encyclopediaofmath.org/index.php/Roses_%28curves%29:
      //    If kn and kd are relatively prime, then the rose consists of 2*kn petals if either kn or kd are even, and kn petals if both kn and kd are odd
      //
      // paraphrased from http://mathworld.wolfram.com/Rose.html:
      //    If k=kn/kd is a rational number, then the curve closes at a polar angle of theta = PI * kd if (kn * kd) is odd, and 2 * PI * kd if (kn * kd) is even
      if ((kn % 2 == 0) || (kd % 2 == 0)) {
        petal_count = 2 * kn;
        cycles_to_close = kd; // 2 * PI * kd
      }
      else {
        petal_count = kn;
        cycles_to_close = kd / 2; // PI * kd
      }
    }

    // additional special cases:
    /*  
    // superceded by relative prime conditional?
    else if (((k * 2) % 1) == 0) { // k is a half-integer (1/2, 3/2, 5/2, etc.), will have 4k petals
    cycles_to_close = 2;  // (4PI)
    petal_count = 4*k;
    }
    */
    /* 
    // superceded by relative prime conditional?
    // If k can be expressed as kn/3, where n is an integer not divisible by 3, the curve will be rose-shaped with n petals if n is odd and 2n petals if n is even.
    //    (case where kn is integer divisible by three is already handled above, where k is integer)
    else if (((k * 3) % 1) == 0) {
    double basekn = k * 3;
    if ((basekn % 2) == 0) {
    cycles_to_close = 3;
    petal_count = 2 * basekn;
    }
    else  {
    cycles_to_close = 1.5;
    petal_count = basekn;
    }
    }
    */
    /*
    // superceded by relative prime conditional????
    If k can be expressed as n±1/6, where n is a nonzero integer, the curve will be rose-shaped with 12k petals.
    else if {
    }
    */
    else {
      //     if one or both of kn and kd are non-integers, then the above may still be true (k may still be [effectively] rational) but haven't 
      //          figured out a way to determine this.
      //    could restrict kn and kd to integers to simplify, but that would exclude a huge space of interesting patterns
      //    could set cycles extremely high, but if k is truly irrational this will just approarch a uniform distribution across a circle, 
      //                and also exclude a large space of interesting patterns with non-closed curves
      //    so for now keep kn and kd as continuous doubles, and just pick a large but not huge number for cycles
      // 
      //    realistically in this case it is better for user to fiddle with manual cycles setting to get a pattern they like
      //
      cycles_to_close = 0;

    }
    if (cycles_param == 0) {
      if (cycles_to_close > 0) {
        // use auto-calculation of cycles (2*PI*radians) to close the curve, 
        //     and metacycles for how many closed curves to cycle through
        // cycles = cycles_to_close * metacycles;
        cycles = cycles_to_close;  // moved metacycle calcs to superclass
      }
      else {
        // cycles_to_close is unknown, or curve never closes
        // setting maximum to (2*kn*kd) arbitrary number relatively higher than kn and kd
        //    and set lower threshold of 16, just want to keep number higher if kn and kd are small but potentially irrational
        cycles = Math.max(2 * kn * kd, 16);
      }
    }
    else {
      // manually set number of cycles (cycles are specified in 2*PI*radians)
      cycles = cycles_param;
    }
    super.recalcCycles();
  }


  /* called from super.transform() */
  public void calcCurvePoint(FlameTransformationContext pContext, double theta, XYZPoint pResult) {
    // pResult should be zero'd out before getting here
    // not factoring in cycle_offset yet, should incorporate into superclass params?)
    //    theta = cycles * (t + (cycle_offset * M_2PI))
    double r = cos(k * theta) + radial_offset; // radius of rose curve
    pResult.x = r * cos(theta);
    pResult.y = r * sin(theta);
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
    return joinArrays(new Object[] { kn, kd, radial_offset }, super.getParameterValues());
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_KNUMER.equalsIgnoreCase(pName))
      kn = pValue;
    else if (PARAM_KDENOM.equalsIgnoreCase(pName))
      kd = pValue;
    else if (PARAM_RADIAL_OFFSET.equalsIgnoreCase(pName))
      radial_offset = pValue;
    else {
      super.setParameter(pName, pValue);
    }
  }

  @Override
  public String getName() {
    return "rhodonea_curve";
  }

}
