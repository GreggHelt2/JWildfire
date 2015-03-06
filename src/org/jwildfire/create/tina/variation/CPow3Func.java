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

/*
     Ported from Apophysis Plugin to JWildfire Variation by CozyG
     original CPow2 Apophysis plugin written by xyrus02:
     http://sourceforge.net/p/apo-plugins/code/HEAD/tree/personal/georgkiehne/updated_for_x64/cpow2.c
 */
package org.jwildfire.create.tina.variation;

import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.exp;
import static org.jwildfire.base.mathlib.MathLib.floor;
import static org.jwildfire.base.mathlib.MathLib.log;
import static org.jwildfire.base.mathlib.MathLib.sin;

import org.jwildfire.base.Tools;
import org.jwildfire.create.tina.base.Layer;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class CPow3Func extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_R = "r";
  private static final String PARAM_CA = "ca";
  private static final String PARAM_DIVISOR = "divisor";
  private static final String PARAM_SPREAD = "spread";

  private static final String[] paramNames = { PARAM_R, PARAM_CA, PARAM_DIVISOR, PARAM_SPREAD };
  
  // Parameters
  private double r = 1.0;
  private double ca = 0.1;
  private int divisor = 1;
  private int spread = 1;

  // Internal fields set in init (based on parameters)
  double c;
  double d;
  double half_c;
  double half_d;
  double ang;
  double inv_spread;
  double full_spread;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
  /* CPow2 
    double sn, cs;
    double a = atan2(FTy, FTx);
    int n = rand() % VAR(cpow2_spread);
    if (a < 0) n++;
    a += 2*M_PI*n;
    if (cos(a*VAR(inv_spread)) < rand()*2.0/RAND_MAX - 1.0)
        a -= VAR(full_spread);
    double lnr2 = log(FTx*FTx + FTy*FTy); // logarithm * 2
    double r = VVAR * exp(VAR(half_c) * lnr2 - VAR(d) * a);
    fsincos(VAR(c) * a + VAR(half_d) * lnr2 + VAR(ang) * rand(),
        &sn, &cs);
    FPx += r * cs;
    FPy += r * sn;
    */

    double a = pAffineTP.getPrecalcAtanYX();
    int n = (int)(pContext.random() * spread);
    if (a < 0) { n++; }
    a += 2 * M_PI * n;
    if (cos(a * inv_spread) < (pContext.random() * 2.0 - 1.0)) {
      a -= full_spread;
    }
    double lnr2 = log(pAffineTP.getPrecalcSumsq());
    double r = pAmount * exp(half_c * lnr2 - d * a);
    
    double ang2 = c * a * half_d * lnr2 * ang * pContext.random();
    pVarTP.x += r * cos(ang2);
    pVarTP.y += r * sin(ang2);

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
    return new Object[] { r, ca, divisor, spread };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    // System.out.println("called CPow2.setParameter(), " + pName + ", " + pValue);
    if (PARAM_R.equalsIgnoreCase(pName))
      r = pValue;
    else if (PARAM_CA.equalsIgnoreCase(pName))
      ca = pValue;
    else if (PARAM_DIVISOR.equalsIgnoreCase(pName))
      divisor = Tools.FTOI(pValue);
    else if (PARAM_SPREAD.equalsIgnoreCase(pName)) {
      spread = Tools.FTOI(pValue);
    }
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "cpow3";
  }

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    /*
    VAR(ang) = 2*M_PI / ((double) VAR(cpow2_divisor));
    VAR(c) = VAR(cpow2_r) * cos(M_PI/2*VAR(cpow2_a)) / ((double) VAR(cpow2_divisor));
    VAR(d) = VAR(cpow2_r) * sin(M_PI/2*VAR(cpow2_a)) / ((double) VAR(cpow2_divisor));
    VAR(half_c) = VAR(c) / 2;
    VAR(half_d) = VAR(d) / 2;
    VAR(inv_spread) = 0.5 / VAR(cpow2_spread);
    VAR(full_spread) = 2*M_PI*VAR(cpow2_spread);
    */
    c = r * cos(M_PI / 2 * ca) / (double)divisor;
    d = r * sin(M_PI / 2 * ca) / (double)divisor;
    half_c = c / 2;
    half_d = d / 2;
    ang = 2.0 * M_PI / (double)divisor;
    inv_spread = 0.5 / spread;
    full_spread = 2 * M_PI * spread;
  }

}
