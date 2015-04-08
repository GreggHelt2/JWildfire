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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.jwildfire.base.Tools;
import static org.jwildfire.base.mathlib.MathLib.sqrt;
import static org.jwildfire.base.mathlib.MathLib.fabs;
import static org.jwildfire.base.mathlib.MathLib.iabs;
import static org.jwildfire.base.mathlib.MathLib.max;
import static org.jwildfire.base.mathlib.MathLib.min;
import static org.jwildfire.base.mathlib.MathLib.sign;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

public class CompoundFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_VARIATION_A = "variationA";
  private static final String PARAM_VARIATION_B = "variationB";
  private static final String PARAM_SUM_A = "sumA";
  private static final String PARAM_SUM_B= "sumB";
  private static final String PARAM_MULT = "mult";
  private static final String PARAM_ARITHMETIC_MEAN = "arithmeticMean";
  private static final String PARAM_GEOMETRIC_MEAN = "geoMean";
  private static final String PARAM_HARMONIC_MEAN = "harmonicMean";
  private static final String PARAM_RADIAL = "radial";
  private static final String PARAM_SIGN_MODE = "signMode";

  private static final String[] paramNames = { PARAM_VARIATION_A,
                                               PARAM_VARIATION_B,
                                               PARAM_SUM_A,
                                               PARAM_SUM_B, 
                                               PARAM_MULT,
                                               PARAM_ARITHMETIC_MEAN,
                                               PARAM_GEOMETRIC_MEAN,
                                               PARAM_HARMONIC_MEAN, 
                                               PARAM_RADIAL,
                                               PARAM_SIGN_MODE };
                                               

  private VariationFunc varA = null;
  private VariationFunc varB = null;
  String defaultVarNameA = "sin";
  String defaultVarNameB = "asteria";
  int variationA = 0;  // position of variationA function in sorted Variation fucntion list
  int variationB = 0; // // position of variationB function in sorted Variation fucntion list
  double sumA = 1;   // separate linear weighting of A (otherwise could remove in favor of arithmetic mean)
  double sumB = 1;   // separate linear weighting of B (otherwise could remove in favor of arithmetic mean)
  double mult = 0;
  double amean = 0;
  double gmean = 0;
  double hmean = 0;
  double radial = 0;
  int signMode = 1;
  int signModeSign;
  int signModeAbs;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    XYZPoint pointA = new XYZPoint(pVarTP);
    XYZPoint pointB = new XYZPoint(pVarTP);
    XYZPoint affineA = new XYZPoint(pAffineTP);
    XYZPoint affineB = new XYZPoint(pAffineTP);
    varA.transform(pContext, pXForm, affineA, pointA, pAmount);
    varB.transform(pContext, pXForm, affineB, pointB, pAmount);
    double rx = sqrt(pointA.x * pointA.x + pointB.x * pointB.x);
    double ry = sqrt(pointA.y * pointA.y + pointB.y * pointB.y);
    int signX, signY;

    if (signModeAbs == 0)  { // random signX and signY (will end up with symmetry?)
      signX = pContext.random() >= 0.5 ? 1 : -1;
      signY = pContext.random() >= 0.5 ? 1 : -1;
    }
    else if (signModeAbs == 1) {  // 1 always 1 (negative inverts)
      signX = signModeSign;
      signY = signModeSign;
    }
    else if (signModeAbs == 2) {   // 2 ==> pointA decides (negative inverts)
      signX = signModeSign * sign(pointA.x);
      signY = signModeSign * sign(pointA.y);
    }
    else if (signModeAbs == 3)  { // 3 ==> pointB decides (negative inverts)
      signX = signModeSign * sign(pointB.x);
      signY = signModeSign * sign(pointB.y);
    }
    else if (signModeAbs == 4) {  // 4 ==> product of pointA and pointB decides (negative inverts)
      signX = signModeSign * sign(pointA.x * pointB.x);
      signY = signModeSign * sign(pointB.y * pointB.y);      
    }
    else if (signModeAbs == 5) {  // 5 ==> diff of pointA and pointB decides (negative inverts)
      signX = signModeSign * sign(pointA.x - pointB.x);
      signY = signModeSign * sign(pointA.y - pointB.y);
    }
    else if (signModeAbs == 6) {  // 6 ==> sum of pointA and pointB decides (negative inverts)
      signX = signModeSign * sign(pointA.x + pointB.x);
      signY = signModeSign * sign(pointA.y + pointB.y);
    }
    else if (signModeAbs == 7) {  // 7 ==> min of pointA and pointB decides (negative inverts)
      signX = signModeSign * sign(min(pointA.x, pointB.x));
      signY = signModeSign * sign(min(pointA.y, pointB.y));
    }
    else if (signModeAbs == 8) {  // 8 ==> max of pointA and pointB decides (negative inverts)
      signX = signModeSign * sign(max(pointA.x, pointB.x));
      signY = signModeSign * sign(max(pointA.y, pointB.y));
    }
    else  {  // same as signModeAbs = 0
      signX = pContext.random() >= 0.5 ? 1 : -1;
      signY = pContext.random() >= 0.5 ? 1 : -1;
    }
    
    pVarTP.x = (mult * signX * (fabs(pointA.x * pointB.x))) +
      (amean * signX * (fabs(pointA.x + pointA.x)/2)) + 
      (gmean * signX * sqrt(fabs(pointA.x * pointB.x))) +
      (hmean * signX * fabs((2 * pointA.x * pointB.x)/(pointA.x + pointB.x))) + 
      (radial * signX * rx) +
      (sumA * pointA.x) + (sumB * pointB.x);
    pVarTP.y = (mult * signY * (fabs(pointA.y * pointB.y))) +
      (amean * signY * (fabs(pointA.y + pointA.y)/2)) + 
      (gmean * signY * sqrt(fabs(pointA.y * pointB.y))) +
      (hmean * signY * fabs((2 * pointA.y * pointB.y)/(pointA.y + pointB.y))) + 
      (radial * signY * ry) + 
      (sumA * pointA.y) + (sumB * pointB.y);
    pVarTP.z = pAmount * pAffineTP.z;

  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { variationA, variationB, sumA, sumB, mult, amean, gmean, hmean, radial, signMode };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_VARIATION_A.equalsIgnoreCase(pName))
      variationA = Tools.FTOI(pValue);
    else if (PARAM_VARIATION_B.equalsIgnoreCase(pName))
      variationB = Tools.FTOI(pValue);
    else if (PARAM_SUM_A.equalsIgnoreCase(pName))
      sumA = pValue;
    else if (PARAM_SUM_B.equalsIgnoreCase(pName))
      sumB = pValue;
    else if (PARAM_MULT.equalsIgnoreCase(pName))
      mult = pValue;
    else if (PARAM_ARITHMETIC_MEAN.equalsIgnoreCase(pName))
      amean = pValue;
    else if (PARAM_GEOMETRIC_MEAN.equalsIgnoreCase(pName))
      gmean = pValue;
    else if (PARAM_HARMONIC_MEAN.equalsIgnoreCase(pName))
      hmean = pValue;
    else if (PARAM_RADIAL.equalsIgnoreCase(pName))
      radial = pValue;
    else if (PARAM_SIGN_MODE.equalsIgnoreCase(pName)) 
      signMode = Tools.FTOI(pValue);
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "compound";
  }


  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    List<String> sortedVarNames = new ArrayList<String>();
    sortedVarNames.addAll(VariationFuncList.getNameList());
    Collections.sort(sortedVarNames);
    String varNameA, varNameB;
    if (variationA <= 0)  { varNameA = defaultVarNameA; }
    else { varNameA = sortedVarNames.get(variationA); }
    if (variationB <= 0)  { varNameB = defaultVarNameB; }
    else { varNameB = sortedVarNames.get(variationB); }
       
    varA = VariationFuncList.getVariationFuncInstance(varNameA);
    varB = VariationFuncList.getVariationFuncInstance(varNameB);
    signModeSign = sign(signMode);
    signModeAbs = iabs(signMode);
    /*
    System.out.println("varA: " + varA.getName());
    System.out.println("varB: " + varB.getName());
    System.out.println("sign mode: " + signMode);
    System.out.println("signModeSign: " + signModeSign);
    System.out.println("signModeAbs: " + signModeAbs);
    */
    varA.init(pContext, pLayer, pXForm, pAmount);
    varB.init(pContext, pLayer, pXForm, pAmount);
    
  }

}
