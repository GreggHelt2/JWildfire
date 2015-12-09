package org.jwildfire.create.tina.variation;

import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.atan2;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.fabs;
import static org.jwildfire.base.mathlib.MathLib.pow;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.sqrt;
import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

/** 
 * based on Apophysis "juni" plugin by Xyrus
 *   converted to JWildfire plugin by CozyG
 *   see JWildfire forum topic for more details: 
 *   http://jwildfire.org/forum/viewtopic.php?f=23&t=1975
 */
public class JuniFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_POWER = "power";
  private static final String PARAM_STRETCH = "stretch";
  private static final String PARAM_INVERT = "invert";

  private static final String[] paramNames = { PARAM_POWER, PARAM_STRETCH, PARAM_INVERT };

  private int power = 2;
  private double stretch = 1.0;
  private boolean invert = false;
  double absn;
  double cn;

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    absn = (int)fabs(power);
    cn = (stretch / (power - 1)) / 2.0;
  }

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double xin = pAffineTP.x;
    double yin = pAffineTP.y;
    double zin = pAffineTP.z;
    double zz = zin / absn;
    double r2 = (xin * xin) + (yin * yin);
    double r = pAmount * pow(r2 + (zz * zz), cn);
    double pasq = pAmount * pAmount;

    if ( ((r2 > pasq || r2 == 0) && !invert) ||
         ((r2 < pasq && r2 != 0) && invert) ) {
      pVarTP.x += pAmount * xin;
      pVarTP.y += pAmount * yin;
    }
    else {
      double theta = (int)(pContext.random() * absn);
      double alpha = (atan2(yin, xin) + (M_2PI * theta)) / power;
      double gamma = r * sqrt(r2);
      pVarTP.x += gamma * cos(alpha);
      pVarTP.y += gamma * sin(alpha);
    }
    pVarTP.z += r * zz;
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { power, stretch, (invert ? 1 : 0) };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_POWER.equalsIgnoreCase(pName))  {
      power = (int)pValue;
      // avoiding two possible divide-by-zero errors
      if (power == 0 || power == 1) { power = -1; }
    }
    else if (PARAM_STRETCH.equalsIgnoreCase(pName))  {
      stretch = pValue;
    }
    else if (PARAM_INVERT.equalsIgnoreCase(pName)) {
      invert = (pValue >= 1);
    }
    else {
      throw new IllegalArgumentException(pName);
    }
  }

  @Override
  public String getName() {
    return "juni_GH_v1";
  }
}
