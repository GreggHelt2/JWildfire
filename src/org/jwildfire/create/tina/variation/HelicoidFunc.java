package org.jwildfire.create.tina.variation;

import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.atan2;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

/*
Helicoid: 
    x = rho * cos (alpha * theta) 
    y = rho * sin (alpha * theta) 
    z = theta, \ 
*/
public class HelicoidFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_RHO_MAX = "rhoMax";
  private static final String PARAM_RHO_MIN = "rhoMin";
  private static final String PARAM_ALPHA = "alpha";

  private static final String[] paramNames = { PARAM_RHO_MAX, PARAM_RHO_MIN, PARAM_ALPHA };

  private double rhoMax = 1.0;
  private double rhoMin = 0.3;
  private double alpha = 1.0;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double theta = atan2(pAffineTP.y, pAffineTP.x);
    double r = pAffineTP.getPrecalcSqrt();
    if (r <= rhoMax && r >= rhoMin) {
      pVarTP.x += pAmount * (rhoMax * r) * cos(alpha * theta);
      pVarTP.y += pAmount * (rhoMax * r) * sin(alpha * theta);
      pVarTP.z += theta;
    }
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { rhoMax, rhoMin, alpha };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_RHO_MAX.equalsIgnoreCase(pName))
      rhoMax = pValue;
    else if (PARAM_RHO_MIN.equalsIgnoreCase(pName))
      rhoMin = pValue;
    else if (PARAM_ALPHA.equalsIgnoreCase(pName))
      alpha = pValue;
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "helicoid";
  }

}
