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

  private static final String PARAM_RHO = "rho";
  private static final String PARAM_ALPHA = "alpha";

  private static final String[] paramNames = { PARAM_RHO, PARAM_ALPHA };

  private double rho = 1.0;
  private double alpha = 1.0;

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double theta = atan2(pAffineTP.y, pAffineTP.x);
    double r = pAffineTP.getPrecalcSqrt();
    if (r <= rho) {
      pVarTP.x += pAmount * (rho * r) * cos(alpha * theta);
      pVarTP.y += pAmount * (rho * r) * sin(alpha * theta);
      pVarTP.z += theta;
    }
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { rho, alpha };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_RHO.equalsIgnoreCase(pName))
      rho = pValue;
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
