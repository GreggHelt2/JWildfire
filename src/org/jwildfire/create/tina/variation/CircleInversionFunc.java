package org.jwildfire.create.tina.variation;

import static java.lang.Math.abs;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sqr;
import static org.jwildfire.base.mathlib.MathLib.pow;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

/**
 * CircleInversionFunc, a variation for inversion in a circle
 * Initially similar to spherical, but with a specific radius parameter
 *     (whereas spherical effectively sets radius = sqrt(pAmount))
 *  and also includes setting xorigin and yorigin
 *     (whereas spherical relies on xtransforms affine transform, but 
 *      then if want circle origin != (0, 0) then need to keep 
 *      pre- and post- coefficients as mirrors of each other 
 *      (in other words "move" circle to origin, do inversion, then reverse move) )
 *  also includes "draw_circles" param
 *      if 0 > draw_circles < 1, then that fraction of incoming points is used to 
 *      draw circle of inversion rather than doing the actual inversion >
 *  In addition to standard circle inversion, can also be used for p-circle inversion as 
 *      described by Ramirez et al in "Generating Fractal Patterns by Using P-Circle Inversion" (2015)
 */
public class CircleInversionFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  public static final String PARAM_RADIUS= "radius";
  public static final String PARAM_XORIGIN = "xorigin";
  public static final String PARAM_YORIGIN = "yorigin";
  public static final String PARAM_P= "p";
  public static final String PARAM_DRAW_CIRCLE = "draw_circle";
  
  double r = 1;
  double a = 0;
  double b = 0;
  double p = 2;
  double draw_circle = 0;

  private static final String[] paramNames = { PARAM_RADIUS, PARAM_XORIGIN, PARAM_YORIGIN, PARAM_P, PARAM_DRAW_CIRCLE };

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double xin = pAffineTP.x;
    double yin = pAffineTP.y;
    double xdiff = xin-a;
    double ydiff = yin-b;
    double iscale;
    if (p == 2) {
      iscale = (r*r)/(sqr(xdiff) + sqr(ydiff));
    }
    else {
      iscale = (r*r)/ pow( (pow(abs(xdiff),p) + pow(abs(ydiff),p)), 2.0/p);
    }
    double xout = a + (iscale * xdiff);
    double yout = b + (iscale * ydiff);
    pVarTP.x += pAmount * xout;
    pVarTP.y += pAmount * yout;
    if (draw_circle > 0) {
      double rnd = pContext.random();
      if (rnd < draw_circle) {
        double theta = pContext.random() * M_2PI;
        pVarTP.x = a + r * cos(theta);
        pVarTP.y = b + r * sin(theta);
      }
    }
  }
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
  }

  @Override
  public String getName() {
    return "circle_inversion";
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { r, a, b, p, draw_circle };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_RADIUS.equalsIgnoreCase(pName)) {
      r = pValue;
    }
    else if (PARAM_XORIGIN.equalsIgnoreCase(pName)) {
      a = pValue;
    }
    else if (PARAM_YORIGIN.equalsIgnoreCase(pName)) {
      b = pValue;
    }
    else if (PARAM_P.equalsIgnoreCase(pName)) {
      p = pValue;
    }
    else if (PARAM_DRAW_CIRCLE.equalsIgnoreCase(pName)) {
      draw_circle = pValue;
    }
    else
      throw new IllegalArgumentException(pName);
  }

}
