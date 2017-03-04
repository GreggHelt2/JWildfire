package org.jwildfire.create.tina.variation;

import static java.lang.Math.abs;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sqr;
import static org.jwildfire.base.mathlib.MathLib.pow;
import static org.jwildfire.base.mathlib.MathLib.sqrt;

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
public class CircleInversionFunc extends VariationFunc implements Guides {
  private static final long serialVersionUID = 1L;

  public static final String PARAM_RADIUS= "radius";
  public static final String PARAM_XORIGIN = "xorigin";
  public static final String PARAM_YORIGIN = "yorigin";
  public static final String PARAM_ZORIGIN = "zorigin";
  
  public static final String PARAM_ZMODE = "zmode";
  public static final String PARAM_INVERSION_MODE = "imode";
  public static final String PARAM_HIDE_UNINVERTED = "hide_uninverted";
  public static final String PARAM_RING_MIN = "ring_min";
  public static final String PARAM_RING_MAX = "ring_max";
  public static final String PARAM_P= "p";
  public static final String PARAM_DRAW_CIRCLE = "draw_circle";
  public static final String PARAM_GUIDES_ENABLED = "guides_enabled";
  
  public static int STANDARD = 0;
  public static int EXTERNAL_INVERSION_ONLY = 1;
  public static int INTERNAL_INVERSION_ONLY = 2;
  // outer inversion is like external inversion, except 
  //    leaves unchanged any  "external" regions that are distance < radius (inner hole)
  public static int EXTERNAL_RING_INVERSION_ONLY = 3;
  public static int INTERNAL_RING_INVERSION_ONLY = 4;

  
  public boolean DESIGN_GUIDE_MODE = false;
  double r = 1;
  double a = 0;
  double b = 0;
  double c = 0;
  double p = 2;
  double ring_min_ratio = 0;
  double ring_max_ratio = 1;
  double ring_min;
  double ring_max;
  double draw_circle = 0;
  boolean guides_enabled = true;
  int inversion_mode = STANDARD;
  boolean hide_uninverted = false;
  boolean zmode = false;

  private static final String[] paramNames = { 
    PARAM_RADIUS, PARAM_XORIGIN, PARAM_YORIGIN, PARAM_ZORIGIN, PARAM_ZMODE, 
    PARAM_INVERSION_MODE, PARAM_HIDE_UNINVERTED, PARAM_RING_MIN, PARAM_RING_MAX, PARAM_P, PARAM_DRAW_CIRCLE, PARAM_GUIDES_ENABLED
  };
  
  public void setDrawGuides(boolean draw) {
    DESIGN_GUIDE_MODE = draw;
  }
  
  public boolean getDrawGuides() { return DESIGN_GUIDE_MODE; }

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double xin = pAffineTP.x;
    double yin = pAffineTP.y;
    double zin = pAffineTP.z;
    double xdiff = xin-a;
    double ydiff = yin-b;
    double zdiff = zin-c;
    double iscale;
    double hide_ratio = 0.99;
    
    if (DESIGN_GUIDE_MODE && guides_enabled) {
//      double rnd = pContext.random() * 4.0;
      double rnd = pContext.random();
      // if (rnd < 2) { // 0 <= rnd < 2
        double theta = rnd * M_2PI;
        // pVarTP.x += xin + a + r * cos(theta);
        // pVarTP.y += yin + b + r * sin(theta);
        // pVarTP.x += xin + r * cos(theta);
        // pVarTP.y += yin + r * sin(theta);
        if (zmode) {
          double split = pContext.random() * 3.1;
          if (split < 1) {
            pVarTP.x += a + r * cos(theta);
            pVarTP.y += b + r * sin(theta);
            pVarTP.z = c;
          }
          else if (split < 2) {
            pVarTP.x += a + r * cos(theta);
            pVarTP.y = b;
            pVarTP.z += c + r * sin(theta);
          }
          else if (split < 3) {
            pVarTP.x = a;
            pVarTP.y += b + r * sin(theta);
            pVarTP.z += c + r * cos(theta);
          }
          else {
            // double offset_rnd = pContext.random() * 0.05;
            // offset_rnd -= 0.05/2;
            //pVarTP.x += a + offset_rnd;
            //pVarTP.y += b + offset_rnd;
            // pVarTP.z += c + offset_rnd;
            pVarTP.x += a;
            pVarTP.y += b;
            pVarTP.z += c;
          }
        }
        else {
          double split = pContext.random() * 1.1;
          if (split < 1) {
            pVarTP.x += a + r * cos(theta);
            pVarTP.y += b + r * sin(theta);
            pVarTP.z += c;
          }
          else {
            pVarTP.x += a;
            pVarTP.y += b;
            pVarTP.z += c;
          }
        }
     // }
  /*    else if (rnd < 3) {  // 2 <= rnd < 3
        // draw rmax circle
        double theta = rnd * M_2PI;
        pVarTP.x += a + ring_max * cos(theta);
        pVarTP.y += b + ring_max * sin(theta);
      }
      else {  // 3 <= rnd < 4
        // draw rmin circle (or point)
        double theta = rnd * M_2PI;
        pVarTP.x += a + ring_min * cos(theta);
        pVarTP.y += b + ring_min * sin(theta);
      }
          */
      return;
    }
    
    if (draw_circle > 0) {
      double rnd = pContext.random();
      if (rnd < draw_circle) {
        double theta = pContext.random() * M_2PI;
        // pVarTP.x = a + r * cos(theta);
        // pVarTP.y = b + r * sin(theta);
        pVarTP.x += a + r * cos(theta);
        pVarTP.y += b + r * sin(theta);
        return;
      }
    }
    boolean do_inversion;
    if (inversion_mode == EXTERNAL_INVERSION_ONLY) {
      // only do inversion if input point is outside of circle
      // calc distance from circle origin
      double d = sqrt(sqr(xdiff) + sqr(ydiff));
      // do_inversion = (d > r);
      do_inversion = (d > ring_max) || (d < ring_min);
      
    }
    else if (inversion_mode == INTERNAL_INVERSION_ONLY) {
      // only do inversion if input point is inside of circle
      double d = sqrt(sqr(xdiff) + sqr(ydiff));
      // do_inversion = (d < r);
      do_inversion = (d < ring_max) && ((d > ring_min) || (ring_min == 0));
    }
    else { // default to STANDARD mode ==> always do inversion
      do_inversion = true;
    }
    if (do_inversion) {
      if (p == 2) {
        // iscale = (r*r)/(sqr(xdiff) + sqr(ydiff));
        iscale = (r*r)/(sqr(xdiff) + sqr(ydiff) + sqr(zdiff));
      }
      else {
        iscale = (r*r)/ pow( (pow(abs(xdiff),p) + pow(abs(ydiff),p)), 2.0/p);
      }
      double xout = a + (iscale * xdiff);
      double yout = b + (iscale * ydiff);
      double zout = c + (iscale * zdiff);
      pVarTP.x += pAmount * xout;
      pVarTP.y += pAmount * yout;
      pVarTP.z += pAmount * zout;
      // if (pContext.isPreserveZCoordinate()) {
      //  pVarTP.z += pAmount * pAffineTP.z;
      // }

      pVarTP.doHide = false;

    }
    else { // if didn't do inversion, check to see if should hide
      pVarTP.x += xin;
      pVarTP.y += yin;
      if (pContext.isPreserveZCoordinate()) {
        pVarTP.z += pAmount * zin;
      }
      pVarTP.doHide = hide_uninverted;

    }
  }
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    ring_min = ring_min_ratio * r;
    ring_max = ring_max_ratio * r;
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
    return new Object[] { r, a, b, c, 
      zmode ? 1 : 0, 
      inversion_mode, hide_uninverted ? 1 : 0, ring_min_ratio, ring_max_ratio, p, draw_circle, guides_enabled ? 1 : 0 };
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
    else if (PARAM_ZORIGIN.equalsIgnoreCase(pName)) {
      c = pValue;
    }
    else if (PARAM_ZMODE.equalsIgnoreCase(pName)) {
      zmode = (pValue == 1) ? true : false;
    }
    else if (PARAM_RING_MIN.equalsIgnoreCase(pName)) {
      ring_min_ratio = pValue;
    }
    else if (PARAM_RING_MAX.equalsIgnoreCase(pName)) {
      ring_max_ratio = pValue;
    }
    else if (PARAM_P.equalsIgnoreCase(pName)) {
      p = pValue;
    }
    else if (PARAM_INVERSION_MODE.equalsIgnoreCase(pName)) {
      inversion_mode = (int)pValue;
    }
    else if (PARAM_HIDE_UNINVERTED.equalsIgnoreCase(pName)) {
      hide_uninverted = (pValue == 1) ? true : false;
    }
    else if (PARAM_DRAW_CIRCLE.equalsIgnoreCase(pName)) {
      draw_circle = pValue;
    }
    else if (PARAM_GUIDES_ENABLED.equalsIgnoreCase(pName)) {
      guides_enabled = (pValue == 1) ? true : false;
    }
    else
      throw new IllegalArgumentException(pName);
  }

}
