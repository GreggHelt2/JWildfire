package org.jwildfire.create.tina.variation;

import static java.lang.Math.abs;
import java.math.BigInteger;
import static org.jwildfire.base.mathlib.MathLib.EPSILON;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.atan2;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.fabs;
import static org.jwildfire.base.mathlib.MathLib.floor;
import static org.jwildfire.base.mathlib.MathLib.sqr;
import static org.jwildfire.base.mathlib.MathLib.pow;
import static org.jwildfire.base.mathlib.MathLib.sqrt;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

/**
 * InversionFunc, a variation for inversion geometry transformations
 * Supports standard circle inversion
 * Includes x0, y0, z0 params for specifying origin translationa specific radius parameter
 *     (efffectively replaces pre- and post- transform coefficients 
 *      so don't need to keep them as mirrors of each other)
 *  also includes "draw_circles" param
 *      if 0 > draw_circles < 1, then that fraction of incoming points is used to 
 *      draw circle of inversion rather than doing the actual inversion >
 *  In addition to standard circle inversion, can also be used for p-circle inversion as 
 *      described by Ramirez et al in "Generating Fractal Patterns by Using P-Circle Inversion" (2015)
 */
public class Inversion3dFunc extends VariationFunc implements Guides {
  private static final long serialVersionUID = 1L;

  public static final String PARAM_XORIGIN = "xorigin";
  public static final String PARAM_YORIGIN = "yorigin";
  public static final String PARAM_ZORIGIN = "zorigin";
  public static final String PARAM_ROTATION = "rotation (pi * n radians)";
  public static final String PARAM_SCALE= "scale";
  public static final String PARAM_SHAPE = "shape";
  public static final String PARAM_A = "a";
  public static final String PARAM_B = "b";
  public static final String PARAM_C = "c";
  public static final String PARAM_D = "d";
  public static final String PARAM_E = "e";
  public static final String PARAM_F = "f";
  
  public static final String PARAM_ZMODE = "zmode";
  public static final String PARAM_INVERSION_MODE = "imode";
  public static final String PARAM_HIDE_UNINVERTED = "hide_uninverted";
  // public static final String PARAM_RING_MIN = "ring_min";
  // public static final String PARAM_RING_MAX = "ring_max";
  public static final String PARAM_P= "p";
  public static final String PARAM_DRAW_CIRCLE = "draw_circle";
  public static final String PARAM_GUIDES_ENABLED = "guides_enabled";
  public static final String PARAM_PASSTHROUGH = "passthrough";
  
  private static final String[] paramNames = { 
    PARAM_SCALE, PARAM_ROTATION, 
    PARAM_SHAPE, 
    PARAM_ZMODE, 
    PARAM_INVERSION_MODE, PARAM_HIDE_UNINVERTED, 
    // PARAM_RING_MIN, PARAM_RING_MAX, 
    PARAM_P, PARAM_DRAW_CIRCLE, PARAM_PASSTHROUGH, PARAM_GUIDES_ENABLED, 
    PARAM_A, PARAM_B, PARAM_C, PARAM_D, PARAM_E, PARAM_F, 
    PARAM_XORIGIN, PARAM_YORIGIN, PARAM_ZORIGIN, 
  };
  
  public static int STANDARD = 0;
  public static int EXTERNAL_INVERSION_ONLY = 1;
  public static int INTERNAL_INVERSION_ONLY = 2;
  // outer inversion is like external inversion, except 
  //    leaves unchanged any  "external" regions that are distance < radius (inner hole)
  // public static int EXTERNAL_RING_INVERSION_ONLY = 3;
  // public static int INTERNAL_RING_INVERSION_ONLY = 4;
  
 /*
  public static int CIRCLE = 0;
  public static int ELLIPSE = 1;
  public static int HYPERBOLA = 2;
  public static int REGULAR_POLYGON = 3;
  public static int RHODONEA = 4;
  public static int SUPERSHAPE = 5;
  public static int WELDED_CIRCLES = 6;
  public static int RING1 = 7;
  public static int RING2 = 8;
  */
  public static int SPHERE = 0;
  public static int ELLIPSOID = 1;
  public static int TETRAHEDRON = 2;
  public static int REGULAR_POLYHEDRA = 3;
  public static int RHODONEA_3D = 4;
  public static int SUPERSHAPE_3D = 5;
//  public static int WELDED_SPHERES = 16;
//  public static int TORUS1 = 17;
//  public static int TORUS2 = 18;
  
  ParametricShape shape;
  boolean draw_guides = false;
  double rotation_pi_fraction = 0;
  double shape_rotation_radians;
  
  public void setDrawGuides(boolean draw) {
    draw_guides = draw;
  }
  
  public boolean getDrawGuides() { return draw_guides; }
  
  /**
   *  only using z coordinate for specific modes
   */
  
  class SphericalPoint3D {
    public double x;
    public double y;
    public double z;
    public double r;
    // public double t;
  }

  
  // abstract class ParametricShape {
  abstract class ParametricShape {
    public SphericalPoint3D ptemp = new SphericalPoint3D();
    public abstract void getCurvePoint(double theta, double phi, SphericalPoint3D point);
    public abstract void getInversePoint(SphericalPoint3D pin, SphericalPoint3D pout);
    public double getThetaPeriod() { return M_2PI; }
    public double getPhiPeriod()  { return M_PI; }
    
    /** simple shapes will not have more than one intersection along a given ray out from shape center */
    public boolean simpleShape() { return true; }

    // tin ==> input theta
    // fin ==> input phi
    public SphericalPoint3D getCurvePoint(double tin, double fin) {
      SphericalPoint3D outpoint = new SphericalPoint3D();
      getCurvePoint(tin, fin, outpoint);
      return outpoint;
    }
    
    // find intersection nearest to point pIn1 of line (defined by point pIn1 and pIn2) and shape 
    //  
    //    shapepoint p to "center" of shape
    //    (where "centroid" is defined by the shape object, 
    //       and is meant to usually be the centroid of the shape)
    /*
    public double getClosestIntersect(PolarPoint2D pin1, PolarPoint2D pin2, PolarPoint2D pout) {
      return 0;
    }
    */
    

    public void getMaxCurvePoint(double tin, double pin, SphericalPoint3D point_out) {
        getCurvePoint(tin, pin, point_out);
    }
  }
  
  // still need a parameter for shape rotation about (x0, y0) ?
  
  class Sphere extends ParametricShape {
    double r;
    double r2;
    
    public Sphere() {
      r = scale;
      r2 = r * r;
    }

    @Override
    public void getCurvePoint(double tin, double fin, SphericalPoint3D point) {
      // for consistency with other shapes, though rotation shouldn't matter for the circle...
     // double t = tin - shape_rotation_radians;  
//       double f = fin - phi_rotation_radians;
      double t = tin;
      double f = fin;
      point.r = r;
      point.x = r * cos(t) * sin(f);
      point.y = r * sin(t) * sin(f);
      point.z = r * cos(f);
    }
    
    public void getInversePoint(SphericalPoint3D pin, SphericalPoint3D pout) {
      // double r2 = r * r;
      double denom = (pin.x * pin.x) + (pin.y * pin.y) + (pin.z * pin.z);

      pout.x = (r2 * pin.x) / denom;
      pout.y = (r2 * pin.y) / denom;
      pout.z = (r2 * pin.z) / denom;
    }
  }
  
   /**
   * uses param a, b as standard a, b params for an ellipse
   */
  class Ellipsoid extends ParametricShape {
    double u, v, w;
    double u2, v2, w2;
    double uvw2;

    public Ellipsoid() {
      u = a;
      v = b;
      w = c;
      u2 = u * u;
      v2 = v * v;
      w2 = w * w;
      uvw2 = u2 * v2 * w2;
    }
    
    @Override
    public void getCurvePoint(double tin, double fin, SphericalPoint3D point) {
      // t is input angle modified by shape rotation
      double t = tin - shape_rotation_radians;
      double r = (a * b) / sqrt(sqr(b*cos(t)) + sqr(a*sin(t)));
      r *= scale;
      point.r = r;
      point.x = r * cos(tin);
      point.y = r * sin(tin);
    }
    
    public void getInversePoint(SphericalPoint3D pin, SphericalPoint3D pout) {
      double denom = (v2 * w2 * pin.x * pin.x) + (u2 * w2 * pin.y * pin.y) + (u2 * v2 * pin.z * pin.z);

      pout.x = (uvw2 * pin.x) / denom;
      pout.y = (uvw2 * pin.y) / denom;
      pout.z = (uvw2 * pin.z) / denom;
    }
  }
  
  /**
   *  uses parameter a as number of sides for polygon
   * 
   */
  /*
  class RegularPolygon extends ParametricShape {
    
    @Override 
    public double getPeriod() { return M_2PI; }
    
    @Override
    // parametric polygon equation derived from: 
    //    http://math.stackexchange.com/questions/41940/is-there-an-equation-to-describe-regular-polygons
    //    http://www.geogebra.org/m/157867
    public void getCurvePoint(double tin, SphericalPoint3D point) {
      double theta = abs((tin - shape_rotation_radians) % M_2PI);
      // double t = theta - shape_rotation_radians;
      double n = Math.floor(a);
      // double r = cos(M_PI/n) / cos(t%(M_2PI/n) - M_PI/n);
      double r = cos(M_PI/n) / cos(theta%(M_2PI/n) - M_PI/n);
      r *= scale;
      point.r = r;
      point.x = r * cos(tin);
      point.y = r * sin(tin);
    }
  }
  */
  
  /*
   * uses param a, b as standard a, b params for a hyperbola
   */
  /*
  class Hyperbola extends ParametricShape {
    
    @Override 
    public double getPeriod() { return M_2PI; }
    
    @Override
    public void getCurvePoint(double tin, SphericalPoint3D point) {
      double t = tin - shape_rotation_radians;
      double r2 = (a * a * b * b) / ((b * b * cos(t) * cos(t)) - (a * a * sin(t) * sin(t)));
      double r = sqrt(r2);
      r *= scale;
      point.r = r;
      point.x = r * cos(tin);
      point.y = r * sin(tin);
    }
  }
  */
  
  /**
   * 
   */
  /*
  class SuperShape extends ParametricShape {
    @Override
    public double getPeriod() {
      return M_2PI;
    }

    @Override
    public void getCurvePoint(double tin, SphericalPoint3D point) {
      double t = tin - shape_rotation_radians;

      // mapping params (a,b,c,d,e,f) to 
      //  naming convention of supershape (a,b,c,n1,n2,n3)
      // a = a param
      // b = b param
      double m = c;
      double n1 = d;
      double n2 = e;
      double n3 = f;
      
      double r = pow(
              (pow( fabs( (cos(m * t / 4))/a), n2) +
                      pow( fabs( (sin(m * t / 4))/b), n3)),
              (-1/n1));
      r *= scale;
      point.r = r;
      point.x = r * cos(tin);
      point.y = r * sin(tin);
    }
  }
  */
  
  // public boolean DESIGN_GUIDE_MODE = false;
  double scale = 1;
  double x0 = 0;
  double y0 = 0;
  double z0 = 0;
  double a = 1;
  double b = 1;
  double c = 0;
  double d = 0;
  double e = 0; 
  double f = 0;
  int shape_mode = SPHERE;
  // if pass_through != 0, then allow (pass_through) fraction of points to pass through unaltered
  double passthrough = 0;
  
  double p = 2;
  // double ring_min_ratio = 0;
  // double ring_max_ratio = 1;
  // double ring_min;
  // double ring_max;
  double draw_shape = 0;
  boolean guides_enabled = true;
  int inversion_mode = STANDARD;
  boolean hide_uninverted = false;
  boolean zmode = false;

  SphericalPoint3D curve_point = new SphericalPoint3D();
  SphericalPoint3D inpoint = new SphericalPoint3D();
  SphericalPoint3D outpoint = new SphericalPoint3D();

  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double xin = pAffineTP.x;
    double yin = pAffineTP.y;
    double zin = pAffineTP.z;
    inpoint.x = xin;
    inpoint.y = yin;
    inpoint.z = zin;

    double iscale;
    if (draw_guides && guides_enabled) { 
      double rnd = pContext.random();
        // double theta = rnd * shape.getPeriod();
      double theta = rnd * M_2PI;
       //  if (zmode) {
          double split = pContext.random() * 3.1;
          if (split < 1) {
            pVarTP.x += scale * cos(theta);
            pVarTP.y += scale * sin(theta);
           // pVarTP.z = z0;
           // pVarTP.z += zin;
          }
          else if (split < 2) {
            pVarTP.x += scale * cos(theta);
            // pVarTP.y = y0;
            // pVarTP.y += yin;
            pVarTP.z += scale * sin(theta);
          }
          else if (split < 3) {
            // pVarTP.x = x0;
            // pVarTP.x += xin;
            pVarTP.y += scale * sin(theta);
            pVarTP.z += scale * cos(theta);
          }
          else {
            // pVarTP.x += x0;
            // pVarTP.x += xin;
            // pVarTP.y += yin;
            // pVarTP.z += zin;
          }
      //  }
    /*
    else {
          double split = pContext.random() * 1.1;
          if (split < 1) {
            shape.getCurvePoint(theta, curve_point);
            pVarTP.x += curve_point.x;
            pVarTP.y += curve_point.y;
            //pVarTP.z += curve_point.z;
            // pVarTP.z += z0;
          }
          else {  // draw point at center of shape
            // pVarTP.x += x0;
            // pVarTP.y += y0;
            // pVarTP.z += z0;
          }
        }
    */
      return;
    }
    
   /* if (draw_shape > 0) {
      double rnd = pContext.random();
      if (rnd < draw_shape) {
        double theta = pContext.random() * shape.getPeriod();
        shape.getCurvePoint(theta, curve_point);
        pVarTP.x += curve_point.x;
        pVarTP.y += curve_point.y;
        return;
      }
    }
    */
    if (passthrough > 0) {
      double rnd = pContext.random(); 
      if (rnd < passthrough) {
        pVarTP.x += xin;
        pVarTP.y += yin;
        pVarTP.z += zin;
        return;
      }
    }
    // to do generalized inversion of input point P, 
    // need two other points:
    // O, the origin of inversion
    // S, the intersection of the line OP and the surface of the shape
    // then output point P' = O + (d1(O,B)^2/d2(O,P)^2) * (P - O)
    // where d1 and d2 are distance metric functions 
    // double tin = atan2(yin, xin);
    double rin = sqrt(xin*xin + yin*yin + zin*zin);
    boolean do_inversion;
    // shape.getCurvePoint(tin, curve_point);
    // shape.getMaxCurvePoint(tin, curve_point);
    // double rcurve = sqrt(curve_point.x * curve_point.x + curve_point.y * curve_point.y);
    double rcurve = curve_point.r;
    if (inversion_mode == EXTERNAL_INVERSION_ONLY) {
      // only do inversion if input point is outside of circle
      do_inversion = rin > rcurve;
    }
    else if (inversion_mode == INTERNAL_INVERSION_ONLY) {
      // only do inversion if input point is inside of circle
      do_inversion = rin < rcurve;
    }
    else { // default to STANDARD mode ==> always do inversion
      do_inversion = true;
    }
    if (do_inversion) {
      shape.getInversePoint(inpoint, outpoint);
      /*
      if (p == 2) {
        iscale = (rcurve * rcurve) / (sqr(xin) + sqr(yin));
      }
      else {
        iscale = (rcurve * rcurve)/ pow( (pow(abs(xin),p) + pow(abs(yin),p)), 2.0/p);
      }
      double xout = iscale * xin;
      double yout = iscale * yin;
      double zout = iscale * zin;
      */
      pVarTP.x += pAmount * outpoint.x;
      pVarTP.y += pAmount * outpoint.y;
      pVarTP.z += pAmount * outpoint.z;
      pVarTP.doHide = false;
    }
    else { // if didn't do inversion, check to see if should hide
      pVarTP.x += xin;
      pVarTP.y += yin;
      pVarTP.z += zin;
      pVarTP.doHide = hide_uninverted;
    }
  }
  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    shape_rotation_radians = M_PI * rotation_pi_fraction;
    // ring_min = ring_min_ratio * r;
    // ring_max = ring_max_ratio * r;
    if (shape_mode == SPHERE) {
      shape = new Sphere();
    }
    else if (shape_mode == ELLIPSOID) {
      shape = new Ellipsoid();
    }
/*
    else if (shape_mode == HYPERBOLA) {
      shape = new Hyperbola();
    }
    else if (shape_mode == REGULAR_POLYGON) {
      shape = new RegularPolygon();
    }
    else if (shape_mode == RHODONEA) {
      shape = new Rhodonea();
    }
    else if (shape_mode == SUPERSHAPE) {
      shape = new SuperShape();
    }
    */
  }

  @Override
  public String getName() {
    return "inversion_3D";
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { 
      scale, 
      rotation_pi_fraction, 
      shape_mode, 
      zmode ? 1 : 0, 
      inversion_mode, hide_uninverted ? 1 : 0, 
      // ring_min_ratio, ring_max_ratio, 
      p, draw_shape, passthrough, guides_enabled ? 1 : 0, 
      a, b, c, d, e, f, 
      x0, y0, z0, 
    };
    
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_SCALE.equalsIgnoreCase(pName)) {
      scale = pValue;
    }
    else if (PARAM_ROTATION.equalsIgnoreCase(pName)) {
      rotation_pi_fraction = pValue;
    }
    else if (PARAM_SHAPE.equalsIgnoreCase(pName)) {
      shape_mode = (int)floor(pValue);
    }
        else if (PARAM_ZMODE.equalsIgnoreCase(pName)) {
      zmode = (pValue == 1) ? true : false;
    }
   /*
    else if (PARAM_RING_MIN.equalsIgnoreCase(pName)) {
      ring_min_ratio = pValue;
    }
    else if (PARAM_RING_MAX.equalsIgnoreCase(pName)) {
      ring_max_ratio = pValue;
    }
    */
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
      draw_shape = pValue;
    }
    else if (PARAM_PASSTHROUGH.equalsIgnoreCase(pName)) {
      passthrough = pValue;
    }
    else if (PARAM_GUIDES_ENABLED.equalsIgnoreCase(pName)) {
      guides_enabled = (pValue == 1) ? true : false;
    }
    else if (PARAM_A.equalsIgnoreCase(pName)) {
      a = pValue;
    }
    else if (PARAM_B.equalsIgnoreCase(pName)) {
      b = pValue;
    }
    else if (PARAM_C.equalsIgnoreCase(pName)) {
      c = pValue;
    }
    else if (PARAM_D.equalsIgnoreCase(pName)) {
      d = pValue;
    }
    else if (PARAM_E.equalsIgnoreCase(pName)) {
      e = pValue;
    }
    else if (PARAM_F.equalsIgnoreCase(pName)) {
      f = pValue;
    }
    else if (PARAM_XORIGIN.equalsIgnoreCase(pName)) {
      x0 = pValue;
    }
    else if (PARAM_YORIGIN.equalsIgnoreCase(pName)) {
      y0 = pValue;
    }
    else if (PARAM_ZORIGIN.equalsIgnoreCase(pName)) {
      z0 = pValue;
    }
    else
      throw new IllegalArgumentException(pName);
  }

}
