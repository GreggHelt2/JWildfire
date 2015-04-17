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

import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.fabs;
import static org.jwildfire.base.mathlib.MathLib.fmod;
import static org.jwildfire.base.mathlib.MathLib.sin;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

/*
        DC-Cylinder Apophysis Plugin
        This is a combo plugin based on the Rod Plugin by Branden Brown, a.k.a. zephyrtronium
        and the DC_linear Plugin by Georg Kiehne a.k.a. Xyruso2
*/
public class DCCylinderFunc extends VariationFunc {

  private static final long serialVersionUID = 1L;

  private static final String PARAM_OFFSET = "offset";
  private static final String PARAM_ANGLE = "angle";
  private static final String PARAM_SCALE = "scale";
  private static final String PARAM_X = "x";
  private static final String PARAM_Y = "y";
  private static final String PARAM_BLUR = "blur";

  private static final String[] paramNames = { PARAM_OFFSET, PARAM_ANGLE, PARAM_SCALE,
                                               PARAM_X, PARAM_Y, PARAM_BLUR  };

  // adjustable params
  private double offset = 0.0;
  private double angle = 0.0;
  private double scale = 0.5;
  private double x = 0.125;
  private double y = 0.125;
  private double blur = 1.0;

  // internal vars
  private int n;
  private double ldcs;
  private double ldca;
  private double[] rndr = new double[4];


  /*
        double r, sr, cr;
     
        fsincos(random01() * M_2PI, &sr, &cr);
        r = VAR(cyl_blur) * (VAR(r)[0] + VAR(r)[1] + VAR(r)[2] + VAR(r)[3] - 2.0);
        VAR(r)[VAR(n)] = random01();
        VAR(n) = VAR(n) + 1 & 3;
       
        FPx += VVAR * sin(FTx + r * sr)* VAR(cyl_x);
        FPy += r + FTy * VAR(cyl_y);
        FPz += VVAR * cos(FTx + r * cr);
        
         double c, s; fsincos(VAR(dc_cyl_angle), &s, &c);
         TC   = fmod( fabs( 0.5 * (VAR(ldcs) * ((c * FPx + s * FPy + VAR(dc_cyl_offset))) + 1.0) ), 1.0 );
  */
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    /* dc_linear by Xyrus02, http://apophysis-7x.org/extensions */
    double rnd_theta = pContext.random() * M_2PI;
    double sr = sin(rnd_theta);
    double cr = cos(rnd_theta);
    double r = blur * (rndr[0] + rndr[1] + rndr[2] + rndr[3] - 2.0);
    rndr[n] = pContext.random();
    n = n + 1 & 3;
    
    
    pVarTP.x += pAmount * sin(pAffineTP.x + (r * sr)) * x;
    pVarTP.y += r + (pAffineTP.y * y);
    pVarTP.z += pAmount * cos(pAffineTP.x + (r * cr));
    
    double sa = sin(angle);
    double ca = cos(angle);
    // pVarTP.color = fmod(fabs(0.5 * (ldcs * (((ca * pVarTP.x) + (sa * pVarTP.y) + offset)) + 1.0)), 1.0);
    pVarTP.rgbColor = true;
    pVarTP.redColor = (ca * pVarTP.x + offset) * 1000;
    pVarTP.greenColor = 0;
    pVarTP.blueColor = (sa * pVarTP.y + offset) * 1000;
    // pAffineTP.color = fmod(fabs(0.5 * (ldcs * (((ca * pVarTP.x) + (sa * pVarTP.y) + offset)) + 1.0)), 1.0);
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { offset, angle, scale, x, y, blur };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_OFFSET.equalsIgnoreCase(pName))
      offset = pValue;
    else if (PARAM_ANGLE.equalsIgnoreCase(pName))
      angle = pValue;
    else if (PARAM_SCALE.equalsIgnoreCase(pName))
      scale = pValue;
    else if (PARAM_X.equalsIgnoreCase(pName))
      x = pValue;
    else if (PARAM_Y.equalsIgnoreCase(pName))
      y = pValue;
    else if (PARAM_BLUR.equalsIgnoreCase(pName))
      blur = pValue;
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "dc_cylinder";
  }



  /* not using any of the random() generation, just a cut&paste holdover from rod plugin?
      int PluginVarPrepare(Variation* vp)
    {
        VAR(n) = 0;
        VAR(r)[0] = random01();
        VAR(r)[1] = random01();
        VAR(r)[2] = random01();
        VAR(r)[3] = random01();
        VAR(ldcs) = 1.0 / (VAR(dc_cyl_scale) == 0.0 ? 10E-6 : VAR(dc_cyl_scale));
    VAR(ldca) = VAR(dc_cyl_offset) * M_PI;
     
        return TRUE;
    }
  */


  
  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    n = 0;
    rndr[0] = pContext.random();
    rndr[1] = pContext.random();
    rndr[2] = pContext.random();
    rndr[3] = pContext.random();
    
    ldcs = 1.0 / (scale == 0.0 ? 10E-6 : scale);
    ldca = offset * M_PI;
  }

}
