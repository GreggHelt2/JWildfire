package org.jwildfire.create.tina.variation;

import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.atan2;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.base.Layer;

/*
 *  Double Helix, using B-DNA double helix measurements 
 *     to set ratios of diameter, helical pitch, minor groove size, major groove size relative to each other
 *
 */
public class DoubleHelixFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;

  private static final String PARAM_THICKNESS = "thickness";
  private static final String PARAM_TURNS = "turns";
  private static final String PARAM_XTHICKNESS = "xthickness";
  private static final String PARAM_YTHICKNESS = "ythickness";
  private static final String PARAM_ZTHICKNESS = "zthickness";

  private static final String[] paramNames = { PARAM_TURNS, PARAM_THICKNESS, PARAM_XTHICKNESS, PARAM_YTHICKNESS, PARAM_ZTHICKNESS };

  private double turns = 2;
  private double thickness = 0.3;
  private double xthickness = 0;
  private double ythickness = 0;
  private double zthickness = 1;

  // defaults are DNA B-structure in angstroms
  private double diameter = 20;
  private double pitch = 33.2;
  private double radius;
  
  private double phi = 1.6180339887;
  private double goldenRatio = phi;
  private double pitchToDiameterRatio;
  private double pitchToRadiusRatio;
  private double pitchScaled;
//  private double helixRadiusScaled;

  // major groove and minor groove size should add up to helixPitch?
  //    measurement I've found don't quite add up exactly
  //    but really just getting the relative sizes of major and minor groove
  //    will later calculate as:
  //        major_groove_size = (majorGroove/(majorGroove + minorGroove)) * helixPitch
  //        minor_groove_size = (minorrGroove/(majorGroove + minorGroove)) * helixPitch
  private double majorGroove = 22;   
  private double minorGroove = 12;
  private double grooveRatio;
  private double minorGroovePercent;
  private double majorGroovePercent;

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    // for B-DNA, this is close to Golden Ratio?
    // radius & pitch from different sources:
    //      32.1 to 33.7, ~20(average), ratio = ~1.6-1.68  (Textbook: Molecular and Cellular Biology, S.L.Wolfe)
    //      34, 19,   ratio = 1.79      (UC Davis course: http://biowiki.ucdavis.edu/Genetics/Unit_I%3A_Genes,_Nucleic_Acids,_Genomes_and_Chromosomes/2%3A_Structures_of_nucleic_acids/B-Form,_A-Form,_Z-Form_of_DNA)
    //      34, 20,   ratio = 1.7       (2002 paper: http://people.bu.edu/mfk/restricted566/dnastructure.pdf)
    //      33.2, 20, ratio = 1.66      (Wikipedia: http://en.wikipedia.org/wiki/Nucleic_acid_double_helix)
    //      34, 23.7, ratio = 1.43      (UIC course: http://tigger.uic.edu/classes/phys/phys461/phys450/ANJUM04/)
    //      Phi (golden ratio): 1.62
    //  so it is certainly within range, and the pitch definitely varies anyway depending on base pair composition (diameter as well?)
    radius = diameter / 2;
    pitchToRadiusRatio = pitch / radius;
    pitchScaled = pitchToRadiusRatio;
    grooveRatio = minorGroove/majorGroove;
    majorGroovePercent = majorGroove/(majorGroove + minorGroove);
    minorGroovePercent = minorGroove/(majorGroove + minorGroove);
  }
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    double theta = atan2(pAffineTP.y, pAffineTP.x);
    //    double pointR = pAffineTP.getPrecalcSqrt();
    double helix = pContext.random();
    double hoffset = (pContext.random()-0.5) * thickness;
    double rotation = 0;
    if (helix > 0.5)  {
      // shift to second helix
      rotation = minorGroovePercent * 2 * M_PI;
    }
    pVarTP.x += (pAmount + (xthickness * hoffset)) * cos((turns * theta) + rotation);
    pVarTP.y += (pAmount + (ythickness * hoffset)) * sin((turns * theta) + rotation);

    pVarTP.z += (pAmount * (pitchScaled / (2 * M_PI)) * turns * theta) + (hoffset * zthickness);
  }

  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { turns, thickness, xthickness, ythickness, zthickness };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_TURNS.equalsIgnoreCase(pName)) {
      turns = pValue;      
    }
    else if (PARAM_THICKNESS.equalsIgnoreCase(pName)) {
      thickness = pValue;
    }
    else if (PARAM_XTHICKNESS.equalsIgnoreCase(pName)) {
      xthickness = pValue;
    }
    else if (PARAM_YTHICKNESS.equalsIgnoreCase(pName)) {
      ythickness = pValue;
    }
    else if (PARAM_ZTHICKNESS.equalsIgnoreCase(pName)) {
      zthickness = pValue;
    }
    else {
      throw new IllegalArgumentException(pName);
    }
  }

  @Override
  public String getName() {
    return "doublehelix";
  }

}
