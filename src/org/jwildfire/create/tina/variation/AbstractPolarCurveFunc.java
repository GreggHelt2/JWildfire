package org.jwildfire.create.tina.variation;

import static java.lang.Math.abs;
import java.util.ArrayList;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.sqrt;
import static org.jwildfire.base.mathlib.MathLib.atan2;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.floor;
import static org.jwildfire.base.mathlib.MathLib.pow;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.sqr;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;

// private enum PointState {INSIDE, OUTSIDE, INISH, OUTISH, ON }  // "state" of point relative to curve

public abstract class AbstractPolarCurveFunc extends VariationFunc {
  protected boolean DEBUG = false;

  public enum RenderMode { ONCURVE, 
                           UNCHANGED, 
                           HIDE,
                           MIRROR_SWAP, 
                           CURVE_SWAP, 
                           WHIRL, 
                           SCALE, 
                           
                                                      
  }
  
  public enum PointState { INSIDE, OUTSIDE, INISH, OUTISH, ON }  // "state" of point relative to curve
  public enum CurveRadiusMode { THETA, RAW_SAMPLING_BIN, INTERPOLATED_SAMPLING_BIN }
  public enum PointRadiusMode { MODIFIED, RAW }
  public enum InsideOutsideRule { EVEN_ODD, MODIFIED_EVEN_ODD, MODIFIED_EVEN_ODD_INISH_OUTISH }
  public enum MergeMode { AUTO, 
                          INNERS_OUTERS, // inner* controls both inner and inish, outer* controls both outer and outish
                          INNERS_OUTER_OUTISH,  // inner* controls both inner and inish
                          INNERSISH_OUTER, // inner* controls inner, inish, outish
                          ALL // inner* controls all (inner/outer/inish/outish)
  }

  protected static final String PARAM_CURVE_SCALE = "curve_scale";
  protected static final String PARAM_MODE_MERGING = "mode_merging";
  protected static final String PARAM_INNER_MODE = "inner_mode";
  protected static final String PARAM_OUTER_MODE = "outer_mode";
  // protected static final String PARAM_INISH_MODE = "inish_mode";
  protected static final String PARAM_OUTISH_MODE = "outish_mode";
  protected static final String PARAM_INNER_SPREAD = "inner_spread";
  protected static final String PARAM_OUTER_SPREAD = "outer_spread";
  // protected static final String PARAM_INISH_SPREAD = "inish_spread";
  protected static final String PARAM_OUTISH_SPREAD = "outish_spread";
  protected static final String PARAM_INNER_SPREAD_RATIO = "inner_spread_ratio";
  protected static final String PARAM_OUTER_SPREAD_RATIO = "outer_spread_ratio";
  // protected static final String PARAM_INISH_SPREAD_RATIO = "inish_spread_ratio";
  protected static final String PARAM_OUTISH_SPREAD_RATIO = "outish_spread_ratio";
  protected static final String PARAM_SPREAD_SPLIT = "spread_split";
  protected static final String PARAM_FILL = "fill";
  protected static final String PARAM_CYCLES = "cycles";
  protected static final String PARAM_CURVE_RADIUS_MODE = "curve_radius_mode";
  protected static final String PARAM_CYCLE_OFFSET = "cycle_offset";
  // metacycles (and metacycle_expansion) are only in effect when cycles = 0 (automatic cycle calculations in effect), 
  //   and cycles to close can be determined (cycles_to_close > 0)
  protected static final String PARAM_METACYCLES = "metacycles";
  protected static final String PARAM_METACYCLE_OFFSET = "metacycle_offset";
  protected static final String PARAM_METACYCLE_SCALE = "metacycle_scale";
 // protected static final String PARAM_METACYCLE_ROTATION = "metacycle_rotation";
  
  protected static final String[] paramNames = { 
                                               PARAM_CURVE_SCALE, 
                                               PARAM_MODE_MERGING, 
                                               PARAM_INNER_MODE, PARAM_OUTER_MODE, PARAM_OUTISH_MODE, 
                                               PARAM_INNER_SPREAD, PARAM_OUTER_SPREAD, PARAM_OUTISH_SPREAD, 
                                               PARAM_INNER_SPREAD_RATIO, PARAM_OUTER_SPREAD_RATIO, PARAM_OUTISH_SPREAD_RATIO, 
                                               PARAM_SPREAD_SPLIT,
                                               PARAM_CYCLES, PARAM_CYCLE_OFFSET, 
                                               PARAM_FILL, PARAM_CURVE_RADIUS_MODE, 
                                               PARAM_METACYCLES, PARAM_METACYCLE_OFFSET, PARAM_METACYCLE_SCALE };
//                                               PARAM_METACYCLE_ROTATION }; 

  protected double curve_scale = 1;
  // mode merging:
  //    0 is auto, code attempts to merge based on chosen inner/outer modes
  //    1 => inner/inish merged, outer/outish merged
  //    2 => inner/inish merged
  protected int mode_merging = 0; 
  protected int inner_mode = 1;
  protected int outer_mode = 1;
  protected int outish_mode = 1;

  protected double inner_spread = 0; // deform based on original x/y
  protected double outer_spread = 0; // deform based on original x/y
  protected double outish_spread = 0;
  
  protected double inner_spread_ratio = 1; // how much inner_spread applies to x relative to y
  protected double outer_spread_ratio = 1; // how much outer_spread applies to x relative to y
  protected double outish_spread_ratio = 1;
  
  // protected int inish_mode = 1;
  // protected double inish_spread = 0;
  // protected double inish_spread_ratio = 1;

  
  protected double spread_split = 1;
  protected double fill = 0;
  protected double cycles_param;
  protected double curve_radius_mode = 2;
  protected double point_radius_mode;
  protected double in_out_mode;
  
  protected double cycles;  // 1 cycle = 2*PI
  // protected double radians; // = 2*PI*cycles
  protected double cycle_length = M_2PI; // 2(PI)
  protected double cycles_to_close = 0; // 0 indicates unknown, -1 indicates curve will never close
  protected double cycle_offset = 0; // cycle offset (in cycles) for incoming points (rotate the cycle)
  protected double metacycles = 1; // if cycles is calculated automatically to close the curve, metacycles is number of times to loop over closed curve
  protected double metacycle_offset = 0.1; // P(m) = P * (1 + mOffset) * ((mScale)^m))  // cumulative offset for metacycles
  protected double metacycle_scale = 1.1;
  // protected double metacycle_rotation = 0; // additional (cumulative) theta offset (in cycles) for each metacycle (rotate the metacycle)

  // vars for determining inner/outer via even-odd rule
  int default_sample_count = 36000;
  int binCount = 720;
  ArrayList<ArrayList<Double>> theta_intersects = null;
  ArrayList<Double> unisects = new ArrayList<Double>(1);  // temp wrapper when needed for single curve radius
  XYZPoint pCurve = new XYZPoint();
  
  CurveRadiusMode curve_rmode = CurveRadiusMode.THETA;
  PointRadiusMode point_rmode = PointRadiusMode.MODIFIED;
  InsideOutsideRule inout_rule = InsideOutsideRule.MODIFIED_EVEN_ODD_INISH_OUTISH;
  MergeMode merger = MergeMode.INNERS_OUTERS;

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    // System.out.println("calling init for VariationFunc: " + this);
    if (curve_radius_mode == 0) {
      // auto calling of curve radius mode?
      // possibly mixed?
    }
    if (curve_radius_mode == 1) {
      curve_rmode = CurveRadiusMode.RAW_SAMPLING_BIN;
    }
    else if (curve_radius_mode == 2) {
      curve_rmode = CurveRadiusMode.THETA;
    }
    else if (curve_radius_mode == 3) {
      curve_rmode =  CurveRadiusMode.INTERPOLATED_SAMPLING_BIN;
    }
    else {
      // auto calling
    }
    unisects = new ArrayList<Double>(1);
    unisects.add(0, 0.0);
    recalcCycles();

    recalcCurveIntersects();
  }
  
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    // input
    // output
    double tin = atan2(pAffineTP.y, pAffineTP.x);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    double theta = cycles * tin;
    
    //    SimplePoint curvePoint = new SimplePoint(0.0, 0.0, 0.0);
    pCurve.clear();
    calcCurvePoint(pContext, theta, pCurve);
    renderByMode(pContext, pXForm, pAffineTP, pVarTP, pAmount, pCurve);
    
  }
  
  /*
  *  calculates number of cycles if not specified (cycles_param == 0);
  *  subclasses need to override this to calculate cycles and cycles_to_close (if it can be calculated)
  *  should end subclass method by calling super.recalcCycles(), to allow handling of metacycles
  */
  public void recalcCycles() {
    // handle metacycles (previously was done in subclass.recalcCycles(), but moved 
    //    to superclass so subclass doesn't have to have any knowledge of metacycles
    if (cycles_param == 0 && cycles_to_close > 0) {
      // auto cycle setting -- add in metacycles
      //   cycles _should_ be set to cycles_to_close in subclass.recalcCycles(), so 
      //   this should be equivalent to (cycles_to_close * metacycles), 
      //   but using cycles just in case a subclass has different ideas

      //cycles = cycles * metacycles; 
      cycles = cycles_to_close * metacycles;
    }
    else {
      // cycles to close is unknown, or curve doesn't close
      // so leave cycles alone as was either specified in cycles_param or calculated in subclass.recalcCycles()
    }
    if (DEBUG) {
      System.out.println("cycles to close: " + cycles_to_close);
      System.out.println("metacycles: " + metacycles);
      System.out.println("total cycles: " + cycles);
    }
  }

  /* 
   *  actual calculation of points on curve -- called by transform() and recalcCurveIntersects()
   *  subclasses need to override this to add curve calculations, 
   *  should end subclass method by calling super.calcCurvePoint(), to allow handling of metacycles
  */
  public void calcCurvePoint(FlameTransformationContext pContext, double theta, XYZPoint pResult) {
    // pResult.clear();  NO, don't clear point! need incoming point from subclass
    if (cycles_param == 0 && cycles_to_close > 0 && metacycles != 1) {
      // double metacycle_count = floor((cycles * (tin + M_PI)) / (cycles_to_close * 2 * M_PI));
      double metacycle_count = floor((theta + (cycles * M_PI)) / (cycles_to_close * 2 * M_PI));
      if (DEBUG && metacycle_count < 0) {
        System.out.println("uhoh: theta = " + theta + ", cycles = " + cycles + ", cycles_to_close = " + cycles_to_close);
      }
      double metacycle_delta = metacycle_count * metacycle_offset;
      pResult.x = pResult.x * (1 + metacycle_delta);
      pResult.y = pResult.y * (1 + metacycle_delta);
      // z unchanged?
    }
  }
  
   public void recalcCurveIntersects() {
    // System.out.println("recalcing curves: " + this);
    theta_intersects = new ArrayList<ArrayList<Double>>(binCount);
    for (int i=0; i<binCount; i++) { 
      theta_intersects.add(new ArrayList<Double>());
    }
    ArrayList<Double> tsects;
    ArrayList<Double> prev_tsects = null;
    int firstbin = -1;
    int lastbin = -1;
    int sampleCount = default_sample_count;
    //if (cycles_param == 0 && cycles_to_close > 0 && metacycles > 1) {
    //  sampleCount = (int)(sampleCount * metacycles);
    //}

    for (int i=0; i<sampleCount; i++) {
//      double theta = ((double)i/(double)sampleCount) * cycles * M_2PI;
      double theta = (((double)i/(double)sampleCount) * cycles * M_2PI);

      // need to compensate to for theta adjusment in calcCurvePoint() ??
      // theta = theta - (cycles * M_PI);
      theta = (cycles * M_PI) - theta;
      
      // SimplePoint pcurve = calcFromTheta(theta);
      pCurve.clear();
      calcCurvePoint(null, theta, pCurve);
      
      double r = sqrt(pCurve.x * pCurve.x + pCurve.y * pCurve.y);
      double angle = atan2(pCurve.y, pCurve.x);
      
      int anglebin =  (int)Math.floor(((angle + M_PI)/M_2PI) * binCount);
     //  int anglebin =  (int)Math.floor(((tin + M_PI)/M_2PI) * binCount);
      if (i == 0) { 
        firstbin = anglebin;
        // System.out.println("anglebin at sample " + i + ": " + anglebin); 
      }
      if (i == (sampleCount-1)) { 
        lastbin = anglebin;
        // System.out.println("anglebin at sample " + i + ": " + anglebin); 
      }
      if (anglebin == binCount) { anglebin--; } // catching any possible cases where angle actually reaches max atan2
      tsects = theta_intersects.get(anglebin);

      // still rotating through same bin, merge results
      if (prev_tsects == tsects) {  
        // try ignoring for now -- should try averaging later?
      }
      else {
        tsects.add(r);  // autoboxing float r to Double object      
      }
      prev_tsects = tsects;
    }
    
    // special-casing of first and last anglebin if they are the same bin:
    ///   want to simulate rotating through same bin to merge "duplicate" intersections
    //    if first and last bin are same, would have merged results, so remove last one 
    if (firstbin > 0 && lastbin > 0 && firstbin == lastbin) {
      tsects = theta_intersects.get(firstbin);
      // remove last result (if start doing averaging, should remove but add to average for first result)
      tsects.remove(tsects.size()-1);
    }
  }
  

  //  public Point 

      /*
    // 1. Determine which inner/outer/outish calling method to use
    //    (so far all are radius-based, but which "radius"?)
    //     a) input point
    //         i) rin
    //         ii) raw_rin
    //     b) curve point(s)
    //         i) r
    //         ii) ray pre-calced bin approximation
    //         iii) ray bin approximation plus linear interpolation
    //         iv) equation solutions???
    //     c) caller
    //         i) input <, =, > max curve
    //         ii) input <, =, > min curve (same as (i) when only one intersection
    //         iii) even-odd rule (similar to (i) and (ii) when only one intersection
    //
    // 2. Make call for input point -- inner/outer/inish/(on?)
    //      INSIDE: point falls on or inside curve 
    //      OUTSIDE: point falls outside of curve 
    //      OUTISH: considered outside by standard even-odd algorithm, but still within max radius of curve at given theta
    //      INISH:  considered inside by standard even-odd algorithm, but not contiguous with point of origin)
    //                ==> not currently making inish calls, too problematic to distinguish from INSIDE
    //       ONCURVE: on the curve ==> not currently making on-the-curve calls
    //       
    // 3. Prep input for feeding to switch statement:
    //       using designated inner/outer/inish/(on?) for
    //       *_mode, *_spread, *_spread_ratio, etc.
    //
    // (3.5 -- insert probablistic mode determination?)
    //
    // 4. Transform input point to get intermediate point: big switch statement, switching on mode
    // 
    // 5. 
    // 

    // assume 1.a is always rin for now
    //
   */
   
  public void renderByMode(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount, XYZPoint curvePoint) {
    PointState pstate = PointState.INSIDE;  
    // input point
    double tin = atan2(pAffineTP.y, pAffineTP.x);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    // theta now only used in transform() method
    // double theta = cycles * (tin + (cycle_offset * M_2PI));
    double raw_rin = sqrt((pAffineTP.x  * pAffineTP.x) + (pAffineTP.y * pAffineTP.y));
    double rin = spread_split * raw_rin;
    double adjustedAmount = pAmount;
    
    double xin, yin;
    double rinx, riny;
    double rout;
    
    // input point mapped to curve
    double x = curvePoint.x;
    double y = curvePoint.y;
    // double t = curvePoint.getPrecalcAtanYX();
    double t = atan2(y, x);
    // double r = curvePoint.getPrecalcSqrt();
    double r = sqrt(x*x + y*y);
    ArrayList<Double> tsects;

    if (curve_rmode == CurveRadiusMode.RAW_SAMPLING_BIN) {
      int anglebin =  (int)Math.floor(((tin + M_PI)/M_2PI) * binCount);
      if (anglebin == binCount) {  // catching any possible cases where tin actually reaches max atan2
        anglebin--; 
      } 
      tsects = theta_intersects.get(anglebin);
    }
    else if (curve_rmode == CurveRadiusMode.THETA)  {
      tsects = unisects;
      tsects.set(0, r);
    }
    // else if (curve_rmode == CurveRadiusMode.INTERPOLATED_SAMPLING_BIN) {
      // same as RAW_SAMPLING_BIN, plus interpolation
    // }
    else {  // default to CurveRadiusMode.THETA ??
      tsects = unisects;
      tsects.set(0, r);
    }

    double rpoint;
    if (point_rmode == PointRadiusMode.RAW) {
      rpoint = raw_rin;
    }
    else if (point_rmode == PointRadiusMode.MODIFIED) {
      rpoint = rin;
    }
    else { // default to PointRadiusMode.MODIFIED ?
      rpoint = rin;
    }
    
    int shorter = 0;
    int longer = 0;
    double longest = Double.NEGATIVE_INFINITY;
    double shortest = Double.POSITIVE_INFINITY;
    double nearest = Double.NEGATIVE_INFINITY;
    double nearest_longer  = Double.POSITIVE_INFINITY;
    double nearest_shorter = Double.NEGATIVE_INFINITY;

    for (double rcurve : tsects) {
        //      if (ir <= raw_rin) { 
      if (rcurve <= rpoint) { 
        shorter++;
        shortest = Math.min(shortest, rcurve);
        if ((rpoint - rcurve) < (rpoint - nearest_shorter)) { nearest_shorter= rcurve; }
      }
      else { // rcurve > rpoint
        longer++; 
        longest = Math.max(longest, rcurve);
        if ((rcurve - rpoint) < (nearest_longer - rpoint)) { nearest_longer = rcurve; }
      }
    }
    if (longer == 0) { nearest = nearest_shorter; }
    else if (shorter == 0) { nearest = nearest_longer; }
    else if ((rpoint - nearest_shorter) < (nearest_longer - rpoint)) { nearest = nearest_shorter; }
    else { nearest = nearest_longer; }

    if (inout_rule == InsideOutsideRule.MODIFIED_EVEN_ODD_INISH_OUTISH) {
      if (longer == 0) { // definitely outside
        pstate = PointState.OUTSIDE;
      }
      else if (shorter == 0) { // definitely inside
        pstate = PointState.INSIDE;
      }
      else { 
        // either INISH or OUTISH, using same rules as MODIFIED_EVEN_ODD uses for INSIDE/OUTSIDE
        //  for now calling INISH as INSIDE -- encountered too many problems trying to distinguish between the two
        if (longer % 2 == 0) { // ray overlaps even number of intersections further away from origin than point
          if (longer == tsects.size()) {
            // pstate = PointState.INISH; // point is "inish" (assumes curve is closed and encloses 0?)
            pstate = PointState.INSIDE; // point is inside (assumes curve is closed and encloses 0?)
          }
          else if (longer == 2 && tsects.size() == 3) {
            pstate = PointState.OUTISH; // point is outside (assumes curve is closed and encloses 0?)
          }
          else { pstate = PointState.OUTISH;  } // point is weird? but usually outside
        }
        else {  // ray overlaps odd number of intersection further away from origin than point
          if (shorter == 2 && longer == 1) { // special handling for problematic case
            pstate = PointState.INSIDE;
          }
          else { 
            //pstate = PointState.INISH; 
            pstate = PointState.INSIDE; 
          }
        }
      }
    }
        
    else if (inout_rule == InsideOutsideRule.MODIFIED_EVEN_ODD) {
      /*
      *  use modified Even-Odd rule for inside/outside determination
      *  asssumes that curve is closed and encloses origin (0, 0)
      *  trying to handle cases where ray touches point on curve but does not actually cross
      *      (because of binning approximations, get this more often than would otherwise expect)
      */
      if (longer % 2 == 0) { // ray overlaps even number of intersections further away from origin than point
        if (longer == 0) { pstate = PointState.OUTSIDE; } // point is outside
        // else if (longer == 2 && tsects.size() == 2) {
        else if (longer == tsects.size()) {
          pstate = PointState.INSIDE; // point is inside (assumes curve is closed and encloses 0?)
        }
        else if (longer == 2 && tsects.size() == 3) {
          pstate = PointState.OUTSIDE; // point is outside (assumes curve is closed and encloses 0?)
        }
        else { pstate = PointState.OUTSIDE;  } // point is weird? but usually outside
      }
      else {  // ray overlaps odd number of intersection further away from origin than point
        pstate = PointState.INSIDE;
      }
    }
    else if (inout_rule == InsideOutsideRule.EVEN_ODD) {
      // use standard Even-Odd rule (well, more standard than the modified one above...):
      //    cast ray from origin through incoming point to infinity
      //    count how many times curve intersects ray further out than incoming point (longer)
      //    if number is odd then point is inside, if number is even then point is outside
      if (longer % 2 == 0) { pstate = PointState.OUTSIDE; } // point is outside
      else { pstate = PointState.INSIDE; } // point is inside
    }

      // 3. Prep input for feeding to switch statement:
      //       using designated inner/outer/inish/(on?) for
      //       *_mode, *_spread, *_spread_ratio, etc.
      double spread;
      double spread_ratio;
      int mode = 0;
      if (pstate == PointState.INSIDE)  {
        mode = inner_mode;
        spread = inner_spread;
        spread_ratio = inner_spread_ratio;
      }
      else if (pstate == PointState.OUTSIDE) {
        mode = outer_mode;
        spread = outer_spread;
        spread_ratio = outer_spread_ratio;
      }
      /*  commenting out INISH mode for now (combined with INSIDE)
      else if (pstate == PointState.INISH) {
        // using same params as INSIDE for now, introduce separate modes later
        // mode = inish_mode;
        // spread = inish_spread;
        // spread_ratio = inish_spread_ratio;
        mode = inish_mode;
        spread = inish_spread;
        spread_ratio = inish_spread_ratio;
      }
      */
      else if (pstate == PointState.OUTISH) {
        mode = outish_mode;
        spread = outish_spread;
        spread_ratio = outish_spread_ratio;
      }
      else if (pstate == PointState.ON) {
        // currently not using ON, so should never get here
        mode = 0;
        spread = 0;
        spread_ratio = 1;
      }
      else {
        // every point must be classified as being within one of the enum'd areas, so should never get here
        mode = 0;
        spread = 0;
        spread_ratio = 0;
      }
      
      // 4. pVarTP transformations, based on mode (for given INSIDE/OUTSIDE/OUTISH classification given above)
      switch (mode) {
        case 0: // ONCURVE -- no spread, place on curve
          pVarTP.x += adjustedAmount * x;
          pVarTP.y += adjustedAmount * y;
          break;
        case 1: // STRETCH2
          rinx = (rin * spread * spread_ratio) - (spread * spread_ratio) + 1;
          riny = (rin * spread) - spread + 1;
          pVarTP.x += adjustedAmount * rinx * x;
          pVarTP.y += adjustedAmount * riny * y;
          break;
        case 2: // STRETCH3
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += adjustedAmount * (x - (spread * spread_ratio * (x-xin)));
          pVarTP.y += adjustedAmount * (y - (spread * (y-yin)));
          break;
        case 3: // STRETCH4
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += adjustedAmount * (x - (spread * spread_ratio * xin));
          pVarTP.y += adjustedAmount * (y - (spread * yin));
          break;
        case 4: // STRETCH5
          rinx = (0.5 * rin) + (spread * spread_ratio);
          riny = (0.5 * rin) + spread;
          pVarTP.x += adjustedAmount * rinx * x;
          pVarTP.y += adjustedAmount * riny * y;
          break;
        case 5: // STRETCH6 -- same as mode 3, but without the sign modifications
          pVarTP.x += adjustedAmount * (x + (spread * spread_ratio * pAffineTP.x));
          pVarTP.y += adjustedAmount * (y + (spread * pAffineTP.y));
          break;
        case 6: // STRETCH7 -- similar to 3, different sign fiddling
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += adjustedAmount * (x + (spread * spread_ratio * xin));
          pVarTP.y += adjustedAmount * (y + (spread * yin));
          break;
        case 7: // STRETCH8 -- same as mode 6, but without the sign modifications
          pVarTP.x += adjustedAmount * (x - (spread * spread_ratio * pAffineTP.x));
          pVarTP.y += adjustedAmount * (y - (spread * pAffineTP.y));
          break;
        case 8: // STRETCH9
          pVarTP.x += adjustedAmount * rin * cos(t) * (spread * spread_ratio);
          pVarTP.y += adjustedAmount * rin * sin(t) * spread;
          break;
        case 9:  // STRETCH1
          // double rout = (rin * spread) + (r * (1 - spread));
          rout = r + ((rin - r) * spread);
          pVarTP.x += adjustedAmount * rout * cos(t);
          pVarTP.y += adjustedAmount * rout * sin(t);
          break;
        case 10:   // UNCHANGED -- leave in place
          // point is inside curve, leave in place
          pVarTP.x += adjustedAmount * pAffineTP.x;
          pVarTP.y += adjustedAmount * pAffineTP.y;
          break;
        case 11: // ONCURVE -- place on curve, same as mode 0 (after rearrangement?)
          // point is inside curve, place on curve
          pVarTP.x += adjustedAmount * x;
          pVarTP.y += adjustedAmount * y;
          break;
        case 12:  // ORIGIN_SWAP (swap around origin [0,0], inspired by RosoniFunc)
          pVarTP.x += adjustedAmount * pAffineTP.x * -1;
          pVarTP.y += adjustedAmount * pAffineTP.y * -1;
          break;
        case 13: // CURVE_XY_OFFSET -- (UNCHANGED + CURVE) combo
          pVarTP.x += adjustedAmount * (pAffineTP.x + x);
          pVarTP.y += adjustedAmount * (pAffineTP.y + y);
          break;
        case 14: // CURVE_RADIUS_OFFSET -- swap around intersect with longest radius, with: P' = C + P
          // only swap "outside" points that are still internal to overall curve
          // true by definition for inside points?), but just going for consistency across inner/outer modes
          if (longer > 0) { 
            double rx = longest * cos(tin);
            double ry = longest * sin(tin);
            pVarTP.x += adjustedAmount * (pAffineTP.x + rx);
            pVarTP.y += adjustedAmount * (pAffineTP.y + ry);
          }
          else { // place on curve 
            pVarTP.x += adjustedAmount * x;
            pVarTP.y += adjustedAmount * y;
          }
          // may want to add a radius clamp -- beyond clamp, don't swap?
          //   if (rin > clamp) { leave in place }
          break;
        case 15:  // CURVE_SWAP_RAW -- swap around intersect with longest radius of curve at angle tin, with P' = C + (C-P)
          // only differs from mode 16 by using raw_rin instead of rin
          // only swap "outside" points that are still internal to overall curve
          // true by definition for inside points?), but just going for consistency across inner/outer modes
          if (longer > 0) { 
             double rdiff = longest - raw_rin;
            rout = longest + rdiff;
            double rx= rout * cos(tin);
            double ry = rout * sin(tin);
            // or equivalently, (rx + rx - x, ry + ry - y) ??? (P' = I + (I-P))
            pVarTP.x += adjustedAmount * rx;
            pVarTP.y += adjustedAmount * ry;
          }
          else { // place on curve 
            pVarTP.x += adjustedAmount * x;
            pVarTP.y += adjustedAmount * y;
          }
          break;
        case 16:  // CURVE_SWAP swap around intersect with longest radius, with P' = C + (C-P)
          // only differs from mode 15 by using rin instead of raw_rin
          // only swap "outside" points that are still internal to overall curve
          // true by definition for inside points?), but just going for consistency across inner/outer modes
          if (longer > 0) { 
            double rdiff = longest - rin;
            rout = longest + rdiff;
            double rx= rout * cos(tin);
            double ry = rout * sin(tin);
            // or equivalently, (rx + rx - x, ry + ry - y) ??? (P' = I + (I-P))
            pVarTP.x += adjustedAmount * rx;
            pVarTP.y += adjustedAmount * ry;
          }
          else { // place on curve 
            pVarTP.x += adjustedAmount * x;
            pVarTP.y += adjustedAmount * y;
          }
          break;
        case 17: // HIDE
          pVarTP.x = pAffineTP.x;
          pVarTP.y = pAffineTP.y;
          pVarTP.doHide = true;
          break;
        case 18: // SCALE (inspired by Circus)
          pVarTP.x = pAffineTP.x * spread * spread_ratio;
          pVarTP.y = pAffineTP.y * spread;
          //          pVarTP.y = pAffineTP.y * (1/spread_ratio);
          break;
        case 19: // WHIRL (inspired by WhorlFunc)
          // a = pAffineTP.getPrecalcAtanYX() + inside / (pAmount - r);
          // a = pAffineTP.getPrecalcAtanYX() + outside / (pAmount - r);
          double a = tin + ((spread/10) / (r - rin));
          double sa = sin(a);
          double ca = cos(a);
          pVarTP.x += pAmount * rin * ca;
          pVarTP.y += pAmount * rin * sa;
          break;
        case 20:  // STRETCH1B
          // double rout = (rin * spread) + (r * (1 - spread));
          double rc = nearest;
          //if (longer > 0) { rc = longest; }
          //else if (shorter > 0) { rc = shortest; }
          // else { rc = r; }
          rout = rc + ((rin - rc) * spread);

          pVarTP.x += adjustedAmount * rout * cos(t);
          pVarTP.y += adjustedAmount * rout * sin(t);
          break; 
        case 21: // POW (inspired by Juliascope)
          rout = r + pow((sqr(x - pAffineTP.x) + sqr(y - pAffineTP.y)), spread) - 1.0;
          pVarTP.x += adjustedAmount * rout * cos(t);
          pVarTP.y += adjustedAmount * rout * sin(t);
          break;
        default:  // if mode specified has no definition, just leave on curve
          pVarTP.x += adjustedAmount * x;
          pVarTP.y += adjustedAmount * y;
          break;
      }
      //    pVarTP.z += adjustedAmount * pAffineTP.z;
      pVarTP.z += pAmount * pAffineTP.z;
    
  }
  
  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { 
                          curve_scale, 
                          mode_merging, 
                          inner_mode, outer_mode, outish_mode, 
                          inner_spread, outer_spread, outish_spread, 
                          inner_spread_ratio, outer_spread_ratio, outish_spread_ratio, 
                          spread_split,
                          cycles_param, cycle_offset, 
                          fill, curve_radius_mode, 
                          metacycles, metacycle_offset, metacycle_scale
    };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_CURVE_SCALE.equalsIgnoreCase(pName)) {
      curve_scale = pValue;
    } 
    else if (PARAM_MODE_MERGING.equalsIgnoreCase(pName)) {
      mode_merging = (pValue == 0 ? 0 : 1);
    }
    else if (PARAM_OUTER_MODE.equalsIgnoreCase(pName)) {
      outer_mode = (int)floor(pValue);
      // if (outer_mode > 15 || outer_mode < 0) { outer_mode = 0; }
    }
    else if (PARAM_INNER_MODE.equalsIgnoreCase(pName)) {
      inner_mode = (int)floor(pValue);
      // if (inner_mode > 15 || inner_mode < 0) { inner_mode = 0; }
    }
    else if (PARAM_OUTISH_MODE.equalsIgnoreCase(pName)) {
      outish_mode = (int)floor(pValue);
      // if (inner_mode > 15 || inner_mode < 0) { inner_mode = 0; }
    }
    else if (PARAM_OUTER_SPREAD.equalsIgnoreCase(pName))
      outer_spread = pValue;
    else if (PARAM_INNER_SPREAD.equalsIgnoreCase(pName))
      inner_spread = pValue;
    else if (PARAM_OUTISH_SPREAD.equalsIgnoreCase(pName))
      outish_spread = pValue;    
    else if (PARAM_OUTER_SPREAD_RATIO.equalsIgnoreCase(pName))
      outer_spread_ratio = pValue;
    else if (PARAM_INNER_SPREAD_RATIO.equalsIgnoreCase(pName))
      inner_spread_ratio = pValue;    
    else if (PARAM_OUTISH_SPREAD_RATIO.equalsIgnoreCase(pName))
      outish_spread_ratio = pValue;    
    else if (PARAM_SPREAD_SPLIT.equalsIgnoreCase(pName))
      spread_split = pValue;
    else if (PARAM_CYCLES.equalsIgnoreCase(pName))
      cycles_param = abs(pValue);
    else if (PARAM_CYCLE_OFFSET.equalsIgnoreCase(pName))
      cycle_offset = pValue;    
    else if (PARAM_FILL.equalsIgnoreCase(pName))
      fill = pValue;
    else if (PARAM_CURVE_RADIUS_MODE.equalsIgnoreCase(pName))
      this.curve_radius_mode = pValue;
    else if (PARAM_METACYCLES.equalsIgnoreCase(pName))
      metacycles = abs(pValue);
    else if (PARAM_METACYCLE_OFFSET.equalsIgnoreCase(pName))
      metacycle_offset = pValue;   
    else if (PARAM_METACYCLE_SCALE.equalsIgnoreCase(pName))
      metacycle_scale = pValue;     
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "abstract_polar_curve";
  }


}
