package org.jwildfire.create.tina.variation;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import static org.jwildfire.base.mathlib.MathLib.EPSILON;
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

public abstract class AbstractPolarCurveFunc extends VariationFunc {
  protected boolean DEBUG = false;
  protected boolean DRAW_DIAGNOSTICS = false;
  protected boolean DEBUG_MODES = false;

  public enum RenderMode { DEFAULT(0), 
                           ONCURVE(1),
                           STRETCH1(2),
                           STRETCH2(3),
                           UNCHANGED(4), 
                           HIDE(5),
                           SCALE(6),
                           MIRROR_SWAP(7), 
                           CURVE_SWAP(8),
                           CURVE_SWAP_RAW(9), 
                           CURVE_XY_OFFSET(10), 
                           CURVE_RADIAL_OFFSET(11), 
                           WHIRL(12), 
                           POW(13),
                           LOOPY(14), 
                           STRETCH3(15), 
                           STRETCH4(16),
                           STRETCH5(17),
                           STRETCH6(18),
                           STRETCH7(19), 
                           STRETCH8(20),
                           STRETCH9(21),
                           STRETCH10(22);
       
     private static final Map<Integer,RenderMode> lookup = new HashMap<Integer,RenderMode>();
     private static int maxInt = 0;
     
     static {
       for(RenderMode mode : values()) {
         int mint = mode.getIntegerMode();
         if (lookup.get(mode.getIntegerMode()) != null) { // for checking during dev
           throw new IllegalArgumentException("duplicate mode int: " + mode.name() + ", " + mint);
         }
         lookup.put(mode.getIntegerMode(), mode);
         maxInt = Math.max(maxInt, mode.getIntegerMode());
       }
     }
     private int modeInt;
     private RenderMode(int i) {
          this.modeInt= i;
     }
     private int getIntegerMode() { return modeInt; }
     
     // will return null if no RenderMode with that modeInt 
     public static RenderMode get(int i) {  return lookup.get(i);   }
     
     public  static int getMaxIntegerMode() { return maxInt; }
  }
  
  // "state" of point relative to curve
  //     currently only using INSIDE, OUTSIDE, OUTISH (and outish is only used with 
  //     CurveRadiusMode *_SAMPLING_BINs, and InsideOutsideRule MODIFIED_*s
  public enum PointState { INSIDE, OUTSIDE, INISH, OUTISH, ON }  
  public enum CurveRadiusMode {  AUTO, TRANSFORMED_INPUT, INTERPOLATED_SAMPLING_BIN, RAW_SAMPLING_BIN;
    // WARNING -- DO NOT CHANGE ORDER OF ENUMS, ADD NEW ONES TO END OF LIST
    public static CurveRadiusMode get(int oindex) {
      if (oindex < 0 || oindex >= values().length) { return null; }
      else  { return values()[oindex]; }
    }
    public int getIntegerMode() { return this.ordinal(); }
  }
  
  public enum PointRadiusMode { AUTO, MODIFIED, RAW }
  public enum InsideOutsideRule { AUTO, EVEN_ODD, EVEN_ODD_OUTISH, MODIFIED_EVEN_ODD, MODIFIED_EVEN_ODD_INISH_OUTISH, MOD2;
      public static InsideOutsideRule get(int oindex) {
      if (oindex < 0 || oindex >= values().length) { return null; }
      else  { return values()[oindex]; }
    }
    public int getIntegerMode() { return this.ordinal(); }
  }
  
  public enum MergeMode { AUTO,  // depends on mode?
                          NONE,  // no merging of modes
                          ALL, // inner* controls all (inner/outer/inish/outish) 
                          INNER_OUTERISH, // outer* controls both outer and outish
                          INNERISH_OUTER,   // inner* controls both inner and outish, 
                          INOUTER_OUTISH;   // inner* controls both inner and outer
    // WARNING -- DO NOT CHANGE ORDER OF ENUMS, ADD NEW ONES TO END OF LIST
    public static MergeMode get(int oindex) {
      if (oindex < 0 || oindex >= values().length) { return null; }
      else  { return values()[oindex]; }
    }
    public int getIntegerMode() { return this.ordinal(); }
  }
  
  protected class LinkedPolarCurvePoint {
    protected double radius;
    protected double angle;
    protected boolean inflection = false;
    protected LinkedPolarCurvePoint prev;
    protected LinkedPolarCurvePoint next;
    protected int bin;
    protected LinkedPolarCurvePoint(double r, double t) {
      radius = r;
      angle = t;
    }
    protected double interpolate(double ain) {
      double rout;
      if (next == null || prev == null) {
        rout = radius;
      }
      else {
        LinkedPolarCurvePoint pother = null;
        if (ain >= angle) {
          if (next.angle >= ain) {      // order is angle < ain < nexta
            pother = next;
          }
          else if (prev.angle >= ain) { // order is angle < ain < preva
            pother = prev;
          }
        }
        else {  // ain < angle
          if (next.angle <= ain) {      // order is nexta < ain < angle
            pother = next;
          }
          else if (prev.angle <= ain) { // order is preva < ain < angle
            pother = prev;
          }
        }
        if (pother == null) {
          rout = radius;
        }
        else {
          rout = this.radius + ((pother.radius - this.radius) * (ain - this.angle) / (pother.angle - this.angle));
        }
      }
      return rout;
    }
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
  protected static final String PARAM_LOCATION_CLASSIFIER = "location_classifier";
  protected static final String PARAM_CYCLE_ROTATION = "cycle_rotation";
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
                                               PARAM_CYCLES, PARAM_CYCLE_ROTATION, 
                                               PARAM_FILL, PARAM_CURVE_RADIUS_MODE, PARAM_LOCATION_CLASSIFIER, 
                                               PARAM_METACYCLES, PARAM_METACYCLE_OFFSET, PARAM_METACYCLE_SCALE };
//                                               PARAM_METACYCLE_ROTATION }; 

  protected double curve_scale = 1;

  protected MergeMode mode_merge = MergeMode.AUTO;
  protected RenderMode inner_mode = RenderMode.ONCURVE;
  protected RenderMode outer_mode = RenderMode.ONCURVE;
  protected RenderMode outish_mode = RenderMode.ONCURVE;
  protected CurveRadiusMode curve_rmode = CurveRadiusMode.AUTO;
  protected InsideOutsideRule location_mode = InsideOutsideRule.AUTO;

  protected MergeMode mode_merge_param = mode_merge;
  protected RenderMode inner_param = inner_mode;
  protected RenderMode outer_param = outer_mode;
  protected RenderMode outish_param = outish_mode;
  protected CurveRadiusMode curve_rmode_param = curve_rmode;
  protected InsideOutsideRule location_mode_param = location_mode;
  
  // point_rmode is currently hardwired, not available as a user param
  PointRadiusMode point_rmode = PointRadiusMode.MODIFIED;

  protected double inner_spread = 0.5; // deform based on original x/y
  protected double outer_spread = 0.5; // deform based on original x/y
  protected double outish_spread = 0.5;
  
  protected double inner_spread_ratio = 1; // how much inner_spread applies to x relative to y
  protected double outer_spread_ratio = 1; // how much outer_spread applies to x relative to y
  protected double outish_spread_ratio = 1;
  
  protected double spread_split = 1;
  protected double fill = 0;
  protected double cycles_param;
  
  // protected double point_radius_mode;
  // protected double in_out_mode;
  
  protected double cycles;  // 1 cycle = 2*PI
  protected double cycle_length = M_2PI; // 2(PI)
  protected double cycles_to_close = 0; // 0 indicates unknown, -1 indicates curve will never close
  protected double cycle_rotation = 0; // cycle rotation (in cycles) for incoming points (rotate the cycle)
  protected double metacycles = 1; // if cycles is calculated automatically to close the curve, metacycles is number of times to loop over closed curve
  protected double metacycle_offset = 0.1; // P(m) = P * (1 + mOffset) * ((mScale)^m))  // cumulative offset for metacycles
  protected double metacycle_scale = 1.1;
  // protected double metacycle_rotation = 0; // additional (cumulative) theta offset (in cycles) for each metacycle (rotate the metacycle)

  // vars for determining inner/outer via even-odd rule
  int default_sample_count = 36000;
  int binCount = 720;
  ArrayList<ArrayList<LinkedPolarCurvePoint>> theta_intersects = null;
  ArrayList<LinkedPolarCurvePoint> unisects;
  LinkedPolarCurvePoint unipolar;
  XYZPoint pCurve = new XYZPoint();

  @Override
  public void init(FlameTransformationContext pContext, Layer pLayer, XForm pXForm, double pAmount) {
    // System.out.println("calling init for VariationFunc: " + this);
    
    if (inner_param == RenderMode.DEFAULT) { inner_mode = RenderMode.ONCURVE; }
    else { inner_mode = inner_param; }
    if (outer_param == RenderMode.DEFAULT) { outer_mode = RenderMode.ONCURVE; }
    else { outer_mode = outer_param; }
    if (outish_param == RenderMode.DEFAULT) { outish_mode = RenderMode.ONCURVE; }
    else { outish_mode = outish_param; }
    
    // if MergeMode.AUTO, for now just set to combine OUTER and OUTISH
    if (mode_merge_param == MergeMode.AUTO) { mode_merge = MergeMode.INNER_OUTERISH; }
    else { mode_merge = mode_merge_param; }
    
    if (mode_merge == MergeMode.ALL) { // combine inner/outer/inish, use inner for all
      outer_mode = inner_mode;
      outish_mode = inner_mode;
    }
    else if (mode_merge == MergeMode.INNERISH_OUTER) {  // combine inner and outish, use inner for outish
      outish_mode = inner_mode;
    }
    else if (mode_merge == MergeMode.INNER_OUTERISH) { // combine outer and outish, use outer for outish
      outish_mode = outer_mode;
    }
    else if (mode_merge == MergeMode.INOUTER_OUTISH) { // combin inner and outer, ouse inner for outer
      outer_mode = inner_mode;
    }
    else if (mode_merge == MergeMode.NONE) {  } // do nothing, already separate by default
    
    boolean uses_sampling_modes = 
            (inner_mode == RenderMode.UNCHANGED || inner_mode == RenderMode.MIRROR_SWAP ||  inner_mode == RenderMode.SCALE ||
              outer_mode == RenderMode.UNCHANGED || outer_mode == RenderMode.MIRROR_SWAP || outer_mode == RenderMode.SCALE ||
              outish_mode == RenderMode.UNCHANGED || outish_mode == RenderMode.MIRROR_SWAP || outish_mode == RenderMode.SCALE);
    
    if (curve_rmode_param == CurveRadiusMode.AUTO) {
      // auto calling of curve radius mode?
      // possibly mixed?
      // if any OUTER/INNER/OUTISH are UNCHANGED or MIRROR_SWAP or SCAle, use linear interpolation
      // otherwise, use input point transformed to curve
      if (uses_sampling_modes) {
        curve_rmode = CurveRadiusMode.INTERPOLATED_SAMPLING_BIN;
      }
      else {
        curve_rmode = CurveRadiusMode.TRANSFORMED_INPUT;
      }
    }
    else {
      curve_rmode = curve_rmode_param;
    }
    
    if (location_mode_param == InsideOutsideRule.AUTO) {
      location_mode = InsideOutsideRule.EVEN_ODD_OUTISH;
    }
    else { location_mode = location_mode_param; }
    
    if (DEBUG_MODES) {
      System.out.println("mode_merge: " + this.mode_merge.name());
      System.out.println("inner_mode: " + inner_mode.name());
      System.out.println("outer_mode: " + outer_mode.name());
      System.out.println("outish_mode: " + outish_mode.name());
      System.out.println("curve_radius_mode: " + this.curve_rmode.name());
      System.out.println("point_radius_mode: " + this.point_rmode.name());
      System.out.println("location_mode: " + this.location_mode.name());
      System.out.println("============================================");
    }

    unipolar = new LinkedPolarCurvePoint(0, 0);
    unisects = new ArrayList<LinkedPolarCurvePoint>(1);
    unisects.add(0, unipolar);
    recalcCycles();

    // System.out.println("curve_rmode: " + curve_rmode.name());
    if (curve_rmode == CurveRadiusMode.RAW_SAMPLING_BIN || curve_rmode == CurveRadiusMode.INTERPOLATED_SAMPLING_BIN) {
      recalcCurveIntersects();
    }
  }
  
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
    // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    //    range of atan2 is from [0 --> PI] for  positive y as x:[+->0->-]  (at x=0, atan2 = PI/2)
    //                 and  from [0 --> -PI] for negative y as x:[+->0->-]  (at y=0, atan2 = -PI/2)
    double tin = atan2(pAffineTP.y, pAffineTP.x);  
    
    // then stretch theta over full number of cycles
    // so range of theta is [0 --> cycles*PI] and [0 --> -cycles*PI], or overall [-cycles*PI --> cycles*PI]
    //      with range delta of cycles*2PI
    double theta = (cycles * tin);  

    // and add cycle_rotation
    theta += cycle_rotation * M_2PI;
    
    // use scratch XYZPoint pCurve to calc points on curve
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
      // double metacycle_count = floor((theta + (cycles * M_PI)) / (cycles_to_close * M_2PI));
      // metacycle_count = floor((theta + (cycles * M_PI) - (cycle_rotation * M_2PI)) / (cycles_to_close * M_2PI));
        double metacycle_progress = (theta + (cycles * M_PI) - (cycle_rotation * M_2PI)) / (cycles_to_close * M_2PI);
        // need to check if metacycles calcs are very slightly too low or too high
        //    (I think due to double rounding errors?)
        if (metacycle_progress < 0) {
          System.out.println("in transform(), metacycle < 0: " + metacycle_progress + " theta: " + theta + " theta/PI: " + theta/M_PI);
          // see if adding epsilon will raise
          metacycle_progress += EPSILON;
        }
        else if (metacycle_progress >= metacycles) {
          // see if subtracting epsilon will lower
          System.out.println("metacycle >= metacycles: " + metacycle_progress + " theta: " + theta + " theta/PI: " + theta/M_PI);
          metacycle_progress -= EPSILON;
        }
        double metacycle_count = (int)floor(metacycle_progress);

        if (metacycle_count > 0) { 
          double metacycle_delta = (metacycle_count * metacycle_offset) + (pow(metacycle_scale, metacycle_count)-1);
          pResult.x = pResult.x * (1 + metacycle_delta);
          pResult.y = pResult.y * (1 + metacycle_delta);
          // z unchanged?
        }
    }
  }  


  int recalcCount = 0; 
  boolean DEBUG_INTERSECTS = false;
   public void recalcCurveIntersects() {
    // System.out.println("recalcing curves");
    theta_intersects = new ArrayList<ArrayList<LinkedPolarCurvePoint>>(binCount);
    for (int i=0; i<binCount; i++) { 
      theta_intersects.add(new ArrayList<LinkedPolarCurvePoint>());
    }
    LinkedPolarCurvePoint prev_point = null;
    // PolarPoint first_point = null;
    LinkedPolarCurvePoint last_point = null;
    // PolarPoint next_point;
    ArrayList<LinkedPolarCurvePoint> tsects;
    ArrayList<LinkedPolarCurvePoint> prev_tsects = null;
    int firstbin = -1;
    int lastbin = -1;
    int sampleCount = default_sample_count;
    if (cycles_param == 0 && cycles_to_close > 0 && metacycles > 1) {
      sampleCount = (int)(sampleCount * metacycles);
    }
    // int prev_metacycle = (int)(this.metacycles - 1);
    int prev_metacycle = -1000;
    ArrayList<LinkedPolarCurvePoint> metacycle_first_points = new ArrayList<LinkedPolarCurvePoint>((int)ceil(metacycles));  // for each metacycle, keep track of first point (for looping)
    ArrayList<LinkedPolarCurvePoint> metacycle_last_points = new ArrayList<LinkedPolarCurvePoint>((int)ceil(metacycles)); // for each metacycle, keep track of last point (for looping)
    for (int m=0; m<(int)ceil(metacycles); m++) {
      metacycle_first_points.add(null);
      metacycle_last_points.add(null);
    }
    // System.out.println("metacycles: " + metacycles + ", meta_points.size(): " + metacycle_first_points.size());

    for (int i=0; i<sampleCount; i++) {
      // want theta to span expected thetas used as input to calcCurvePoint(), so range needs to be       
      //    same as theta calculated in transform() (that is used as input to calcCurvePoint())
      // double theta = ((((double)i/(double)sampleCount) - 0.5) * cycles * M_2PI) + (cycle_rotation * M_2PI);
      double theta = (((double)i/(double)sampleCount) * cycles * M_2PI)  - (cycles * M_PI) + (cycle_rotation * M_2PI);
      // or equivalently, 
      // double theta = ((((double)i/(double)sampleCount) - 0.5) * cycles * M_2PI) + (cycle_rotation * M_2PI);
      double cycles_for_calc;
      if (cycles_to_close > 0)  { cycles_for_calc = cycles_to_close; }
      else { cycles_for_calc = cycles; }

      pCurve.clear();
      calcCurvePoint(null, theta, pCurve);
      
      double r = sqrt(pCurve.x * pCurve.x + pCurve.y * pCurve.y);
      double angle = atan2(pCurve.y, pCurve.x);
      int anglebin =  (int)Math.floor(((angle + M_PI)/M_2PI) * binCount);
      // catching any possible cases where angle actually reaches max atan2
      //   (actually, maybe should loop around instead to bin 0?
      if (anglebin == binCount) { anglebin--; } 
      tsects = theta_intersects.get(anglebin);

      // still rotating through same bin, merge results
      if (prev_tsects == tsects) {  
        // try ignoring for now -- should try averaging later?
      }
      else {
        // tsects.add(r);
        LinkedPolarCurvePoint point = new LinkedPolarCurvePoint(r, angle);
        point.bin = anglebin;
        
        double metacycle_progress = (theta + (cycles * M_PI) - (cycle_rotation * M_2PI)) / (cycles_for_calc * M_2PI);
        // need to check if metacycles calcs are very slightly too low or too high
        //    (I think due to double rounding errors?)
        if (metacycle_progress < 0) {
          if (DEBUG_INTERSECTS) { System.out.println("metacycle < 0: " + metacycle_progress + ", samplecount: " + i + ", theta: " + theta + " theta/PI: " + theta/M_PI); }
          // see if adding epsilon will raise
          metacycle_progress += EPSILON;
        }
        else if (metacycle_progress >= metacycles) {
          if (DEBUG_INTERSECTS) { System.out.println("metacycle >= metacycles: " + metacycle_progress + ", samplecount: " + i + ", theta: " + theta + " theta/PI: " + theta/M_PI); }
          // see if subtracting epsilon will lower
          metacycle_progress -= EPSILON;
        }
        int metacycle_count = (int)floor(metacycle_progress);
        
        // int metacycle = (int)(ceil((theta + (cycles * M_PI) + (cycle_rotation * M_2PI))/ (cycles_for_calc * M_2PI)) - 1);
        if (metacycle_count == prev_metacycle) { // still in same metacycle
          point.prev = prev_point;
          if (prev_point != null) { prev_point.next = point; }
        }
        else {
          if (DEBUG_INTERSECTS) { System.out.println("new metacycle: " + metacycle_count + ", bin: " + anglebin + ", radius: " + r + ", angle: " + angle); }
          // initialize new metacycle, close previous metacycle
          if (metacycle_count < 0 || metacycle_count > metacycles) {
            if (DEBUG_INTERSECTS) { System.out.println("error: " + ((theta + (cycles * M_PI) - (cycle_rotation * M_2PI)) / (cycles_to_close * M_2PI))); }
          }
          metacycle_first_points.set(metacycle_count, point);
          if (prev_point != null && prev_metacycle >= 0) {
            metacycle_last_points.set(prev_metacycle, prev_point);
          }
        }
        tsects.add(point);  // autoboxing float r to Double object
        if (i == 0) {
          firstbin = anglebin;
          // first_point = point;
        }
        prev_point = point;
        prev_metacycle = metacycle_count;
      }

      if (i == (sampleCount-1)) { 
        lastbin = anglebin;
        last_point = prev_point;
      }
      prev_tsects = tsects;
    }
    
    metacycle_last_points.set(prev_metacycle, prev_point);
    for (int m=0; m<metacycle_last_points.size(); m++) {
      LinkedPolarCurvePoint fp = metacycle_first_points.get(m);
      LinkedPolarCurvePoint lp = metacycle_last_points.get(m);
      if (DEBUG_INTERSECTS && (fp.bin == 0 || lp.bin == 0)) {
        System.out.println("metacycle boundary falls on zero bin: ");
        printBin(0);
        printBin(binCount-1);
      }
      LinkedPolarCurvePoint prev = new LinkedPolarCurvePoint(lp.radius, lp.angle);
      LinkedPolarCurvePoint next = new LinkedPolarCurvePoint(fp.radius, fp.angle);
      prev.bin = lp.bin;
      next.bin = fp.bin;
      fp.prev = prev;
      lp.next = next;
    }

    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN PRE:"); printBin(0); System.out.println("MAXBIN PRE:");     printBin(binCount-1); }
    //    if (last_point != null)  {  // connect to first point of current metacycle?? }

    // bin angles rotate through single circle, [-PI ==> PI] (actual point angles in min/max bin  will be near bin min/max val]
    //       (regardless of theta, metacycles, etc.)
    //    so need to fix prev/next point angles in both min and max bin, 
    //           since one of prev/next will be 2PI off when it crosses line
    //    which of prev/next needs to be fixed will depend on curveCalc
    //        (and _should_ be determinable by bin number comparison, 
    //             but appears to be bug where at least for bin 0, prev/next can point back to point in same bin?)
    //        but don't need to know this, can test simply by closeness to point.angle, since if cross line 
    //             then prev/next angle should be ~2PI away. 
    //    using PI to check instead of 2PI, could use somethiing much closer to 2PI but PI allows should catch everything
    //        without worrying about exact size of bin (except for _really_ extreme circumstances)
    //    UPDATE: fixed bug mentioned above, possibly switch to only changing prev/next that is crosses bin count loop boundary?
    //        for now, still checking both prev/next and changing based on andle delta to point angle
    ArrayList<LinkedPolarCurvePoint> minBin = theta_intersects.get(0);
    ArrayList<LinkedPolarCurvePoint> maxBin = theta_intersects.get(binCount-1);
    ArrayList<LinkedPolarCurvePoint> minmax[] = new ArrayList[]{minBin, maxBin};

    for (ArrayList<LinkedPolarCurvePoint> bin : minmax) {
      for (LinkedPolarCurvePoint point : bin) {
        // prev and null should already be cloned in metacycle step above 
        //    (since metacycle boundary will always fall on bin loop boundary? -- hmm, not sure about this)
        // but cloning again just in case, don't want to mess up others
        //      (after calcIntersects(), prev/next are only used for interpolation, 
        //      and past this point in calcIntersects() they are not recursively crawled, so don't need deep clone, 
        //      just radius and corrected angle (also setting bin just for consistency))
        if (point.prev != null) {
          if (abs(point.angle-point.prev.angle) > M_PI) {
            double nangle; 
            if (point.prev.angle < point.angle) { nangle = point.prev.angle + M_2PI; }
            else { nangle = point.prev.angle - M_2PI; }
            LinkedPolarCurvePoint new_prev = new LinkedPolarCurvePoint(point.prev.radius, nangle);
            new_prev.bin = point.prev.bin;
            point.prev = new_prev;
          }
        }
        if (point.next != null) {
          if (abs(point.angle-point.next.angle) > M_PI) {
            double nangle;
            if (point.next.angle < point.angle) { nangle = point.next.angle + M_2PI; }
            else { nangle = point.next.angle - M_2PI; }
            LinkedPolarCurvePoint new_next = new LinkedPolarCurvePoint(point.next.radius, nangle);
            new_next.bin = point.next.bin;
            point.next = new_next;
          }
        }
      }
    }
    
    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN MID:"); printBin(0); System.out.println("MAXBIN MID:");     printBin(binCount-1); }
    
    // cleanup of PolarPoint prev/next pointers -- 
    //     cyclic, but also need to factor in metacylces
    //     if first of metacycle, keep (meta_first)
    //     if last of metacycle (meta_last), set metal_last.next = meta_first, meta_first.prev = meta_last
/*    int metacycle = metacycles / cycles;
    if (cycles_param == 0 && cycles_to_close > 0 && metacycles != 1) {
      // metacycle_count does not include first cycles, so goes from 0 to metacycles - 1 (since metacycles includes first cycle)
      double metacycle_count = ceil((theta + (cycles * M_PI)) / (cycles_to_close * 2 * M_PI)) - 1;
    }  
    */
    
    // special-casing of first and last anglebin if they are the same bin:
    ///   want to simulate rotating through same bin to merge "duplicate" intersections
    //    if first and last bin are same, would have merged results, so remove last one 
   /*
    if (firstbin > 0 && lastbin > 0 && firstbin == lastbin) {
      tsects = theta_intersects.get(firstbin);
      // remove last result (if start doing averaging, should remove but add to average for first result)
      tsects.remove(tsects.size()-1);
    }
    */
    
    // WARNING: still need to factor in metacycle
    for (int k=0; k<theta_intersects.size(); k++) {
      ArrayList<LinkedPolarCurvePoint> bin = theta_intersects.get(k);
            
      for (LinkedPolarCurvePoint p : bin) {
        if (p.prev == null || p.next == null) { 
          p.inflection = false;
        }
        else {
          double aprev = p.prev.angle;
          double anext = p.next.angle;
          if (p.angle <= aprev && p.angle <= anext) {
            p.inflection = true;
            if (DEBUG_INTERSECTS) { System.out.println("- inflection at bin: " + k + ", " + p.angle); }
          }
          else if (p.angle >= aprev && p.angle >= anext) {
            p.inflection= true;
            if (DEBUG_INTERSECTS) { System.out.println("+ inflection at bin: " + k + ", " + p.angle); }
          }
          else { p.inflection = false; }
        }
      }
    }

    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN POST:"); printBin(0); System.out.println("MAXBIN POST:");     printBin(binCount-1); }
  }
  
public void printBin(int index) {
  ArrayList<LinkedPolarCurvePoint> bin = theta_intersects.get(index);
  System.out.println("bin: " + index + ", points: " + bin.size());
  for (int i=0; i< bin.size(); i++) {
    LinkedPolarCurvePoint p = bin.get(i);
    System.out.println("    point: " + i + ", inflect: " + p.inflection + ", a = " + p.angle + ", r = " + p.radius + ", prev = " + (p.prev!=null) + ", next = " + (p.next!=null));
    if (p.prev != null) { 
      System.out.println("         prev: " + p.prev.bin + ", i: " + p.prev.inflection + ", a: " + p.prev.angle + ", r: " + p.prev.radius);
    }      
    if (p.next != null) { 
      System.out.println("         next: " + p.next.bin + ", i: " + p.next.inflection + ", a: " + p.next.angle + ", r: " + p.next.radius);
    }      
  }
}

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
    if (fill != 0) {
      r = r + (fill * (pContext.random() - 0.5));
      x = r * cos(t);
      y = r * sin(t);
    }
    ArrayList<LinkedPolarCurvePoint> tsects;

    int anglebin =  (int)Math.floor(((tin + M_PI)/M_2PI) * binCount);
 
    if (curve_rmode == CurveRadiusMode.RAW_SAMPLING_BIN || curve_rmode == CurveRadiusMode.INTERPOLATED_SAMPLING_BIN) {
      if (anglebin == binCount) {  // catching any possible cases where tin actually reaches max atan2
        anglebin--; 
      } 
      tsects = theta_intersects.get(anglebin);
    }
    else if (curve_rmode == CurveRadiusMode.TRANSFORMED_INPUT)  {
      tsects = unisects;
      unipolar.radius = r;
      unipolar.angle = t;
    }
    else {  // default to CurveRadiusMode.THETA ??
      tsects = unisects;
      unipolar.radius = r;
      unipolar.angle = t;
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

    for (LinkedPolarCurvePoint curve : tsects) {
      if (curve.inflection) { continue; } // if inflection, don't count as an intersection
      double rcurve; 
      if (curve_rmode == CurveRadiusMode.INTERPOLATED_SAMPLING_BIN) {
        rcurve = curve.interpolate(tin);
      }
      else { rcurve = curve.radius; }
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
    
    if (location_mode == InsideOutsideRule.EVEN_ODD) {
      // use standard Even-Odd rule (well, more standard than the modified one above...):
      //    cast ray from origin through incoming point to infinity
      //    count how many times curve intersects ray further out than incoming point (longer)
      //    if number is odd then point is inside, if number is even then point is outside
      if (longer % 2 == 0) { pstate = PointState.OUTSIDE; } // point is outside
      else { pstate = PointState.INSIDE; } // point is inside
    }
    else if (location_mode == InsideOutsideRule.EVEN_ODD_OUTISH) {
      // use standard Even-Odd rule as above (well, more standard than the modified ones below...), 
      //    but distinguish between OUTER an OUTISH:
      //    OUTISH if standard Even-Odd calls as OUTER, but still have at least one curve intersect on ray 
      //      further out from origin than point is
      if (longer % 2 == 0) { 
        if (longer > 0) {
          pstate = PointState.OUTISH;
        }
        else { 
          pstate = PointState.OUTSIDE; } // point is outside
      }
      else { pstate = PointState.INSIDE; } // point is inside
    }
    else if (location_mode == InsideOutsideRule.MOD2) {
      if (longer == 0) { // definitely outside
        pstate = PointState.OUTSIDE;
      }
      else if (shorter == 0) { // definitely inside
        pstate = PointState.INSIDE;
      }
      else if (longer % 2 == 0) { 
        pstate = PointState.OUTISH;
      }
    }

    else if (location_mode == InsideOutsideRule.MODIFIED_EVEN_ODD_INISH_OUTISH) {
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
        
    else if (location_mode == InsideOutsideRule.MODIFIED_EVEN_ODD) {
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

      // 3. Prep input for feeding to switch statement:
      //       using designated inner/outer/outish/(on?) for
      //       *_mode, *_spread, *_spread_ratio, etc.
      double spread;
      double spread_ratio;
      RenderMode mode = RenderMode.DEFAULT;
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
        mode = RenderMode.DEFAULT;
        spread = 0;
        spread_ratio = 1;
      }
      else {
        // every point must be classified as being within one of the enum'd areas, so should never get here
        mode = RenderMode.DEFAULT;
        spread = 0;
        spread_ratio = 0;
      }
      
      if (mode == RenderMode.DEFAULT) {
        mode = RenderMode.ONCURVE;
      }
      
      // 4. pVarTP transformations, based on mode (for given INSIDE/OUTSIDE/OUTISH classification given above)
      switch (mode) {
        case ONCURVE: // ONCURVE -- no spread, place on curve
          pVarTP.x += adjustedAmount * x;
          pVarTP.y += adjustedAmount * y;
          break;
        case STRETCH1:  // STRETCH1
          // double rout = (rin * spread) + (r * (1 - spread));
          rout = r + ((rin - r) * spread);
          pVarTP.x += adjustedAmount * rout * cos(t);
          pVarTP.y += adjustedAmount * rout * sin(t);
          break;
        case STRETCH2: // STRETCH2
          rinx = (rin * spread * spread_ratio) - (spread * spread_ratio) + 1;
          riny = (rin * spread) - spread + 1;
          pVarTP.x += adjustedAmount * rinx * x;
          pVarTP.y += adjustedAmount * riny * y;
          break;
        case UNCHANGED:   // UNCHANGED -- leave in place
          // point is inside curve, leave in place
          pVarTP.x += adjustedAmount * pAffineTP.x;
          pVarTP.y += adjustedAmount * pAffineTP.y;
          break;
        case MIRROR_SWAP:  // MIRROR_SWAP (swap around origin [0,0], inspired by RosoniFunc)
          pVarTP.x += adjustedAmount * pAffineTP.x * -1;
          pVarTP.y += adjustedAmount * pAffineTP.y * -1;
          break;
        case CURVE_XY_OFFSET: // CURVE_XY_OFFSET -- (UNCHANGED + ONCURVE) combo
          pVarTP.x += adjustedAmount * (pAffineTP.x + x);
          pVarTP.y += adjustedAmount * (pAffineTP.y + y);
          break;
        case CURVE_RADIAL_OFFSET: // CURVE_RADIUS_OFFSET -- offset with curve intersect with longest radius at input angle, with: P' = C + P
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
        case CURVE_SWAP_RAW:  // CURVE_SWAP_RAW -- swap around curve intersect with longest radius at input angle,  with P' = C + (C-P)
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
        case CURVE_SWAP:  // CURVE_SWAP swap around intersect with longest radius, with P' = C + (C-P)
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
        case HIDE: // HIDE
          pVarTP.x = pAffineTP.x;
          pVarTP.y = pAffineTP.y;
          pVarTP.doHide = true;
          break;
        case SCALE: // SCALE (inspired by Circus)
          pVarTP.x = pAffineTP.x * spread * spread_ratio;
          pVarTP.y = pAffineTP.y * spread;
          //          pVarTP.y = pAffineTP.y * (1/spread_ratio);
          break;
        case WHIRL: // WHIRL (inspired by WhorlFunc)
          // a = pAffineTP.getPrecalcAtanYX() + inside / (pAmount - r);
          // a = pAffineTP.getPrecalcAtanYX() + outside / (pAmount - r);
          double a = tin + ((spread/10) / (r - rin));
          double sa = sin(a);
          double ca = cos(a);
          pVarTP.x += adjustedAmount * rin * ca;
          pVarTP.y += adjustedAmount * rin * sa;
          break;
        case STRETCH10:  // STRETCH1B (only differs from STRETCH1 when using binning)
          // double rout = (rin * spread) + (r * (1 - spread));
          double rc = nearest;
          //if (longer > 0) { rc = longest; }
          //else if (shorter > 0) { rc = shortest; }
          // else { rc = r; }
          rout = rc + ((rin - rc) * spread);

          pVarTP.x += adjustedAmount * rout * cos(t);
          pVarTP.y += adjustedAmount * rout * sin(t);
          break; 
        case POW: // POW (inspired by Juliascope)
          rout = r + pow((sqr(x - pAffineTP.x) + sqr(y - pAffineTP.y)), spread) - 1.0;
          pVarTP.x += adjustedAmount * rout * cos(t);
          pVarTP.y += adjustedAmount * rout * sin(t);
          break;
        case LOOPY:  // LOOPY (inspired by loonie
          rout = pAmount * sqrt((curvePoint.getPrecalcSumsq() / pAffineTP.getPrecalcSumsq()) - 1.0);
          pVarTP.x += adjustedAmount * rout * pAffineTP.x;
          pVarTP.y += adjustedAmount * rout * pAffineTP.y;
          break;
        case STRETCH3: // STRETCH3
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += adjustedAmount * (x - (spread * spread_ratio * (x-xin)));
          pVarTP.y += adjustedAmount * (y - (spread * (y-yin)));
          break;
        case STRETCH4: // STRETCH4
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += adjustedAmount * (x - (spread * spread_ratio * xin));
          pVarTP.y += adjustedAmount * (y - (spread * yin));
          break;
        case STRETCH5: // STRETCH5
          rinx = (0.5 * rin) + (spread * spread_ratio);
          riny = (0.5 * rin) + spread;
          pVarTP.x += adjustedAmount * rinx * x;
          pVarTP.y += adjustedAmount * riny * y;
          break;
        case STRETCH6: // STRETCH6 -- same as mode 3, but without the sign modifications
          pVarTP.x += adjustedAmount * (x + (spread * spread_ratio * pAffineTP.x));
          pVarTP.y += adjustedAmount * (y + (spread * pAffineTP.y));
          break;
        case STRETCH7: // STRETCH7 -- similar to 3, different sign fiddling
          xin = Math.abs(pAffineTP.x);
          yin = Math.abs(pAffineTP.y);
          if (x<0) { xin = xin * -1; }
          if (y<0) { yin = yin * -1; }
          pVarTP.x += adjustedAmount * (x + (spread * spread_ratio * xin));
          pVarTP.y += adjustedAmount * (y + (spread * yin));
          break;
        case STRETCH8: // STRETCH8 -- same as mode 6, but without the sign modifications
          pVarTP.x += adjustedAmount * (x - (spread * spread_ratio * pAffineTP.x));
          pVarTP.y += adjustedAmount * (y - (spread * pAffineTP.y));
          break;
        case STRETCH9: // STRETCH9
          pVarTP.x += adjustedAmount * rin * cos(t) * (spread * spread_ratio);
          pVarTP.y += adjustedAmount * rin * sin(t) * spread;
          break;
        default:  // if mode specified has no definition, just leave on curve
          pVarTP.x += adjustedAmount * x;
          pVarTP.y += adjustedAmount * y;
          break;
      }
      //    pVarTP.z += adjustedAmount * pAffineTP.z;
      pVarTP.z += pAmount * pAffineTP.z;

      if (DRAW_DIAGNOSTICS) {
        drawDiagnostics(pContext, pVarTP);
      }
  }
  
  protected void drawDiagnostics(FlameTransformationContext pContext, XYZPoint pVarTP) {
        double diagnostic = pContext.random() * 200;
      // draw diagnostic unit circles
      if (diagnostic == 0) {
        // ignore zero
      }
      if (diagnostic <= 4) { // diagnostic = (0-4]
        double radius = ceil(diagnostic)/2; // radius = 0.5, 1, 1.5, 2
        double angle = diagnostic * 2 * M_PI; // in radians, ensures coverage of unit circles
        pVarTP.x = radius * cos(angle);
        pVarTP.y = radius * sin(angle);
      }
      // draw diagnostic unit squares
      else if (diagnostic <= 8) { // diagnostic = (4-8]
        double unit = (ceil(diagnostic) - 4)/2; // unit = 0.5, 1, 1.5, 2
        int side = (int) ceil(5 * (ceil(diagnostic) - diagnostic)); // side = 1, 2, 3, 4
        double varpos = (pContext.random() * unit * 2) - unit;
        double sx = 0, sy = 0;
        if (side == 1) {
          sx = unit;
          sy = varpos;
        }
        else if (side == 2) {
          sx = varpos;
          sy = unit;
        }
        else if (side == 3) {
          sx = -1 * unit;
          sy = varpos;
        }
        else if (side == 4) {
          sx = varpos;
          sy = -1 * unit;
        }
        pVarTP.x = sx;
        pVarTP.y = sy;
      }
      else if (diagnostic <= 9) {
        // x = 0 gridline
        pVarTP.x = 0;
        pVarTP.y = (4 * (diagnostic - 8)) - 2; // line where y = [-2, 2]
      }
      else if (diagnostic <= 10) {
        // y = 0 gridline
        pVarTP.x = (4 * (diagnostic - 9)) - 2; // line where y = [-2, 2]
        pVarTP.y = 0;
      }
  }
  
  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] { 
                          curve_scale, 
                          mode_merge_param.getIntegerMode(),  
                          inner_param.getIntegerMode(), outer_param.getIntegerMode(), outish_param.getIntegerMode(), 
                          inner_spread, outer_spread, outish_spread, 
                          inner_spread_ratio, outer_spread_ratio, outish_spread_ratio, 
                          spread_split,
                          cycles_param, cycle_rotation, 
                          fill, curve_rmode_param.getIntegerMode(), location_mode_param.getIntegerMode(), 
                          metacycles, metacycle_offset, metacycle_scale
    };
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_CURVE_SCALE.equalsIgnoreCase(pName)) {
      curve_scale = pValue;
    } 
    else if (PARAM_MODE_MERGING.equalsIgnoreCase(pName)) {
      this.mode_merge_param = MergeMode.get((int)floor(pValue));
      if (mode_merge_param == null) { mode_merge_param = MergeMode.AUTO; }
    }
    else if (PARAM_OUTER_MODE.equalsIgnoreCase(pName)) {
      outer_param = RenderMode.get((int)floor(pValue));
      if (outer_param == null) { outer_param = RenderMode.DEFAULT; }
    }
    else if (PARAM_INNER_MODE.equalsIgnoreCase(pName)) {
      inner_param = RenderMode.get((int)floor(pValue));
      if (inner_param == null) { inner_param = RenderMode.DEFAULT; }
    }
    else if (PARAM_OUTISH_MODE.equalsIgnoreCase(pName)) {
      outish_param = RenderMode.get((int)floor(pValue));
      if (outish_param == null)  { outish_param = RenderMode.DEFAULT; }
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
    else if (PARAM_CYCLE_ROTATION.equalsIgnoreCase(pName))
      cycle_rotation = pValue;    
    else if (PARAM_FILL.equalsIgnoreCase(pName))
      fill = pValue;
    else if (PARAM_CURVE_RADIUS_MODE.equalsIgnoreCase(pName)) {
      this.curve_rmode_param = CurveRadiusMode.get((int)floor(pValue));
      if (curve_rmode_param == null) { curve_rmode_param = CurveRadiusMode.AUTO; }
    }
    else if (PARAM_LOCATION_CLASSIFIER.equalsIgnoreCase(pName)) {
      this.location_mode_param = InsideOutsideRule.get((int)floor(pValue));
      if (location_mode_param == null) { location_mode_param = InsideOutsideRule.AUTO; }
    }
    else if (PARAM_METACYCLES.equalsIgnoreCase(pName)) {
      metacycles = abs(pValue);
      if (abs(metacycles) < 0.01) { metacycles = 0.01; }
    }
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
