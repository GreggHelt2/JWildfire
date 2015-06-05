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
  boolean DEBUG_INTERSECTS = false;

  public enum RenderMode { DEFAULT(0), 
                           ONCURVE(1),
                           STRETCH1(2),

                           STRETCH2(3),
                           UNCHANGED(4), 
                           HIDE(5),
                           SCALE(6),
                           MIRROR_SWAP(7), 
                           CURVE_SWAP(8),
                           ROTATE(9), 
                           CURVE_XY_OFFSET(10), 
                           CURVE_RADIAL_OFFSET(11), 
                           WHIRL(12), 
                           POW(13),
                           LOOPY(14), 
                           STRETCH3(15), 
                           STRETCH4(16),
                           STRETCH5(17),
                           STRETCH7(19), 
                           STRETCH8(20),
                           STRETCH9(21), 
                           // TEST(30), 
                           TEST2(31),
                           XY_SWAP(18), 
                           ROTATE_RADIAL(30), 
                           END(50);
       
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

  public enum CurveProximityMode { AUTO, TRANSFORMED_RT, TRANSFORMED_R, LONGEST, NEAREST_LONGER, NEAREST, NEAREST_SHORTER, SHORTEST, NEAREST_NONZERO, SHORTEST_NONZERO;
    public static CurveProximityMode get(int oindex) {
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
  protected static final String PARAM_OUTISH_MODE = "outish_mode";
  // tried an "inish" mode too, but was problematic...
  // protected static final String PARAM_INISH_MODE = "inish_mode";  
  protected static final String PARAM_INNER_RSPREAD = "inner_rspread";
  protected static final String PARAM_OUTER_RSPREAD = "outer_rspread";
  protected static final String PARAM_OUTISH_RSPREAD = "outish_rspread";
  protected static final String PARAM_INNER_XSPREAD = "inner_xspread";
  protected static final String PARAM_OUTER_XSPREAD = "outer_xspread";
  protected static final String PARAM_OUTISH_XSPREAD = "outish_xspread";
  protected static final String PARAM_INNER_YSPREAD = "inner_yspread";
  protected static final String PARAM_OUTER_YSPREAD = "outer_yspread";
  protected static final String PARAM_OUTISH_YSPREAD = "outish_yspread";
  protected static final String PARAM_INNER_RBLUR = "inner_rblur";
  protected static final String PARAM_OUTER_RBLUR = "outer_rblur";
  protected static final String PARAM_OUTISH_RBLUR = "outish_rblur";
  protected static final String PARAM_INNER_ABLUR = "inner_ablur";
  protected static final String PARAM_OUTER_ABLUR = "outer_ablur";
  protected static final String PARAM_OUTISH_ABLUR = "outish_ablur";
  
  protected static final String PARAM_SPREAD_SPLIT = "spread_split";
  protected static final String PARAM_FILL = "fill";
  protected static final String PARAM_CURVE_THICKNESS = "curve_thickness";
  protected static final String PARAM_CYCLES = "cycles";
  protected static final String PARAM_CURVE_RADIUS_MODE = "curve_radius_mode";
  protected static final String PARAM_LOCATION_CLASSIFIER = "location_classifier";
  protected static final String PARAM_CYCLE_ROTATION = "cycle_rotation";
  // metacycles (and metacycle_expansion) are only in effect when cycles = 0 (automatic cycle calculations in effect), 
  //   and cycles to close can be determined (cycles_to_close > 0)
  protected static final String PARAM_METACYCLES = "metacycles";
  protected static final String PARAM_METACYCLE_OFFSET = "metacycle_offset";
  protected static final String PARAM_METACYCLE_SCALE = "metacycle_scale";
  protected static final String PARAM_METACYCLE_ROTATION = "metacycle_rotation";
  protected static final String PARAM_PROXIMITY_MODE = "proximity_mode";
  protected static final String PARAM_ANGLE_BINS = "angle_bins";
  
  protected static final String[] paramNames = { 
                                               PARAM_CURVE_SCALE, 
                                               PARAM_MODE_MERGING, 
                                               PARAM_INNER_MODE, PARAM_OUTER_MODE, PARAM_OUTISH_MODE, 
                                               PARAM_INNER_RSPREAD, PARAM_OUTER_RSPREAD, PARAM_OUTISH_RSPREAD, 
                                               PARAM_INNER_XSPREAD, PARAM_OUTER_XSPREAD, PARAM_OUTISH_XSPREAD, 
                                               PARAM_INNER_YSPREAD, PARAM_OUTER_YSPREAD, PARAM_OUTISH_YSPREAD, 
                                               PARAM_INNER_RBLUR, PARAM_OUTER_RBLUR, PARAM_OUTISH_RBLUR, 
                                               PARAM_INNER_ABLUR, PARAM_OUTER_ABLUR, PARAM_OUTISH_ABLUR, 
                                               PARAM_SPREAD_SPLIT,
                                               PARAM_CYCLES, PARAM_CYCLE_ROTATION, 
                                               PARAM_CURVE_THICKNESS, PARAM_FILL, 
                                               PARAM_PROXIMITY_MODE, PARAM_CURVE_RADIUS_MODE, PARAM_LOCATION_CLASSIFIER, 
                                               PARAM_ANGLE_BINS, 
                                               PARAM_METACYCLES, PARAM_METACYCLE_OFFSET, PARAM_METACYCLE_SCALE, PARAM_METACYCLE_ROTATION }; 

  protected double curve_scale = 1;

  protected MergeMode mode_merge = MergeMode.AUTO;
  protected RenderMode inner_mode = RenderMode.ONCURVE;
  protected RenderMode outer_mode = RenderMode.ONCURVE;
  protected RenderMode outish_mode = RenderMode.ONCURVE;
  protected CurveRadiusMode curve_rmode = CurveRadiusMode.AUTO;
  protected InsideOutsideRule location_mode = InsideOutsideRule.AUTO;
  protected CurveProximityMode proximity_mode = CurveProximityMode.AUTO;

  protected MergeMode mode_merge_param = mode_merge;
  protected RenderMode inner_param = inner_mode;
  protected RenderMode outer_param = outer_mode;
  protected RenderMode outish_param = outish_mode;
  protected CurveRadiusMode curve_rmode_param = curve_rmode;
  protected InsideOutsideRule location_mode_param = location_mode;
  protected CurveProximityMode proximity_param = proximity_mode;
  
  // point_rmode is currently hardwired, not available as a user param
  PointRadiusMode point_rmode = PointRadiusMode.MODIFIED;

  protected double inner_rspread = 0.5; 
  protected double outer_rspread = 0.5; 
  protected double outish_rspread = 0.5;
  
  protected double inner_xspread = 1;  // rspread multiplier for x
  protected double outer_xspread = 1;  // rspread multiplier for x
  protected double outish_xspread = 1; // rspread multiplier for x
  
  protected double inner_yspread = 1;  // rspread multiplier for y
  protected double outer_yspread = 1;  // rspread multiplier for y
  protected double outish_yspread = 1; // rspread multiplier for y
  
  protected double inner_rblur = 0;
  protected double outer_rblur = 0;
  protected double outish_rblur = 0;
  
  protected double inner_ablur = 0;
  protected double outer_ablur = 0;
  protected double outish_ablur = 0;
  
  protected double spread_split = 1;
  protected double fill = 0;

  protected double curve_thickness = 0;
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
  protected double metacycle_rotation = 0; // additional (cumulative) theta offset (in cycles) for each metacycle (rotate the metacycle)

  // vars for determining inner/outer via even-odd rule
  int default_sample_count = 36000;
  int angle_bin_count = 720;
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
    else if (mode_merge == MergeMode.NONE) { } // do nothing, already separate by default
    
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
    
    if (proximity_param == CurveProximityMode.AUTO) { 
      proximity_mode = CurveProximityMode.TRANSFORMED_R;
    }
    else {
      proximity_mode = proximity_param;
    }
    
    if (DEBUG_MODES) {
      System.out.println("mode_merge: " + this.mode_merge.name());
      System.out.println("inner_mode: " + inner_mode.name());
      System.out.println("outer_mode: " + outer_mode.name());
      System.out.println("outish_mode: " + outish_mode.name());
      System.out.println("curve_radius_mode: " + this.curve_rmode.name());
      System.out.println("point_radius_mode: " + this.point_rmode.name());
      System.out.println("location_mode: " + this.location_mode.name());
      System.out.println("proximity_mode: " + this.proximity_mode.name());
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

    // and add cycle_rotation (metacycle_rotation is added in calcCurvePoint())
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
    // DO NOT clear point! need initial point from subclass.calcCurvePoint()
    if (metacycles != 1 && cycles_param == 0 && cycles_to_close > 0) {
      int metacycle_count = calcMetacycle(theta);
      if (metacycle_count > 0) { 
        double metacycle_delta = (metacycle_count * metacycle_offset) + (pow(metacycle_scale, metacycle_count)-1);
        pResult.x = pResult.x * (1 + metacycle_delta);
        pResult.y = pResult.y * (1 + metacycle_delta);
        // z unchanged
        
        // add in metacycle_rotation (cycle_rotation is added in transform())
        if (metacycle_rotation != 0) {
          // rotation around origin: 
          //    x = x*cos(t) - y*sin(t)
          //    y = x*sin(t) + y*cos(t)
          double tdelta = (metacycle_count * metacycle_rotation * M_2PI);
          double cost = Math.cos(tdelta);
          double sint = Math.sin(tdelta);
          double x = pResult.x;
          double y = pResult.y;
          pResult.x = (x * cost) - (y * sint);
          pResult.y = (x * sint) + (y * cost);
        }
      }
    }
  }  
  
  public int calcMetacycle(double theta) {
    // double metacycle_count = floor((theta + (cycles * M_PI)) / (cycles_to_close * M_2PI));
    // metacycle_count = floor((theta + (cycles * M_PI) - (cycle_rotation * M_2PI)) / (cycles_to_close * M_2PI));
    double metacycle_progress = (theta + (cycles * M_PI) - (cycle_rotation * M_2PI)) / (cycles_to_close * M_2PI);
    // need to check if metacycles calcs are very slightly too low or too high
    //    (I think due to double rounding errors?)
    if (metacycle_progress < 0) {
      if (DEBUG_INTERSECTS) { System.out.println("in transform(), metacycle < 0: " + metacycle_progress + " theta: " + theta + " theta/PI: " + theta/M_PI); }
      // see if adding epsilon will raise
      metacycle_progress += EPSILON;
    }
    else if (metacycle_progress >= metacycles) {
      // see if subtracting epsilon will lower
      if (DEBUG_INTERSECTS) { System.out.println("metacycle >= metacycles: " + metacycle_progress + " theta: " + theta + " theta/PI: " + theta/M_PI); }
      metacycle_progress -= EPSILON;
    }
    int metacycle_count = (int)floor(metacycle_progress);
    return metacycle_count;
  }


  int recalcCount = 0; 

   public void recalcCurveIntersects() {
    // System.out.println("recalcing curves");
    theta_intersects = new ArrayList<ArrayList<LinkedPolarCurvePoint>>(angle_bin_count);
    for (int i=0; i<angle_bin_count; i++) { 
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
    /*if (cycles_param == 0 && cycles_to_close > 0 && metacycles > 1) {
      sampleCount = (int)(sampleCount * metacycles);
    }
    */
    // int prev_metacycle = (int)(this.metacycles - 1);
    int prev_metacycle = -1000;
    ArrayList<LinkedPolarCurvePoint> metacycle_first_points = new ArrayList<LinkedPolarCurvePoint>((int)ceil(metacycles));  // for each metacycle, keep track of first point (for looping)
    ArrayList<LinkedPolarCurvePoint> metacycle_last_points = new ArrayList<LinkedPolarCurvePoint>((int)ceil(metacycles)); // for each metacycle, keep track of last point (for looping)
    for (int m=0; m<(int)ceil(metacycles); m++) {
      metacycle_first_points.add(null);
      metacycle_last_points.add(null);
    }

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
      int anglebin =  (int)Math.floor(((angle + M_PI)/M_2PI) * angle_bin_count);
      // catching any possible cases where angle actually reaches max atan2
      //   (actually, maybe should loop around instead to bin 0?
      if (anglebin == angle_bin_count) { anglebin--; } 
      tsects = theta_intersects.get(anglebin);

      // still rotating through same bin, merge results
      if (prev_tsects == tsects) {  
        // try ignoring for now -- should try averaging later?
      }
      else {
        // tsects.add(r);
        LinkedPolarCurvePoint point = new LinkedPolarCurvePoint(r, angle);
        point.bin = anglebin;
        int metacycle_count = calcMetacycle(theta);
        if (metacycle_count == prev_metacycle) { // still in same metacycle
          point.prev = prev_point;
          if (prev_point != null) { prev_point.next = point; }
        }
        else {
          if (DEBUG_INTERSECTS) { System.out.println("new metacycle: " + metacycle_count + ", bin: " + anglebin + ", radius: " + r + ", angle: " + angle);  }
          
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

    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN PRE-PRE:"); printBin(0); System.out.println("MAXBIN PRE-PRE:");     printBin(angle_bin_count-1); }

    // for interpolation, want to close metacycle loop by redoing first_point.prev as last_point, and last_point.next as first_point
    for (int m=0; m<metacycle_last_points.size(); m++) {
      LinkedPolarCurvePoint first_meta_point = metacycle_first_points.get(m);
      LinkedPolarCurvePoint last_meta_point = metacycle_last_points.get(m);

      if (first_meta_point.bin == last_meta_point.bin) {
        if (DEBUG_INTERSECTS) { System.out.println("got first/last metacycle points in same bin: " + first_meta_point.bin); }
        // if first and last point are in same bin, replace last_meta_point with last_meta_point.prev
        tsects = theta_intersects.get(first_meta_point.bin);
        tsects.remove(last_meta_point);  
        last_meta_point = last_meta_point.prev;
        // metacycle_last_points.set(m, last_meta_point); // shouldn't be necessary, metacycle_last_points not used again
        
      }
      first_meta_point.prev = new LinkedPolarCurvePoint(last_meta_point.radius, last_meta_point.angle);
      first_meta_point.prev.bin = last_meta_point.bin;
      last_meta_point.next = new LinkedPolarCurvePoint(first_meta_point.radius, first_meta_point.angle);
      last_meta_point.next.bin = first_meta_point.bin;
      
      if (first_meta_point.next.bin == first_meta_point.bin || first_meta_point.prev.bin == first_meta_point.bin) {
        System.out.println("got first point prev/next in same bin: " + first_meta_point.bin);
      }
      else if (last_meta_point.next.bin == last_meta_point.bin || last_meta_point.prev.bin == last_meta_point.bin) {
        System.out.println("got last point prev/next in same bin: " + last_meta_point.bin);
      }
    }

    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN PRE:"); printBin(0); System.out.println("MAXBIN PRE:");     printBin(angle_bin_count-1); }
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
    ArrayList<LinkedPolarCurvePoint> maxBin = theta_intersects.get(angle_bin_count-1);
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
    
    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN MID:"); printBin(0); System.out.println("MAXBIN MID:");     printBin(angle_bin_count-1); }
    
    // special-casing of first and last anglebin if they are the same bin:
    ///   want to simulate rotating through same bin to merge "duplicate" intersections
    //    if first and last bin are same, would have merged results, so remove last one 
    if (firstbin > 0 && lastbin > 0 && firstbin == lastbin) {
      if (DEBUG_INTERSECTS) { System.out.println("WARNING: firstbin and lastbin are same: " + firstbin); }
      if (DEBUG_INTERSECTS) { printBin(firstbin); } 
      // tsects = theta_intersects.get(firstbin);
      // remove last result (if start doing averaging, should remove but add to average for first result)
      // tsects.remove(tsects.size()-1);
    }
    
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

    if (this.DEBUG_INTERSECTS) { System.out.println("MINBIN POST:"); printBin(0); System.out.println("MAXBIN POST:");     printBin(angle_bin_count-1); }
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
    // (5. apply mask modes?)
    // 
   */
public void renderByMode(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount, XYZPoint calcPoint) {
    PointState pstate = PointState.INSIDE;  
    double xin, yin, rin, tin;
    double xcalc, ycalc, rcalc, tcalc;
    double xcurve, ycurve, rcurve, tcurve;
    double xout, yout, rout, tout;

    // input point
    xin = pAffineTP.x;
    yin = pAffineTP.y;
    tin = atan2(yin, xin);  // atan2 range is [-PI, PI], so covers 2PI, or 1 cycle
    rin = sqrt((xin  * xin) + (yin * yin));
    if (point_rmode == PointRadiusMode.MODIFIED) { rin = rin * spread_split; }
    
    double rinx, riny;
    
    // input point mapped to curve
    xcalc = calcPoint.x;
    ycalc = calcPoint.y;
    // double t = curvePoint.getPrecalcAtanYX();
    tcalc = atan2(ycalc, xcalc);
    // double r = curvePoint.getPrecalcSqrt();
    rcalc = sqrt(xcalc*xcalc + ycalc*ycalc);
    if (fill != 0) {
      rcalc = rcalc + (fill * (pContext.random() - 0.5));
      xcalc = rcalc * cos(tcalc);
      ycalc = rcalc * sin(tcalc);
    }

    ArrayList<LinkedPolarCurvePoint> tsects;

    int anglebin =  (int)Math.floor(((tin + M_PI)/M_2PI) * angle_bin_count);
 
    if (curve_rmode == CurveRadiusMode.RAW_SAMPLING_BIN || curve_rmode == CurveRadiusMode.INTERPOLATED_SAMPLING_BIN) {
      if (anglebin == angle_bin_count) {  // catching any possible cases where tin actually reaches max atan2
        anglebin--; 
      } 
      tsects = theta_intersects.get(anglebin);
    }
    else  {  // curve_rmode == CurveRadiusMode.TRANSFORMED_INPUT
      tsects = unisects;
      unipolar.radius = rcalc;
      unipolar.angle = tcalc;
    }

    // count curve intersections with a ray from origin that passes through input point
    int shorter = 0;
    int longer = 0;
    double longest = Double.NEGATIVE_INFINITY;
    double shortest = Double.POSITIVE_INFINITY;
    double nearest = Double.NEGATIVE_INFINITY;
    double nearest_longer  = Double.POSITIVE_INFINITY;
    double nearest_shorter = Double.NEGATIVE_INFINITY;
    for (LinkedPolarCurvePoint curve : tsects) {
      if (curve.inflection) { continue; } // if inflection, don't count as an intersection
      double curve_radius; 
      if (curve_rmode == CurveRadiusMode.INTERPOLATED_SAMPLING_BIN) {
        curve_radius = curve.interpolate(tin);
      }
      else { curve_radius = curve.radius; }
      shortest = Math.min(shortest, curve_radius);
      longest = Math.max(longest, curve_radius);
      if (curve_radius <= rin) { 
        shorter++;
        if ((rin - curve_radius) < (rin - nearest_shorter)) { nearest_shorter = curve_radius; }
      }
      else { // rcurve > rpoint
        longer++; 
        if ((curve_radius - rin) < (nearest_longer - rin)) { nearest_longer = curve_radius; }
      }
    }
    if (longer == 0) { nearest = nearest_shorter; }
    else if (shorter == 0) { nearest = nearest_longer; }
    else if ((rin - nearest_shorter) < (nearest_longer - rin)) { nearest = nearest_shorter; }
    else { nearest = nearest_longer; }
    
    if (location_mode == InsideOutsideRule.EVEN_ODD) {
      // use standard Even-Odd rule (well, more standard than the modified ones below...):
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
      double spread, xspread, yspread, radial_blur, angular_blur;

      RenderMode mode = RenderMode.DEFAULT;
      if (pstate == PointState.INSIDE)  {
        mode = inner_mode;
        spread = inner_rspread;
        xspread = inner_xspread;
        yspread = inner_yspread;
        radial_blur = inner_rblur;
        angular_blur = inner_ablur;
      }
      else if (pstate == PointState.OUTSIDE) {
        mode = outer_mode;
        spread = outer_rspread;
        xspread = outer_xspread;
        yspread = outer_yspread;
        radial_blur = outer_rblur;
        angular_blur = outer_ablur;
      }
      else if (pstate == PointState.OUTISH) {
        mode = outish_mode;
        spread = outish_rspread;
        xspread = outish_xspread;
        yspread = outish_yspread;
        radial_blur = outish_rblur;
        angular_blur = outish_ablur;
      }
      //  else if (pstate == PointState.ON) {
      // currently not using ON, so should never get here
      else {
        // every point must be classified as being within one of the enum'd areas, so should never get here
        mode = RenderMode.DEFAULT;
        spread = 0;
        xspread = 0;
        yspread = 0;
        radial_blur = 0;
        angular_blur = 0;
      }
      
      if (mode == RenderMode.DEFAULT) {
        mode = RenderMode.ONCURVE;
      }
      
      tcurve = tin;
      switch (proximity_mode) {
        case TRANSFORMED_R:
          rcurve = rcalc;
          break;
        case TRANSFORMED_RT:
          rcurve = rcalc;
          tcurve = tcalc;
          break;
        case LONGEST:
          rcurve = longest;
          break;
        case NEAREST_LONGER:
          if (longer > 0) { rcurve = nearest_longer; }
          else { rcurve = longest; } // if no longer curve intersect, default to longest
          break;
        case NEAREST:
          rcurve = nearest;
          break;
        case NEAREST_SHORTER:
          if (shorter > 0) { rcurve = nearest_shorter; } 
          else { rcurve = 0; } // if no shorter curve intersect, default to 0
          break;
        case SHORTEST:
          rcurve = shortest;
          break;
        default:
          rcurve = rcalc;
          break;
      }
      // need to figure out what to do when 
      //       using nearest_longer for rcurve but have no longer points (nearest_longer = POSITIVE_INFINITY) 
      //       using nearest_shorter for rcurve but have no shorter points (nearest_shorter = NEGATIVE_INFINITY) 
      // currently these will often end up "offscreen" at +/- infinity? could also cause problems?
       
      if (mode == RenderMode.ONCURVE) {
        // override proximity mode to set rcurve to calculated
        rcurve = rcalc;
        tcurve = tcalc;
        if (curve_thickness != 0) {
          rcurve = rcurve + (curve_thickness * (pContext.random() - 0.5));
          // recalced via xcurve calc
          xcalc = rcurve * cos(tcalc);
          ycalc = rcurve * sin(tcalc);
        }
      }
      else if (radial_blur == 0 && angular_blur == 0) {
          // do nothing
      }
      else { 
        if (radial_blur != 0) {
          rin = rin + (radial_blur * longest * (pContext.random() - 0.5));
        }
        if (angular_blur != 0) {
          tin = tin + (angular_blur * M_2PI * (pContext.random() - 0.5));
        }
        xin = rin * cos(tin);
        yin = rin * sin(tin);
      }
      xcurve = rcurve * cos(tcurve);
      ycurve = rcurve * sin(tcurve);
      
      double spreadx = spread * xspread;
      double spready = spread * yspread;
      // double spready = spread / spread_ratio;
      
      // 4. pVarTP transformations, based on mode (for given INSIDE/OUTSIDE/OUTISH classification given above)
      switch (mode) {
        case ONCURVE: // ONCURVE -- no spread, place on curve
          pVarTP.x += pAmount * xcalc;
          pVarTP.y += pAmount * ycalc;
          // pVarTP.x += pAmount * xcurve;
          // pVarTP.y += pAmount * ycurve;
          break;
        case STRETCH1: // stretch1 with modal curve radius, and using angle from input point
          // default to TRANSFORMED_RT
          rinx = rcurve + ((rin - rcurve) * spreadx);
          riny = rcurve + ((rin - rcurve) * spready);
          pVarTP.x += pAmount * rinx * cos(tcurve);
          pVarTP.y += pAmount * riny * sin(tcurve);
          break;
        case STRETCH2: // STRETCH2
          // default to xcalc/ycalc (TRANSFORMED_RT)
          rinx = (rin * spreadx) - spreadx + 1;
          riny = (rin * spready) - spready + 1;
          pVarTP.x += pAmount * rinx * xcurve;
          pVarTP.y += pAmount * riny * ycurve;
          break;       
        case UNCHANGED:   // UNCHANGED -- leave in place (degenerate case of SCALE where SCALE where scale = 1
          // point is inside curve, leave in place
          pVarTP.x += pAmount * xin;
          pVarTP.y += pAmount * yin;
          break;
        case MIRROR_SWAP:  // MIRROR_SWAP (swap around origin [0,0], inspired by RosoniFunc) [and degenerate case of ROTATE where rotation = 1PI (180 degrees)
          pVarTP.x += pAmount * xin * -1;
          pVarTP.y += pAmount * yin * -1;
          break;
        case XY_SWAP:
          pVarTP.x += pAmount * yin * spreadx;
          pVarTP.y += pAmount * xin * spready;
          break;
        case ROTATE: // ROTATE 
          // rotation around origin: 
          //    x = x*cos(t) - y*sin(t)
          //    y = x*sin(t) + y*cos(t)
          // double tdelta = spread * M_2PI;
          if (true) {
          double cost = Math.cos(spreadx * M_2PI);
          double sint = Math.sin(spready * M_2PI);
          pVarTP.x += (xin * cost) - (yin * sint);
          pVarTP.y += (xin * sint) + (yin * cost);
          }
          break;
        case ROTATE_RADIAL: // ROTATE_RADIAL: similar to ROTATE, but angle increases as extend 
          // default to 
          // rotation around origin: 
          //    x = x*cos(t) - y*sin(t)
          //    y = x*sin(t) + y*cos(t)
          // double tdelta = spread * M_2PI;
          if (true) {
          double cost = Math.cos(spreadx * rcurve);
          double sint = Math.sin(spready * rcurve);
          // possibly better to normalize if can reliably get rmax?
          //   (so if rcurve = rmax, will get spreadx * M_2PI), 
          //    (otherwise rotation increases with curve_scale)
          //   double cost = Math.cos(spreadx * (rcurve/rmax) * M_2PI);
          //   double sint = Math.sin(spready * (rcurve/rmax) * M_2PI);
          pVarTP.x += (xin * cost) - (yin * sint);
          pVarTP.y += (xin * sint) + (yin * cost);
          }
          break;          
        case CURVE_XY_OFFSET: // CURVE_XY_OFFSET -- (UNCHANGED + ONCURVE) combo (with spread)
          // default to xcalc/ycalc (TRANSFORMED_RT)
          pVarTP.x += pAmount * (xcurve + (spreadx * xin));
          pVarTP.y += pAmount * (ycurve + (spready * yin));
          break;
        case CURVE_RADIAL_OFFSET: // CURVE_RADIUS_OFFSET -- offset with curve intersect point at input angle, with: P' = C + P
          // default to longest? / tin (LONGEST)
          // default to longest?
          double rx = rcurve * cos(tcurve);
          double ry = rcurve * sin(tcurve);
          pVarTP.x += pAmount * ((xin * spreadx) + rx);
          pVarTP.y += pAmount * ((yin * spready) + ry);
          // may want to add a radius clamp -- beyond clamp, don't swap?
          //   if (rin > clamp) { leave in place }
          break;
        case CURVE_SWAP:  // CURVE_SWAP swap around curve intersect point, with P' = C + (C-P)
          // default to longest? / tin (LONGEST)
          // need to figure out what to do when 
          //       using nearest_longer for rcurve but have no longer points (nearest_longer = POSITIVE_INFINITY) 
          //       using nearest_shorter for rcurve but have no shorter points (nearest_shorter = NEGATIVE_INFINITY) 
          // currently these will end up "offscreen" at +/- infinity?
          rinx = rcurve + ((rcurve - rin) * spreadx);
          riny = rcurve + ((rcurve - rin) * spready);
          pVarTP.x += pAmount * rinx * cos(tcurve);
          pVarTP.y += pAmount * riny * sin(tcurve);
          break;
        case HIDE: // HIDE
          pVarTP.x = xin;
          pVarTP.y = yin;
          pVarTP.doHide = true;
          break;
        case SCALE: // SCALE (inspired by Circus)
          pVarTP.x = xin * spreadx;
          pVarTP.y = yin * spready;
          break;
        case WHIRL: // WHIRL (inspired by WhorlFunc)
          // default to rcalc and tin? (TRANSFORMED_R)
          double rdiff = rcurve - rin;
          if (rdiff == 0) { rdiff = EPSILON; }
          double ax = tcurve + ((spreadx/10) / rdiff);
          double ay = tcurve + ((spready/10) / rdiff);
          double cax = cos(ax);
          double say = sin(ay);
          pVarTP.x += pAmount * rin * cax;
          pVarTP.y += pAmount * rin * say;
          break;
        case POW: // POW (inspired by Juliascope)
          // default to rcalc, tcalc? (TRANSFOMRED_RT)
          double sqr = sqr(xcurve - xin) + sqr(ycurve - yin);
          rinx = rcurve + pow(sqr, spreadx) - 1.0;
          riny = rcurve + pow(sqr, spready) - 1.0;
          pVarTP.x += pAmount * rinx * cos(tcurve);
          pVarTP.y += pAmount * riny * sin(tcurve);
          break;
        case LOOPY:  // LOOPY (inspired by loonie)
          // default to longest, tin? (LONGEST)
          // rout = sqrt((calcPoint.getPrecalcSumsq() / pAffineTP.getPrecalcSumsq()) - 1.0);
          rout = sqrt(((xcurve*xcurve + ycurve*ycurve)/ pAffineTP.getPrecalcSumsq()) - 1.0);
          pVarTP.x += pAmount * rout * spreadx * xin;
          pVarTP.y += pAmount * rout * spready * yin;
          break;          
        case STRETCH3: // STRETCH3
          // default xcalc/ycalc (TRANSFORM_RT)
          // or posibly LONGEsT
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (xcurve<0) { xin = xin * -1; }
          if (ycurve<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (xcurve - (spreadx * (xcurve-xin)));
          pVarTP.y += pAmount * (ycurve - (spready * (ycurve-yin)));
          break;          
        case STRETCH4: // STRETCH4
          // default xcalc/ycalc? (TRANSFORM_RT)
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (xcurve<0) { xin = xin * -1; }
          if (ycurve<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (xcurve - (spreadx * xin));
          pVarTP.y += pAmount * (ycurve - (spready * yin));
          break;
        case STRETCH5: // STRETCH5
          // default xcalc/ycalc? (TRANSFORM_RT)
          rinx = (0.5 * rin) + spreadx;
          riny = (0.5 * rin) + spready;
          pVarTP.x += pAmount * rinx * xcurve;
          pVarTP.y += pAmount * riny * ycurve;
          break;
        case STRETCH7: // STRETCH7 -- similar to 3, different sign fiddling
          // default xcalc/ycalc? (TRANSFORM_RT)
          xin = Math.abs(xin);
          yin = Math.abs(yin);
          if (xcurve<0) { xin = xin * -1; }
          if (ycurve<0) { yin = yin * -1; }
          pVarTP.x += pAmount * (xcurve + (spreadx * xin));
          pVarTP.y += pAmount * (ycurve + (spready * yin));
          break;
        case STRETCH8: // STRETCH8 -- same as mode 6, but without the sign modifications
          // default xcalc/ycalc? (TRANSFORM_RT)
          pVarTP.x += pAmount * (xcurve - (spreadx * xin));
          pVarTP.y += pAmount * (ycurve - (spready * yin));
          break;
        case STRETCH9: // STRETCH9
          pVarTP.x += pAmount * rin * cos(tcalc) * spreadx;
          pVarTP.y += pAmount * rin * sin(tcalc) * spready;
          break;
        default:  // if mode specified has no definition, just leave on curve
          pVarTP.x += pAmount * xcalc;
          pVarTP.y += pAmount * ycalc;
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
                          inner_rspread, outer_rspread, outish_rspread, 
                          inner_xspread, outer_xspread, outish_xspread, 
                          inner_yspread, outer_yspread, outish_yspread, 
                          inner_rblur, outer_rblur, outish_rblur, 
                          inner_ablur, outer_ablur, outish_ablur, 
                          spread_split,
                          cycles_param, cycle_rotation, 
                          curve_thickness, fill, 
                          proximity_param.getIntegerMode(), curve_rmode_param.getIntegerMode(), location_mode_param.getIntegerMode(), 
                          angle_bin_count, 
                          metacycles, metacycle_offset, metacycle_scale, metacycle_rotation
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
    else if (PARAM_INNER_RSPREAD.equalsIgnoreCase(pName))
      inner_rspread = pValue;
    else if (PARAM_OUTER_RSPREAD.equalsIgnoreCase(pName))
      outer_rspread = pValue;
    else if (PARAM_OUTISH_RSPREAD.equalsIgnoreCase(pName))
      outish_rspread = pValue;    
    else if (PARAM_INNER_XSPREAD.equalsIgnoreCase(pName))
      inner_xspread = pValue;
    else if (PARAM_OUTER_XSPREAD.equalsIgnoreCase(pName))
      outer_xspread = pValue;
    else if (PARAM_OUTISH_XSPREAD.equalsIgnoreCase(pName))
      outish_xspread = pValue;  
    else if (PARAM_INNER_YSPREAD.equalsIgnoreCase(pName))
      inner_yspread = pValue;
    else if (PARAM_OUTER_YSPREAD.equalsIgnoreCase(pName))
      outer_yspread = pValue;
    else if (PARAM_OUTISH_YSPREAD.equalsIgnoreCase(pName))
      outish_yspread = pValue;
    else if (PARAM_INNER_RBLUR.equalsIgnoreCase(pName))
      inner_rblur = pValue;
    else if (PARAM_OUTER_RBLUR.equalsIgnoreCase(pName))
      outer_rblur = pValue;
    else if (PARAM_OUTISH_RBLUR.equalsIgnoreCase(pName))
      outish_rblur = pValue;
    else if (PARAM_INNER_ABLUR.equalsIgnoreCase(pName))
      inner_ablur = pValue;
    else if (PARAM_OUTER_ABLUR.equalsIgnoreCase(pName))
      outer_ablur = pValue;
    else if (PARAM_OUTISH_ABLUR.equalsIgnoreCase(pName))
      outish_ablur = pValue;    
    
    else if (PARAM_SPREAD_SPLIT.equalsIgnoreCase(pName))
      spread_split = pValue;
    else if (PARAM_CYCLES.equalsIgnoreCase(pName))
      cycles_param = abs(pValue);
    else if (PARAM_CYCLE_ROTATION.equalsIgnoreCase(pName))
      cycle_rotation = pValue; 
    else if (PARAM_CURVE_THICKNESS.equalsIgnoreCase(pName))
      curve_thickness = pValue;
    else if (PARAM_FILL.equalsIgnoreCase(pName))
      fill = pValue;
    else if (PARAM_CURVE_RADIUS_MODE.equalsIgnoreCase(pName)) {
      this.curve_rmode_param = CurveRadiusMode.get((int)floor(pValue));
      if (curve_rmode_param == null) { curve_rmode_param = CurveRadiusMode.AUTO; }
    }
    else if (PARAM_PROXIMITY_MODE.equalsIgnoreCase(pName)) {
      proximity_param = CurveProximityMode.get((int)floor(pValue));
      if (proximity_param == null) { proximity_param = CurveProximityMode.AUTO; }
    }
    else if (PARAM_LOCATION_CLASSIFIER.equalsIgnoreCase(pName)) {
      this.location_mode_param = InsideOutsideRule.get((int)floor(pValue));
      if (location_mode_param == null) { location_mode_param = InsideOutsideRule.AUTO; }
    }
    else if (PARAM_ANGLE_BINS.equalsIgnoreCase(pName)) {
      this.angle_bin_count = (int)abs(ceil(pValue));
    }
    else if (PARAM_METACYCLES.equalsIgnoreCase(pName)) {
      metacycles = abs(pValue);
      if (abs(metacycles) < 0.01) { metacycles = 0.01; }
    }
    else if (PARAM_METACYCLE_OFFSET.equalsIgnoreCase(pName))
      metacycle_offset = pValue;   
    else if (PARAM_METACYCLE_SCALE.equalsIgnoreCase(pName))
      metacycle_scale = pValue;   
    else if (PARAM_METACYCLE_ROTATION.equalsIgnoreCase(pName))
      metacycle_rotation = pValue; 
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "abstract_polar_curve";
  }


}
