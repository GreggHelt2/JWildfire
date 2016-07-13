/*
  JWildfire - an image and animation processor written in Java 
  Copyright (C) 1995-2015 Andreas Maschke

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
package org.jwildfire.create.tina.render;

import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.log;
import static org.jwildfire.base.mathlib.MathLib.log10;
import static org.jwildfire.base.mathlib.MathLib.sin;

import org.jwildfire.base.Tools;
import org.jwildfire.base.mathlib.MathLib;
import org.jwildfire.create.tina.base.Flame;
import org.jwildfire.create.tina.base.raster.AbstractRaster;
import org.jwildfire.create.tina.base.raster.RasterPoint;
import org.jwildfire.create.tina.random.AbstractRandomGenerator;
import org.jwildfire.create.tina.random.MarsagliaRandomGenerator;
import org.jwildfire.create.tina.swing.ChannelMixerCurves;

public class LogDensityFilter extends FilterHolder {
  private final ColorFunc colorFunc;

  private AbstractRaster raster;
  private int rasterWidth, rasterHeight, rasterSize;
  private final int PRECALC_LOG_ARRAY_SIZE = 512;
  private double precalcLogArray[];
  private double k1, k2;
  private double motionBlurScl;
  private final AbstractRandomGenerator randGen;
  private final boolean jitter;
  private final int colorOversampling;

  private double logda;
  
  private long transformCalls;
  private long innerLoops;

  private double maxOrigR;
  private double maxRawR;
  private double maxMapR;
  private double maxTransR;
  private double maxCalcR;
  private double maxCalcG;
  private double maxCalcB;
  private long maxOrigCount;
  private long minOrigCount; // min non-zero count
  private double minPreIntensity;
  private double maxPreIntensity;
  private double minPostIntensity;
  private double maxPostIntensity;
  private double maxIntensity;  
  private double minIntensity; // min non-zero intensity
  
  private double totalIntensity;
  private double totalNonZero;
  
  private double area;
  private boolean DEBUG = false;
  
  private DensityMapper densityMap;
  
  abstract class DensityMapper {
    abstract public double map(double count);
  }
  
  class Log10Mapper extends DensityMapper {
    public double map(double count) {
      double mcount = count * motionBlurScl;
      return (k1 * log10(1.0 + mcount * k2)) / (flame.getWhiteLevel() * mcount);
    }
  }
  
  class LogNaturalMapper extends DensityMapper {
    public double map(double count) {
      double mcount = count * motionBlurScl;
      return (k1 * log(1.0 + mcount * k2)) / (flame.getWhiteLevel() * mcount);
    }
  }
    
  class LinearMapper extends DensityMapper {
    public double map(double count) {
      double mcount = count * motionBlurScl;
      return (k1 * (1.0 + mcount * k2)) / (flame.getWhiteLevel() * mcount);
    }
  }
  
  class LogBaseNMapper extends DensityMapper {
    double base;
    double logb;

    public LogBaseNMapper(double b) {
      base = b;
      logb = log(b);
    }
    // log_b(x) = log_k(x) / log_k(b)  ==> log(x)/log(b) ==> log(x)/logb
    public double map(double count) {
      double mcount = count * motionBlurScl;
      return (k1 * (log(1.0 + mcount * k2)/logb)) / (flame.getWhiteLevel() * mcount);
    }
  }

  public LogDensityFilter(Flame pFlame, AbstractRandomGenerator pRandGen) {
    super(pFlame);
    colorFunc = pFlame.getChannelMixerMode().getColorFunc(pFlame, pRandGen);
    motionBlurScl = flame.getMotionBlurLength() <= 0 ? 1.0 : 1.0 / (flame.getMotionBlurLength() + 1.0);
    jitter = pFlame.isSampleJittering();
    colorOversampling = jitter ? pFlame.getColorOversampling() : 1;
    if (jitter) {
      randGen = new MarsagliaRandomGenerator();
      randGen.randomize(pFlame.hashCode());
    }
    else {
      randGen = null;
    }
    double log_density_base = flame.getLogDensityBase();
    if (log_density_base == 10.0) {
      densityMap = new Log10Mapper();
    }
    else if (Math.abs(log_density_base - Math.E) < 0.01) {
      densityMap = new LogNaturalMapper();
    }
    else { 
      densityMap = new LogBaseNMapper(log_density_base);
    }
  }
  
  public ColorFunc getColorFunction() {
    return colorFunc;
  }
  
  public void clearCounts() {
    transformCalls = 0;
    innerLoops = 0; 
    
    maxOrigR = 0;
    maxRawR = 0;
    maxMapR = 0;
    maxTransR = 0;
    maxCalcR = 0;
    maxCalcG = 0;
    maxCalcB = 0;
    minOrigCount = 100000000;
    maxOrigCount = 0;
    minPreIntensity = 100000000;
    maxPreIntensity = 0;
    minPreIntensity = 100000000;
    maxPostIntensity = 0;
    minIntensity = 100000000;
    maxIntensity = 0; 
    totalIntensity = 0;
    totalNonZero = 0;
  }
  
  public void reportCounts() {
    /*
            System.out.println("maxOrigR: " + maxOrigR);
    System.out.println("maxFinalR: " + maxCalcR);    
    System.out.println("maxFinalG: " + maxCalcG);    
    System.out.println("maxFinalB: " + maxCalcB);    
    System.out.println("minOrigCount:" + minOrigCount);
    System.out.println("maxOrigCount:" + maxOrigCount);
    System.out.println("minPreIntensity: " + minPreIntensity);
    System.out.println("maxPreIntensity: " + maxPreIntensity);
    System.out.println("minPostIntensity: " + minPostIntensity);
    System.out.println("maxPostIntensity: " + maxPostIntensity);
            */
    double density = flame.getSampleDensity();  
    double avgIntensity = totalIntensity/totalNonZero;
    System.out.println("sample density: " + flame.getSampleDensity() + ", area: " + area + ", d*a: " + (area * density) + 
            ", log10(d*a): " + logda + ", oversample: " + oversample + ", noiseSize: " + noiseFilterSize);
    System.out.println("minIntensity: " + minIntensity + ", *lda: " + minIntensity*logda + ", *ldao: " + minIntensity*logda/oversample);
    System.out.println("maxIntensity: " + maxIntensity + ", *lda: " + maxIntensity*logda + ", *ldao " + maxIntensity*logda/oversample);
    if (totalNonZero > 0) { System.out.println("avgIntensity: " + avgIntensity + ", *lda: " + avgIntensity*logda + ", /lda: " + avgIntensity/logda); }
  }

  public void setRaster(AbstractRaster pRaster, int pRasterWidth, int pRasterHeight, int pImageWidth, int pImageHeight) {
    raster = pRaster;
    rasterWidth = pRasterWidth;
    rasterHeight = pRasterHeight;
    rasterSize = rasterWidth * rasterHeight;
    // k1 has oversample in denominator
    k1 = flame.getContrast() * 2.0 * flame.getBrightness() / (double) (oversample);
    switch (flame.getPostSymmetryType()) {
      case POINT:
        k1 /= (double) flame.getPostSymmetryOrder();
        break;
      case X_AXIS:
      case Y_AXIS:
        k1 /= 2.0;
    }
    double pixelsPerUnit = flame.getPixelsPerUnit() * flame.getCamZoom();
    // double area = ((double) pImageWidth * (double) pImageHeight) / (pixelsPerUnit * pixelsPerUnit);
    area = ((double) pImageWidth * (double) pImageHeight) / (pixelsPerUnit * pixelsPerUnit);
    logda = log10(area * flame.getSampleDensity());
    
    //      flame sample density calculated in calcDensity ==>  pSampleCount / pRasterSize * oversample;
    // therefore k2 has oversample in numerator
    k2 = 1.0 / (flame.getContrast() * area * flame.getSampleDensity());

    precalcLogArray = new double[PRECALC_LOG_ARRAY_SIZE + 1];
    for (int i = 1; i <= PRECALC_LOG_ARRAY_SIZE; i++) {
      double x = i * motionBlurScl;
      // precalcLogArray[i] = (k1 * log10(1 + x * k2)) / (flame.getWhiteLevel() * x);
      precalcLogArray[i] = densityMap.map(i);
    }
  }
  

/*
  public void transformPointSimple(LogDensityPoint pFilteredPnt, int pX, int pY) {
    pFilteredPnt.red = pFilteredPnt.green = pFilteredPnt.blue = 0;
    pFilteredPnt.intensity = 0;
    for (int px = 0; px < oversample; px++) {
      for (int py = 0; py < oversample; py++) {
        getSample(pFilteredPnt, pX * oversample + px, pY * oversample + py);
        double logScale;
        long pCount = pFilteredPnt.rp.count;
        // mods for point_count_negation
        //     since counts can go negative, need to zero out before calculating logScale
        if (pCount < 0) { pCount = 0; }
        if (pCount < precalcLogArray.length) {
          logScale = precalcLogArray[(int) pCount];
        }
        else {
          // logScale = (k1 * log10(1.0 + pCount * motionBlurScl * k2)) / (flame.getWhiteLevel() * pCount * motionBlurScl);
          logScale = densityMap.map(pCount);
        }
        if (pCount > 0) {
          if (colorFunc == ColorFunc.NULL) {
            pFilteredPnt.red += logScale * pFilteredPnt.rp.red;
            pFilteredPnt.green += logScale * pFilteredPnt.rp.green;
            pFilteredPnt.blue += logScale * pFilteredPnt.rp.blue;
          }
          else {
            final double scale = ChannelMixerCurves.FILTER_SCALE;
            double rawR = pFilteredPnt.rp.red * scale / pCount;
            double rawG = pFilteredPnt.rp.green * scale / pCount;
            double rawB = pFilteredPnt.rp.blue * scale / pCount;

            pFilteredPnt.red += logScale * colorFunc.mapRGBToR(rawR, rawG, rawB) * pCount / scale;
            pFilteredPnt.green += logScale * colorFunc.mapRGBToG(rawR, rawG, rawB) * pCount / scale;
            pFilteredPnt.blue += logScale * colorFunc.mapRGBToB(rawR, rawG, rawB) * pCount / scale;
          }
          pFilteredPnt.intensity += logScale * pCount * flame.getWhiteLevel();
        }
      }
    }
  }
  */

  private void getSample(LogDensityPoint pFilteredPnt, int pX, int pY) {
    if (jitter) {
      final double epsilon = 0.0001;
      final double radius = 0.25;
      double dr = log(randGen.random() + 0.1) + 1;
      if (dr < epsilon) {
        dr = epsilon;
      }
      else if (dr > 1.0 - epsilon) {
        dr = 1.0 - epsilon;
      }
      double da = epsilon + (randGen.random() - 2 * epsilon) * M_PI * 2.0;
      double x = dr * cos(da) * radius;
      int xi = x < 0 ? -1 : 1;
      x = MathLib.fabs(x);
      double y = dr * sin(da) * radius;
      int yi = y < 0 ? -1 : 1;
      y = MathLib.fabs(y);

      raster.readRasterPointSafe(pX, pY, pFilteredPnt.lu);
      raster.readRasterPointSafe(pX + xi, pY, pFilteredPnt.ru);
      raster.readRasterPointSafe(pX, pY + yi, pFilteredPnt.lb);
      raster.readRasterPointSafe(pX + xi, pY + yi, pFilteredPnt.rb);
      pFilteredPnt.rp.red = Tools.blerp(pFilteredPnt.lu.red, pFilteredPnt.ru.red, pFilteredPnt.lb.red, pFilteredPnt.rb.red, x, y);
      pFilteredPnt.rp.green = Tools.blerp(pFilteredPnt.lu.green, pFilteredPnt.ru.green, pFilteredPnt.lb.green, pFilteredPnt.rb.green, x, y);
      pFilteredPnt.rp.blue = Tools.blerp(pFilteredPnt.lu.blue, pFilteredPnt.ru.blue, pFilteredPnt.lb.blue, pFilteredPnt.rb.blue, x, y);
      pFilteredPnt.rp.count = Math.round(Tools.blerp(pFilteredPnt.lu.count, pFilteredPnt.ru.count, pFilteredPnt.lb.count, pFilteredPnt.rb.count, x, y));
      //      System.out.println(pFilteredPnt.rp.red + " (" + pFilteredPnt.lu.red + " " + pFilteredPnt.ru.red + " " + pFilteredPnt.lb.red + "," + pFilteredPnt.rb.red + ") at (" + x + " " + y + ")");
    }
    else {
      raster.readRasterPointSafe(pX, pY, pFilteredPnt.rp);
    }
  }

  public double calcDensity(long pSampleCount, long pRasterSize) {
    return (double) pSampleCount / (double) pRasterSize * oversample;
  }

  public double calcDensity(long pSampleCount) {
    if (rasterSize == 0) {
      throw new IllegalStateException();
    }
    return (double) pSampleCount / (double) rasterSize * oversample;
  }
  
  int pointTotal = 0;
  
  public void transformPointSimple(LogDensityPoint pFilteredPnt, int pX, int pY) {
    transformPoint(pFilteredPnt, pX, pY, true);
  }
  
  public void transformPoint(LogDensityPoint pFilteredPnt, int pX, int pY) { 
    transformPoint(pFilteredPnt, pX, pY, false);
  }
  
  public void transformPoint(LogDensityPoint pFilteredPnt, int pX, int pY, boolean simple) {
    pFilteredPnt.clear();
    pointTotal++;
    int innerLoopMax;
    int colorSampling;
    if (simple) {
      colorSampling = 1;
      innerLoopMax = oversample;
    }
    else {
      colorSampling = colorOversampling;
      if (noiseFilterSize <= 1) { innerLoopMax = oversample; }
      else { innerLoopMax = noiseFilterSize; }
    }

    for (int c = 0; c < colorSampling; c++) {
      for (int i = 0; i < innerLoopMax; i++) {
        for (int j = 0; j < innerLoopMax; j++) {
          double filterScale;
          if (noiseFilterSize <= 1 || simple) { filterScale = 1.0; }
          else { filterScale = filter[i][j]; }
          getSample(pFilteredPnt, pX * oversample + j, pY * oversample + i);
          RasterPoint rpoint = pFilteredPnt.rp;
          double rawR, rawG, rawB;
          rawR = rpoint.red;
          rawG = rpoint.green;
          rawB = rpoint.blue;
          long count = rpoint.count;
          if (count > 0) {
            double logScale;
            if (count < precalcLogArray.length) {
              logScale = precalcLogArray[(int)count];
            }
            else {
              //logScale = (k1 * log10(1.0 + count * motionBlurScl * k2)) / (flame.getWhiteLevel() * count * motionBlurScl);
              logScale = densityMap.map(count);
            }
            // double currIntensity = logScale * count * flame.getWhiteLevel() / (double) colorOversampling;
            double calcIntensity = logScale * count * flame.getWhiteLevel();
            double calcR, calcG, calcB;
            if (colorFunc == ColorFunc.NULL) {
              calcR = logScale * rawR;
              calcG = logScale * rawG;
              calcB = logScale * rawB;
            }
            else {
              double mixerScale = ChannelMixerCurves.FILTER_SCALE;
              double modR, modG, modB;
              modR = rawR * mixerScale / count;
              modG = rawG * mixerScale / count;
              modB = rawB * mixerScale / count;
              double transR, transG, transB;
              transR = colorFunc.mapRGBToR(modR, modG, modB) * count / mixerScale;
              transG = colorFunc.mapRGBToG(modR, modG, modB) * count / mixerScale;
              transB = colorFunc.mapRGBToB(modR, modG, modB) * count / mixerScale;
              calcR = logScale * transR;
              calcG = logScale * transG;
              calcB = logScale * transB;
            }
            // pFilteredPnt.red += filterScale * logScale * pFilteredPnt.rp.red / (double) colorOversampling;
            // pushing colorOversampling scaling to outside of loops
            pFilteredPnt.red   += filterScale * calcR;
            pFilteredPnt.green += filterScale * calcG;
            pFilteredPnt.blue  += filterScale * calcB;
            // pFilteredPnt.intensity += filterScale * logScale * count * flame.getWhiteLevel() / colorOversampling;
            pFilteredPnt.intensity += filterScale * calcIntensity;
            /*  
            // GAH: having trouble getting this to work in sampling loop
            // so wrote similar method with different signature that works outside loop, after 
            //    all other adjustments have been made.
            // TODO: would really like to figure out how to make this work in-loop, presumably would get better quality that way
            //        (though what I'm getting post-loop is looking pretty good)
            if (colorFunc instanceof IntensityColorFunc) {
              adjustIntensity(pFilteredPnt, calcIntensity, calcR, calcG, calcB, filterScale);
            }
            */
          }
        }
      }
    }
    pFilteredPnt.red /= colorSampling;
    pFilteredPnt.green /= colorSampling;
    pFilteredPnt.blue /= colorSampling;
    pFilteredPnt.intensity /= colorSampling;
    if (pFilteredPnt.intensity > maxIntensity) { maxIntensity = pFilteredPnt.intensity; }
    if (pFilteredPnt.intensity < minIntensity) { minIntensity = pFilteredPnt.intensity; }
    if (pFilteredPnt.intensity > 0.0 && colorFunc instanceof IntensityColorFunc) {
      adjustIntensity(pFilteredPnt);
    }
    pFilteredPnt.clip();
    totalIntensity += pFilteredPnt.intensity;
    transformCalls++;
  }
  
  protected void adjustIntensity(LogDensityPoint pFilteredPnt) {
    IntensityColorFunc intensityFunc = (IntensityColorFunc)colorFunc;
    double inputIntensity = pFilteredPnt.intensity;
    double inputR = pFilteredPnt.red;
    double inputG = pFilteredPnt.green;
    double inputB = pFilteredPnt.blue;
    
    double modScale = logda * ChannelMixerCurves.FILTER_SCALE * 10;
    double scaledIntensity = inputIntensity * modScale;
    double modScaledIntensity = intensityFunc.mapIntensity(scaledIntensity);

    double intensity_ratio;
    // hack to work around bug in > 25000 always mapping to 0
    if (scaledIntensity > 25000 && modScaledIntensity == 0.0) {
      modScaledIntensity = intensityFunc.mapIntensity(scaledIntensity);    // strictly for debugging use
      modScaledIntensity = 25000;
    }  
    double modIntensity = modScaledIntensity / modScale;
    if (modScaledIntensity == 0.0) {
      intensity_ratio = 0.0;
    }
    else { 
      intensity_ratio = modIntensity / inputIntensity;
    }
    if (modIntensity < minPostIntensity) { minPostIntensity = modIntensity; }
    if (modIntensity > maxPostIntensity) { maxPostIntensity = modIntensity; }
    
    double addR = inputR * (intensity_ratio - 1);  // already added 1 * inputR (in transform call), so now add/subtract diff ratio
    double addG = inputG * (intensity_ratio - 1);
    double addB = inputB * (intensity_ratio - 1);
    double addIntensity = modIntensity - inputIntensity; // already added inputIntensity (in transorm call), so now add/subtract diff
    
    pFilteredPnt.red +=  addR;
    pFilteredPnt.green += addG;
    pFilteredPnt. blue += addB;
    pFilteredPnt.intensity += addIntensity;
  }
   
  /*
  protected void adjustIntensity(LogDensityPoint pFilteredPnt, double inputIntensity, double inputR, double inputG, double inputB, double filterScale) {
    IntensityColorFunc intensityFunc = (IntensityColorFunc)colorFunc;
    double modScale = logda * ChannelMixerCurves.FILTER_SCALE * 10;
    double scaledIntensity = inputIntensity * modScale;
    double modScaledIntensity = intensityFunc.mapIntensity(scaledIntensity);
    // if (modScaledIntensity > 25000) { modScaledIntensity = 25000; }
    double modIntensity = modScaledIntensity / modScale;
    double intensity_ratio;
    if (inputIntensity == 0 || modIntensity == 0) {
      intensity_ratio = 0;
      // modScaledIntensity = intensityFunc.mapIntensity(scaledIntensity);    // strictly for debugging use
    }  // shouldn't happen with count > 0 ??
    else { intensity_ratio = modIntensity / inputIntensity; }
    if (modIntensity < minPostIntensity) { minPostIntensity = modIntensity; }
    if (modIntensity > maxPostIntensity) { maxPostIntensity = modIntensity; }
    
    double addR = inputR * (intensity_ratio - 1);  // already added 1 * inputR (in transform call), so now add/subtract diff ratio
    double addG = inputG * (intensity_ratio - 1);
    double addB = inputB * (intensity_ratio - 1);
    double addIntensity = modIntensity - inputIntensity; // already added inputIntensity (in transorm call), so now add/subtract diff
    
    pFilteredPnt.red +=   filterScale * addR;
    pFilteredPnt.green += filterScale * addG;
    pFilteredPnt. blue += filterScale * addB;
    pFilteredPnt.intensity += filterScale * addIntensity;
    // pFilteredPnt.intensity += addIntensity;
 
  }
  */

 
}
