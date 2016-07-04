/*
  JWildfire - an image and animation processor written in Java 
  Copyright (C) 1995-2014 Andreas Maschke

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

import static org.jwildfire.base.mathlib.MathLib.pow;

import org.jwildfire.base.Tools;
import org.jwildfire.base.mathlib.MathLib;
import org.jwildfire.create.tina.base.Flame;
import org.jwildfire.create.tina.variation.RessourceManager;
import org.jwildfire.image.Pixel;
import org.jwildfire.image.SimpleImage;

public class GammaCorrectionFilter {
  private final Flame flame;
  private int vibInt;
  private int inverseVibInt;
  private double gamma;
  private double sclGamma;
  private SimpleImage bgImage;
  private boolean withAlpha;
  private double modSaturation;
  private final int rasterWidth, rasterHeight;
  private final double alphaScale;
  private final int oversample;
  private boolean binary_transparency;
  private double luminosity_threshold;

  public static class ColorF {
    public double r, g, b;
  }

  public static class ColorI {
    public int r, g, b;
  }

  public GammaCorrectionFilter(Flame pFlame, boolean pWithAlpha, int pRasterWidth, int pRasterHeight) {
    flame = pFlame;
    withAlpha = pWithAlpha;
    rasterWidth = pRasterWidth;
    rasterHeight = pRasterHeight;
    alphaScale = 1.0 - MathLib.atan(3.0 * (pFlame.getForegroundOpacity() - 1.0)) / 1.25;
    oversample = pFlame.getSpatialOversampling();

    initFilter();
  }

  private void initFilter() {
    gamma = (flame.getGamma() == 0.0) ? flame.getGamma() : 1.0 / flame.getGamma();

    vibInt = (int) (flame.getVibrancy() * 256.0 + 0.5);
    if (vibInt < 0) {
      vibInt = 0;
    }
    else if (vibInt > 256) {
      vibInt = 256;
    }
    inverseVibInt = 256 - vibInt;

    sclGamma = 0.0;
    if (flame.getGammaThreshold() != 0.0) {
      sclGamma = pow(flame.getGammaThreshold(), gamma - 1);
    }

    if (flame.getBGImageFilename().length() > 0) {
      try {
        bgImage = (SimpleImage) RessourceManager.getImage(flame.getBGImageFilename());
        if (bgImage.getImageWidth() < 2 || bgImage.getImageHeight() < 2) {
          bgImage = null;
        }
      }
      catch (Exception ex) {
        ex.printStackTrace();
      }
    }

    modSaturation = flame.getSaturation() - 1.0;
    if (modSaturation < -1.0)
      modSaturation = -1.0;
        
    // binary transparency
    // if binary_transparency and withAlpha, then render each pixel as either fully transparent or fully opaque
    //    (alpha channel = 0 or 255, no intermediate values)
    binary_transparency = flame.isBinaryTransparency();
    // intensity threshold
    // if intensity of incoming point is <= threshold then ignore point color, just use background
    luminosity_threshold = flame.getLuminosityThresh();
  }

  public void transformPoint(LogDensityPoint logDensityPnt, GammaCorrectedRGBPoint pRGBPoint, int pX, int pY) {
    calculateBGColor(pRGBPoint, pX, pY);
    double logScl;
    int inverseAlphaInt;

    if (logDensityPnt.intensity > 0.0) {
      double alpha;
      if (logDensityPnt.intensity <= flame.getGammaThreshold()) {
        double frac = logDensityPnt.intensity / flame.getGammaThreshold();
        alpha = (1.0 - frac) * logDensityPnt.intensity * sclGamma + frac * pow(logDensityPnt.intensity, gamma);
      }
      else {
        alpha = pow(logDensityPnt.intensity, gamma);
      }
      logScl = vibInt * alpha / logDensityPnt.intensity;
      int alphaInt = (int) (alpha * 255 * alphaScale + 0.5);
      if (alphaInt < 0)
        alphaInt = 0;
      else if (alphaInt > 255)
        alphaInt = 255;
      // alphaInt: 0 => totally transparent, 255 => totally opaque
      // inverseAlphaInt: 0 ==> totally opaque, 255 => totally transparent

      // binary_transparency overrides withAlpha and makes any populated pixels fully opaque (while leaving background fully transparent)
      if (binary_transparency) { 
        pRGBPoint.alpha = 255;
        inverseAlphaInt = 0;
      }
      else {
        pRGBPoint.alpha = withAlpha ? alphaInt : 255;
        inverseAlphaInt = 255 - alphaInt;
      }

      ColorF transfColor = applyLogScale(logDensityPnt, logScl);
      ColorI finalColor = addBackground(pRGBPoint, transfColor, inverseAlphaInt);

      pRGBPoint.red = finalColor.r;
      pRGBPoint.green = finalColor.g;
      pRGBPoint.blue = finalColor.b;
      // luminosity ranges from 0 => 1 and takes into account alpha
      double luminosity = (pRGBPoint.red + pRGBPoint.green + pRGBPoint.blue) * ((double)pRGBPoint.alpha/255.0) / (255.0 * 3.0);
      // if point luminosity is below luminosity threshold then ignore point color and use background
      //   (TODO: switch to calculating point color similarity to background color and basing thresholding on this rather than luminosity?)
      if (luminosity < luminosity_threshold) {
        pRGBPoint.red = pRGBPoint.bgRed;
        pRGBPoint.green = pRGBPoint.bgGreen;
        pRGBPoint.blue = pRGBPoint.bgBlue;
        pRGBPoint.alpha = withAlpha ? 0 : 255;
      }
    }
    else {
      pRGBPoint.red = pRGBPoint.bgRed;
      pRGBPoint.green = pRGBPoint.bgGreen;
      pRGBPoint.blue = pRGBPoint.bgBlue;
      pRGBPoint.alpha = withAlpha ? 0 : 255;
    }

    if (modSaturation != 0) {
      applyModSaturation(pRGBPoint, modSaturation);
    }
  }

  private void calculateBGColor(PointWithBackgroundColor pBGColor, int pX, int pY) {
    if (bgImage != null) {
      Pixel toolPixel = pBGColor.toolPixel;
      if (rasterWidth == bgImage.getImageWidth() * oversample && rasterHeight == bgImage.getImageHeight() * oversample) {
        toolPixel.setARGBValue(bgImage.getARGBValue(pX, pY));
        pBGColor.bgRed = toolPixel.r;
        pBGColor.bgGreen = toolPixel.g;
        pBGColor.bgBlue = toolPixel.b;
      }
      else {
        double xCoord = (double) pX * (double) (bgImage.getImageWidth() * oversample - 1) / (double) (rasterWidth - 1);
        double yCoord = (double) pY * (double) (bgImage.getImageHeight() * oversample - 1) / (double) (rasterHeight - 1);

        toolPixel.setARGBValue(bgImage.getARGBValueIgnoreBounds((int) xCoord, (int) yCoord));
        int luR = toolPixel.r;
        int luG = toolPixel.g;
        int luB = toolPixel.b;

        toolPixel.setARGBValue(bgImage.getARGBValueIgnoreBounds(((int) xCoord) + 1, (int) yCoord));
        int ruR = toolPixel.r;
        int ruG = toolPixel.g;
        int ruB = toolPixel.b;
        toolPixel.setARGBValue(bgImage.getARGBValueIgnoreBounds((int) xCoord, ((int) yCoord) + 1));
        int lbR = toolPixel.r;
        int lbG = toolPixel.g;
        int lbB = toolPixel.b;
        toolPixel.setARGBValue(bgImage.getARGBValueIgnoreBounds(((int) xCoord) + 1, ((int) yCoord) + 1));
        int rbR = toolPixel.r;
        int rbG = toolPixel.g;
        int rbB = toolPixel.b;

        double x = MathLib.frac(xCoord);
        double y = MathLib.frac(yCoord);
        pBGColor.bgRed = Tools.roundColor(Tools.blerp(luR, ruR, lbR, rbR, x, y));
        pBGColor.bgGreen = Tools.roundColor(Tools.blerp(luG, ruG, lbG, rbG, x, y));
        pBGColor.bgBlue = Tools.roundColor(Tools.blerp(luB, ruB, lbB, rbB, x, y));
      }
    }
    else {
      pBGColor.bgRed = flame.getBGColorRed();
      pBGColor.bgGreen = flame.getBGColorGreen();
      pBGColor.bgBlue = flame.getBGColorBlue();
    }
  }

  private ColorI addBackground(GammaCorrectedRGBPoint pRGBPoint, ColorF pTransfColor, int pInverseAlphaInt) {
    ColorI res = new ColorI();

    res.r = (int) (pTransfColor.r + 0.5) + ((pInverseAlphaInt * pRGBPoint.bgRed) >> 8);
    if (res.r < 0)
      res.r = 0;
    else if (res.r > 255)
      res.r = 255;

    res.g = (int) (pTransfColor.g + 0.5) + ((pInverseAlphaInt * pRGBPoint.bgGreen) >> 8);
    if (res.g < 0)
      res.g = 0;
    else if (res.g > 255)
      res.g = 255;

    res.b = (int) (pTransfColor.b + 0.5) + ((pInverseAlphaInt * pRGBPoint.bgBlue) >> 8);
    if (res.b < 0)
      res.b = 0;
    else if (res.b > 255)
      res.b = 255;
    return res;
  }

  final static double ALPHA_RANGE = 256.0;

  private ColorF addBackgroundF(PointWithBackgroundColor pBGColor, ColorF pTransfColor, double pInverseAlphaInt) {
    ColorF res = new ColorF();

    res.r = pTransfColor.r + (pInverseAlphaInt * pBGColor.bgRed) / ALPHA_RANGE;
    if (res.r < 0.0)
      res.r = 0.0;

    res.g = pTransfColor.g + (pInverseAlphaInt * pBGColor.bgGreen) / ALPHA_RANGE;
    if (res.g < 0.0)
      res.g = 0.0;

    res.b = pTransfColor.b + (pInverseAlphaInt * pBGColor.bgBlue) / ALPHA_RANGE;
    if (res.b < 0.0)
      res.b = 0.0;
    return res;
  }

  private ColorF applyLogScale(LogDensityPoint pLogDensityPnt, double pLogScl) {
    ColorF res = new ColorF();
    double rawRed, rawGreen, rawBlue;
    if (inverseVibInt > 0) {
      rawRed = pLogScl * pLogDensityPnt.red + inverseVibInt * pow(pLogDensityPnt.red, gamma);
      rawGreen = pLogScl * pLogDensityPnt.green + inverseVibInt * pow(pLogDensityPnt.green, gamma);
      rawBlue = pLogScl * pLogDensityPnt.blue + inverseVibInt * pow(pLogDensityPnt.blue, gamma);
    }
    else {
      rawRed = pLogScl * pLogDensityPnt.red;
      rawGreen = pLogScl * pLogDensityPnt.green;
      rawBlue = pLogScl * pLogDensityPnt.blue;
    }
    res.r = rawRed;
    res.g = rawGreen;
    res.b = rawBlue;
    return res;
  }

  private void applyModSaturation(GammaCorrectedRGBPoint pRGBPoint, double currModSaturation) {
    HSLRGBConverter hslrgbConverter = pRGBPoint.hslrgbConverter;
    hslrgbConverter.fromRgb(pRGBPoint.red / COLORSCL, pRGBPoint.green / COLORSCL, pRGBPoint.blue / COLORSCL);
    hslrgbConverter.fromHsl(hslrgbConverter.getHue(), hslrgbConverter.getSaturation() + currModSaturation, hslrgbConverter.getLuminosity());
    pRGBPoint.red = Tools.roundColor(hslrgbConverter.getRed() * COLORSCL);
    pRGBPoint.green = Tools.roundColor(hslrgbConverter.getGreen() * COLORSCL);
    pRGBPoint.blue = Tools.roundColor(hslrgbConverter.getBlue() * COLORSCL);
  }

  private static final double COLORSCL = 255.0;

  public void transformPointHDR(LogDensityPoint logDensityPnt, GammaCorrectedHDRPoint pHDRPoint, int pX, int pY) {
    calculateBGColor(pHDRPoint, pX, pY);

    double logScl;
    double inverseAlphaInt;
    if (logDensityPnt.intensity > 0.0) {
      double alpha;
      if (logDensityPnt.intensity <= flame.getGammaThreshold()) {
        double frac = logDensityPnt.intensity / flame.getGammaThreshold();
        alpha = (1.0 - frac) * logDensityPnt.intensity * sclGamma + frac * pow(logDensityPnt.intensity, gamma);
      }
      else {
        alpha = pow(logDensityPnt.intensity, gamma);
      }

      logScl = vibInt * alpha / logDensityPnt.intensity;
      double alphaInt = alpha * ALPHA_RANGE;
      if (alphaInt < 0.0)
        alphaInt = 0.0;
      else if (alphaInt > ALPHA_RANGE)
        alphaInt = ALPHA_RANGE;
      inverseAlphaInt = ALPHA_RANGE - alphaInt;

      ColorF transfColor = applyLogScale(logDensityPnt, logScl);
      ColorF finalColor = addBackgroundF(pHDRPoint, transfColor, inverseAlphaInt);

      pHDRPoint.red = (float) (finalColor.r / ALPHA_RANGE);
      pHDRPoint.green = (float) (finalColor.g / ALPHA_RANGE);
      pHDRPoint.blue = (float) (finalColor.b / ALPHA_RANGE);
    }
    else {
      pHDRPoint.red = (float) (pHDRPoint.bgRed / ALPHA_RANGE);
      pHDRPoint.green = (float) (pHDRPoint.bgGreen / ALPHA_RANGE);
      pHDRPoint.blue = (float) (pHDRPoint.bgBlue / ALPHA_RANGE);
    }
  }

  public static class HSLRGBConverter {
    private double red, green, blue;
    private double hue, saturation, luminosity;

    public void fromHsl(double pHue, double pSaturation, double pLuminosity) {
      hue = limitVal(pHue, 0.0, 1.0);
      saturation = limitVal(pSaturation, 0.0, 1.0);
      luminosity = limitVal(pLuminosity, 0.0, 1.0);
      double v = (luminosity <= 0.5) ? (luminosity * (1.0 + saturation))
          : (luminosity + saturation - luminosity * saturation);
      if (v <= 0) {
        red = 0.0;
        green = 0.0;
        blue = 0.0;
        return;
      }
      hue *= 6.0;
      if (hue < 0.0)
        hue = 0.0;
      else if (hue > 6.0)
        hue = 6.0;
      double y = luminosity + luminosity - v;
      double x = y + (v - y) * (hue - (int) hue);
      double z = v - (v - y) * (hue - (int) hue);

      switch ((int) hue) {
        case 0:
          red = v;
          green = x;
          blue = y;
          break;
        case 1:
          red = z;
          green = v;
          blue = y;
          break;
        case 2:
          red = y;
          green = v;
          blue = x;
          break;
        case 3:
          red = y;
          green = z;
          blue = v;
          break;
        case 4:
          red = x;
          green = y;
          blue = v;
          break;
        case 5:
          red = v;
          green = y;
          blue = z;
          break;
        default:
          red = v;
          green = x;
          blue = y;
          //          red = v;
          //          green = y;
          //          blue = z;
      }
    }

    public void fromRgb(double pRed, double pGreen, double pBlue) {
      hue = 1.0;
      saturation = 0.0;
      red = limitVal(pRed, 0.0, 1.0);
      green = limitVal(pGreen, 0.0, 1.0);
      blue = limitVal(pBlue, 0.0, 1.0);
      double max = Math.max(red, Math.max(green, blue));
      double min = Math.min(red, Math.min(green, blue));
      luminosity = (min + max) / 2.0;
      if (Math.abs(luminosity) <= MathLib.EPSILON)
        return;
      saturation = max - min;

      if (Math.abs(saturation) <= MathLib.EPSILON)
        return;

      saturation /= ((luminosity) <= 0.5) ? (min + max) : (2.0 - max - min);
      if (Math.abs(red - max) < MathLib.EPSILON) {
        hue = ((green == min) ? 5.0 + (max - blue) / (max - min) : 1.0 - (max - green) / (max - min));
      }
      else {
        if (Math.abs(green - max) < MathLib.EPSILON) {
          hue = ((blue == min) ? 1.0 + (max - red) / (max - min) : 3.0 - (max - blue) / (max - min));
        }
        else {
          hue = ((red == min) ? 3.0 + (max - green) / (max - min) : 5.0 - (max - red) / (max - min));
        }
      }
      hue /= 6.0;
    }

    private double limitVal(double pValue, double pMin, double pMax) {
      return pValue < pMin ? pMin : pValue > pMax ? pMax : pValue;
    }

    public double getRed() {
      return red;
    }

    public double getGreen() {
      return green;
    }

    public double getBlue() {
      return blue;
    }

    public double getHue() {
      return hue;
    }

    public double getSaturation() {
      return saturation;
    }

    public double getLuminosity() {
      return luminosity;
    }
  }

}
