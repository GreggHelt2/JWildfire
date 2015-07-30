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

import static org.jwildfire.base.mathlib.MathLib.EPSILON;
import static org.jwildfire.base.mathlib.MathLib.M_2PI;
import static org.jwildfire.base.mathlib.MathLib.M_PI;
import static org.jwildfire.base.mathlib.MathLib.cos;
import static org.jwildfire.base.mathlib.MathLib.exp;
import static org.jwildfire.base.mathlib.MathLib.fabs;
import static org.jwildfire.base.mathlib.MathLib.log;
import static org.jwildfire.base.mathlib.MathLib.sin;
import static org.jwildfire.base.mathlib.MathLib.sqrt;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.jwildfire.base.Tools;
import org.jwildfire.base.mathlib.MathLib;
import org.jwildfire.create.GradientCreator;
import org.jwildfire.create.tina.base.Constants;
import org.jwildfire.create.tina.base.DrawMode;
import org.jwildfire.create.tina.base.Flame;
import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.base.XYZProjectedPoint;
import org.jwildfire.create.tina.base.raster.AbstractRaster;
import org.jwildfire.create.tina.palette.RGBPalette;
import org.jwildfire.create.tina.palette.RenderColor;
import org.jwildfire.create.tina.random.AbstractRandomGenerator;
import org.jwildfire.create.tina.variation.FlameTransformationContext;
import org.jwildfire.create.tina.variation.RessourceManager;
import org.jwildfire.image.Pixel;
import org.jwildfire.image.SimpleImage;

public class DefaultRenderIterationState extends RenderIterationState {

  private static final long serialVersionUID = 2L;

  protected XYZPoint affineT;
  protected XYZPoint varT;
  protected XYZPoint p;
  protected XYZPoint q;
  protected XForm xf;
  protected final XYZProjectedPoint prj = new XYZProjectedPoint();
  protected PointProjector projector;
  protected final ColorProvider colorProvider;

  public DefaultRenderIterationState(AbstractRenderThread pRenderThread, FlameRenderer pRenderer, RenderPacket pPacket, Layer pLayer, FlameTransformationContext pCtx, AbstractRandomGenerator pRandGen) {
    super(pRenderThread, pRenderer, pPacket, pLayer, pCtx, pRandGen);
    projector = new DefaultPointProjector();
    if (pLayer.getGradientMapFilename() != null && pLayer.getGradientMapFilename().length() > 0) {
      colorProvider = new GradientMapColorProvider(pLayer.getGradientMapFilename());
    }
    else {
      colorProvider = pLayer.isSmoothGradient() ? new SmoothColorProvider() : new DefaultColorProvider();
    }

    Flame flame = pPacket.getFlame();
    switch (flame.getPostSymmetryType()) {
      case POINT: {
        if (flame.getPostSymmetryOrder() > 1) {
          int order = flame.getPostSymmetryOrder() <= 64 ? flame.getPostSymmetryOrder() : 64;
          projector = new PointSymmetryProjector(projector, order, flame.getPostSymmetryCentreX(), flame.getPostSymmetryCentreY());
        }
        break;
      }
      case X_AXIS:
        projector = new XAxisSymmetryProjector(projector, flame.getPostSymmetryDistance(), flame.getPostSymmetryCentreX(), flame.getPostSymmetryCentreY(), flame.getPostSymmetryRotation());
        break;
      case Y_AXIS:
        projector = new YAxisSymmetryProjector(projector, flame.getPostSymmetryDistance(), flame.getPostSymmetryCentreX(), flame.getPostSymmetryCentreY(), flame.getPostSymmetryRotation());
        break;
    }
  }

  public void init() {
  }

  public void preFuseIter() {
    affineT = new XYZPoint(); // affine part of the transformation
    varT = new XYZPoint(); // complete transformation
    p = new XYZPoint();
    q = new XYZPoint();
    p.x = 2.0 * randGen.random() - 1.0;
    p.y = 2.0 * randGen.random() - 1.0;
    p.z = 0.0;
    p.color = randGen.random();
    p.modGamma = 0.0;
    p.modContrast = 0.0;
    p.modSaturation = 0.0;

    xf = layer.getXForms().get(0);
    transformPoint();
    for (int i = 0; i <= Constants.INITIAL_ITERATIONS; i++) {
      xf = selectNextXForm(xf);
      if (xf == null) {
        xf = layer.getXForms().get(0);
        return;
      }
      transformPoint();
    }
  }

  protected XForm selectNextXForm(XForm pFrom) {
    return pFrom.getNextAppliedXFormTable()[randGen.random(Constants.NEXT_APPLIED_XFORM_TABLE_SIZE)];
  }

  public void validateState() {
    if (Double.isInfinite(p.x) || Double.isInfinite(p.y) || Double.isInfinite(p.z) || Double.isNaN(p.x) || Double.isNaN(p.y) || Double.isNaN(p.z) || xf == null) {
      preFuseIter();
    }
  }

  public void iterateNext() {
    xf = selectNextXForm(xf);
    transformPoint();
    if (xf.getDrawMode() == DrawMode.HIDDEN)
      return;
    else if ((xf.getDrawMode() == DrawMode.OPAQUE) && (randGen.random() > xf.getOpacity()))
      return;
    List<XForm> finalXForms = layer.getFinalXForms();
    if (finalXForms.size() > 0) {
      applyFinalTransforms(finalXForms);
    }
    else {
      applyEmptyFinalTransform();
    }
    projector.projectPoint(q);
  }

  public void iterateNext(List<RenderSlice> pSlices, double pThicknessMod, int pTicknessSamples) {
    int nextXForm = randGen.random(Constants.NEXT_APPLIED_XFORM_TABLE_SIZE);
    xf = xf.getNextAppliedXFormTable()[nextXForm];
    if (xf == null) {
      return;
    }
    transformPoint();
    if (xf.getDrawMode() == DrawMode.HIDDEN)
      return;
    else if ((xf.getDrawMode() == DrawMode.OPAQUE) && (randGen.random() > xf.getOpacity()))
      return;
    List<XForm> finalXForms = layer.getFinalXForms();
    if (finalXForms.size() > 0) {
      applyFinalTransforms(finalXForms);
    }
    else {
      applyEmptyFinalTransform();
    }

    if (pThicknessMod > EPSILON && pTicknessSamples > 0) {
      XYZPoint w = new XYZPoint();
      w.assign(q);
      try {
        for (int i = 0; i < pTicknessSamples; i++) {
          addNoise(w, q, pThicknessMod);
          List<AbstractRaster> sliceRasters = collectSliceRasters(q, pSlices);
          for (AbstractRaster sliceRaster : sliceRasters) {
            raster = sliceRaster;
            projector.projectPoint(q);
          }
        }
      }
      finally {
        q.assign(w);
      }
    }
    else {
      if (setupSliceRaster(q, pSlices)) {
        projector.projectPoint(q);
      }
    }
  }

  private void addNoise(XYZPoint pSrc, XYZPoint pDst, double pRadius) {
    double angle = randGen.random() * 2 * M_PI;
    double sina = sin(angle);
    double cosa = cos(angle);
    angle = randGen.random() * M_PI;
    double sinb = sin(angle);
    double cosb = cos(angle);
    pDst.x = pSrc.x + pRadius * sinb * cosa;
    pDst.y = pSrc.y + pRadius * sinb * sina;
    pDst.z = pSrc.z + pRadius * cosb;
  }

  private boolean setupSliceRaster(XYZPoint pPoint, List<RenderSlice> pSlices) {
    int sliceIdx = getSliceIndex(pPoint, pSlices, 0);
    if (sliceIdx >= 0) {
      raster = pSlices.get(sliceIdx).getRaster();
      return true;
    }
    else {
      raster = null;
      return false;
    }
  }

  private List<AbstractRaster> collectSliceRasters(XYZPoint pPoint, List<RenderSlice> pSlices) {
    List<AbstractRaster> res = new ArrayList<AbstractRaster>();
    int sliceIdx = -1;
    while (true) {
      sliceIdx = getSliceIndex(pPoint, pSlices, sliceIdx + 1);
      if (sliceIdx >= 0) {
        res.add(pSlices.get(sliceIdx).getRaster());
      }
      else {
        break;
      }
    }
    return res;
  }

  private int getSliceIndex(XYZPoint pPoint, List<RenderSlice> pSlices, int pStartIndex) {
    double value = pPoint.z;
    for (int i = pStartIndex; i < pSlices.size(); i++) {
      RenderSlice slice = pSlices.get(i);
      if (value >= slice.getZmin() && value < slice.getZmax()) {
        return i;
      }
    }
    return -1;
  }

  protected void applyEmptyFinalTransform() {
    q.assign(p);
  }

  protected void applyFinalTransforms(List<XForm> finalXForms) {
    finalXForms.get(0).transformPoint(ctx, affineT, varT, p, q);
    for (int i = 1; i < finalXForms.size(); i++) {
      finalXForms.get(i).transformPoint(ctx, affineT, varT, q, q);
    }
  }

  protected void transformPoint() {
    xf.transformPoint(ctx, affineT, varT, p, p);
  }

  protected double plotRed, plotGreen, plotBlue, plotContribution;

  protected int plotBufferIdx = 0;

  protected PlotSample[] plotBuffer = initPlotBuffer();

  interface ColorProvider extends Serializable {
    RenderColor getColor(XYZPoint pTPoint, XYZPoint pFPoint);
  }

  private class DefaultColorProvider implements ColorProvider {
    private static final long serialVersionUID = 1L;

    @Override
    public RenderColor getColor(XYZPoint pTPoint, XYZPoint pFPoint) {
      int colorIdx = (int) (pTPoint.color * paletteIdxScl + 0.5);
      if (colorIdx < 0)
        colorIdx = 0;
      else if (colorIdx >= RGBPalette.PALETTE_SIZE)
        colorIdx = RGBPalette.PALETTE_SIZE - 1;
      return colorMap[colorIdx];
    }

  }

  private class SmoothColorProvider implements ColorProvider {
    private static final long serialVersionUID = 1L;

    private final RenderColor rc = new RenderColor();

    @Override
    public RenderColor getColor(XYZPoint pTPoint, XYZPoint pFPoint) {
      double colorIdx = pTPoint.color * paletteIdxScl;
      int lIdx = (int) colorIdx;
      double lR, lG, lB;
      if (lIdx >= 0 && lIdx < RGBPalette.PALETTE_SIZE) {
        lR = colorMap[lIdx].red;
        lG = colorMap[lIdx].green;
        lB = colorMap[lIdx].blue;
      }
      else {
        lR = lG = lB = 0.0;
      }

      double rR, rG, rB;
      int rIdx = lIdx + 1;
      if (rIdx >= 0 && rIdx < RGBPalette.PALETTE_SIZE) {
        rR = colorMap[rIdx].red;
        rG = colorMap[rIdx].green;
        rB = colorMap[rIdx].blue;
      }
      else {
        rR = rG = rB = 0.0;
      }

      double t = MathLib.frac(colorIdx);
      rc.red = Tools.lerp(lR, rR, t);
      rc.green = Tools.lerp(lG, rG, t);
      rc.blue = Tools.lerp(lB, rB, t);
      return rc;
    }
  }

  private class GradientMapColorProvider implements ColorProvider {
    private static final long serialVersionUID = 1L;
    private final SimpleImage map;

    private SimpleImage createDfltImage() {
      GradientCreator creator = new GradientCreator();
      return creator.createImage(256, 256);
    }

    private final Pixel toolPixel = new Pixel();
    private final RenderColor rc = new RenderColor();

    public GradientMapColorProvider(String pGradientMap) {
      SimpleImage image;
      try {
        image = (SimpleImage) RessourceManager.getImage(pGradientMap);
      }
      catch (Exception e) {
        e.printStackTrace();
        image = createDfltImage();
      }
      map = image;
    }

    @Override
    public RenderColor getColor(XYZPoint pTPoint, XYZPoint pFPoint) {
      double localColorAdd = layer.getGradientMapLocalColorAdd();
      double localColorMultiply = layer.getGradientMapLocalColorScale();

      double x = (pFPoint.x * (1.0 - localColorMultiply) + pFPoint.x * localColorMultiply * pFPoint.color + localColorAdd * pFPoint.color) * layer.getGradientMapHorizScale() + layer.getGradientMapHorizOffset();
      double y = (pFPoint.y * (1.0 - localColorMultiply) + pFPoint.y * localColorMultiply * pFPoint.color + localColorAdd * pFPoint.color) * layer.getGradientMapVertScale() + layer.getGradientMapVertOffset();

      double width = map.getImageWidth() - 2;
      double height = map.getImageHeight() - 2;
      double fx = MathLib.fabs(x);
      double imageX = MathLib.fmod(fx * width, width);
      if (((int) fx) % 2 == 1) {
        imageX = width - imageX;
      }
      double fy = MathLib.fabs(y);
      double imageY = MathLib.fmod(fy * height, height);
      if (((int) fy) % 2 == 1) {
        imageY = height - imageY;
      }

      toolPixel.setARGBValue(map.getARGBValueIgnoreBounds((int) imageX, (int) imageY));
      int luR = toolPixel.r;
      int luG = toolPixel.g;
      int luB = toolPixel.b;

      toolPixel.setARGBValue(map.getARGBValueIgnoreBounds(((int) imageX) + 1, (int) imageY));
      int ruR = toolPixel.r;
      int ruG = toolPixel.g;
      int ruB = toolPixel.b;
      toolPixel.setARGBValue(map.getARGBValueIgnoreBounds((int) imageX, ((int) imageY) + 1));
      int lbR = toolPixel.r;
      int lbG = toolPixel.g;
      int lbB = toolPixel.b;
      toolPixel.setARGBValue(map.getARGBValueIgnoreBounds(((int) imageX) + 1, ((int) imageY) + 1));
      int rbR = toolPixel.r;
      int rbG = toolPixel.g;
      int rbB = toolPixel.b;

      double localX = MathLib.frac(imageX);
      double localY = MathLib.frac(imageY);
      rc.red = Tools.blerp(luR, ruR, lbR, rbR, localX, localY);
      rc.green = Tools.blerp(luG, ruG, lbG, rbG, localX, localY);
      rc.blue = Tools.blerp(luB, ruB, lbB, rbB, localX, localY);
      return rc;
    }
  }

  protected void plotPoint(int xIdx, int yIdx, double intensity) {
    if (p.rgbColor) {
      plotRed = p.redColor;
      plotGreen = p.greenColor;
      plotBlue = p.blueColor;
    }
    else if (q.rgbColor) {
      plotRed = q.redColor;
      plotGreen = q.greenColor;
      plotBlue = q.blueColor;
    }
    else {
      RenderColor color = colorProvider.getColor(p, q);
      plotRed = color.red;
      plotGreen = color.green;
      plotBlue = color.blue;
    }
    transformPlotColor(p);

    plotBuffer[plotBufferIdx++].set(xIdx, yIdx, plotRed * intensity, plotGreen * intensity, plotBlue * intensity);
    if (plotBufferIdx >= plotBuffer.length) {
      applySamplesToRaster();
    }
    if (observers != null && observers.size() > 0) {
      for (IterationObserver observer : observers) {
        observer.notifyIterationFinished(renderThread, xIdx, yIdx);
      }
    }
  }

  long t0 = System.currentTimeMillis();

  protected void applySamplesToRaster() {
    raster.addSamples(plotBuffer, plotBufferIdx);
    plotBufferIdx = 0;
    long t1 = System.currentTimeMillis();
    if (t1 - t0 > 125) {
      try {
        Thread.sleep(1);
      }
      catch (InterruptedException e) {
        e.printStackTrace();
      }
      t0 = t1;
    }
  }

  private PlotSample[] initPlotBuffer() {
    PlotSample[] res = new PlotSample[Tools.PLOT_BUFFER_SIZE];
    for (int i = 0; i < res.length; i++) {
      res[i] = new PlotSample();
    }
    return res;
  }

  protected void transformPlotColor(XYZPoint p) {
    if (fabs(p.modGamma) > EPSILON) {
      double gamma = 4.2 / (4.2 - p.modGamma);
      double alpha = plotRed * 0.299 + plotGreen * 0.588 + plotBlue * 0.113;
      if (alpha > EPSILON) {
        double modAlpha = Math.pow(alpha, gamma);
        plotRed *= modAlpha / alpha;
        plotGreen *= modAlpha / alpha;
        plotBlue *= modAlpha / alpha;
      }
    }
    if (fabs(p.modContrast) > EPSILON) {
      double gamma = 1.2 / (1.2 - p.modContrast * 0.5);
      plotRed = (plotRed - 127.5) * gamma + 127.5;
      plotGreen = (plotGreen - 127.5) * gamma + 127.5;
      plotBlue = (plotBlue - 127.5) * gamma + 127.5;
    }
    if (fabs(p.modSaturation) > EPSILON) {
      double avg = plotRed * 0.299 + plotGreen * 0.588 + plotBlue * 0.113;
      plotRed += (plotRed - avg) * p.modSaturation;
      plotGreen += (plotGreen - avg) * p.modSaturation;
      plotBlue += (plotBlue - avg) * p.modSaturation;
    }
  }

  public interface PointProjector {
    void projectPoint(XYZPoint q);
  }

  public abstract class AxisSymmetryProjector implements PointProjector {
    protected final PointProjector parent;
    protected final double distance;
    protected final double centreX;
    protected final double centreY;
    protected final double rotation;
    protected final XYZPoint a;
    protected final XYZPoint b;
    protected final double sina, cosa;
    protected final double halve_dist;
    protected final boolean doRotate;

    public AxisSymmetryProjector(PointProjector pParent, double pDistance, double pCentreX, double pCentreY, double pRotation) {
      parent = pParent;
      distance = pDistance;
      centreX = pCentreX;
      centreY = pCentreY;
      rotation = pRotation;
      a = new XYZPoint();
      b = new XYZPoint();

      double a = rotation * M_2PI / 180.0 / 2.0;
      doRotate = fabs(a) > EPSILON;

      sina = sin(a);
      cosa = cos(a);
      halve_dist = distance / 2.0;
    }
  }

  public class XAxisSymmetryProjector extends AxisSymmetryProjector implements PointProjector {

    public XAxisSymmetryProjector(PointProjector pParent, double pDistance, double pCentreX, double pCentreY, double pRotation) {
      super(pParent, pDistance, pCentreX, pCentreY, pRotation);
    }

    @Override
    public void projectPoint(XYZPoint q) {
      if (q.doHide)
        return;
      a.assign(q);
      b.assign(q);
      double dx, dy;
      dx = q.x - centreX;
      a.x = centreX + dx + halve_dist;
      b.x = centreX - dx - halve_dist;
      if (doRotate) {
        dx = a.x - centreX;
        dy = a.y - centreY;
        a.x = centreX + dx * cosa + dy * sina;
        a.y = centreY + dy * cosa - dx * sina;

        dx = b.x - centreX;
        dy = b.y - centreY;
        b.x = centreX + dx * cosa - dy * sina;
        b.y = centreY + dy * cosa + dx * sina;
      }

      parent.projectPoint(a);
      parent.projectPoint(b);
    }

  }

  public class YAxisSymmetryProjector extends AxisSymmetryProjector implements PointProjector {

    public YAxisSymmetryProjector(PointProjector pParent, double pDistance, double pCentreX, double pCentreY, double pRotation) {
      super(pParent, pDistance, pCentreX, pCentreY, pRotation);
    }

    @Override
    public void projectPoint(XYZPoint q) {
      if (q.doHide)
        return;
      a.assign(q);
      b.assign(q);
      double dx, dy;
      dy = q.y - centreY;
      a.y = centreY + dy + halve_dist;
      b.y = centreY - dy - halve_dist;
      if (doRotate) {
        dx = a.x - centreX;
        dy = a.y - centreY;
        a.x = centreX + dx * cosa + dy * sina;
        a.y = centreY + dy * cosa - dx * sina;

        dx = b.x - centreX;
        dy = b.y - centreY;
        b.x = centreX + dx * cosa - dy * sina;
        b.y = centreY + dy * cosa + dx * sina;
      }

      parent.projectPoint(a);
      parent.projectPoint(b);
    }

  }

  public class PointSymmetryProjector implements PointProjector {
    private final PointProjector parent;
    private final double centreX, centreY;
    private final int order;
    private final double sina[], cosa[];
    private final XYZPoint ps;

    public PointSymmetryProjector(PointProjector pParent, int pOrder, double pCentreX, double pCentreY) {
      parent = pParent;
      order = pOrder >= 1 ? pOrder : 1;
      centreX = pCentreX;
      centreY = pCentreY;
      ps = new XYZPoint();
      sina = new double[order];
      cosa = new double[order];
      double da = M_2PI / (double) order;
      double angle = 0.0;
      for (int i = 0; i < order; i++) {
        sina[i] = sin(angle);
        cosa[i] = cos(angle);
        angle += da;
      }
    }

    @Override
    public void projectPoint(XYZPoint q) {
      if (q.doHide)
        return;
      double dx = q.x - centreX;
      double dy = q.y - centreY;
      parent.projectPoint(q);
      for (int i = 0; i < order; i++) {
        ps.assign(q);
        ps.x = centreX + dx * cosa[i] + dy * sina[i];
        ps.y = centreY + dy * cosa[i] - dx * sina[i];
        parent.projectPoint(ps);
      }
    }
  }

  public class DefaultPointProjector implements PointProjector {

    @Override
    public void projectPoint(XYZPoint q) {
      if (q.doHide || !view.project(q, prj) || q.isNaN())
        return;
      int xIdx, yIdx;
      if ((flame.getAntialiasAmount() > EPSILON) && (flame.getAntialiasRadius() > EPSILON) && (randGen.random() > 1.0 - flame.getAntialiasAmount())) {
        double dr = exp(flame.getAntialiasRadius() * sqrt(-log(randGen.random()))) - 1.0;
        double da = randGen.random() * 2.0 * M_PI;
        xIdx = (int) (view.bws * prj.x + dr * cos(da) + 0.5);
        if (xIdx < 0 || xIdx >= renderer.rasterWidth)
          return;
        yIdx = (int) (view.bhs * prj.y + dr * sin(da) + 0.5);
        if (yIdx < 0 || yIdx >= renderer.rasterHeight)
          return;
      }
      else {
        xIdx = (int) (view.bws * prj.x + 0.5);
        if (xIdx < 0 || xIdx >= renderer.rasterWidth)
          return;
        yIdx = (int) (view.bhs * prj.y + 0.5);
        if (yIdx < 0 || yIdx >= renderer.rasterHeight)
          return;
      }
      double intensity = prj.intensity * layer.getWeight();
      plotPoint(xIdx, yIdx, intensity);
    }

  }

  public void cleanup() {
    if (plotBufferIdx > 0) {
      applySamplesToRaster();
    }
  }
}
