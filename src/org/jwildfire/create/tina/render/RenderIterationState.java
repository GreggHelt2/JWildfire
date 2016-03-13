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

import java.io.Serializable;
import java.util.List;

import org.jwildfire.create.tina.base.Flame;
import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.raster.AbstractRaster;
import org.jwildfire.create.tina.io.MapGradientReader;
import org.jwildfire.create.tina.palette.RGBPalette;
import org.jwildfire.create.tina.palette.RenderColor;
import org.jwildfire.create.tina.random.AbstractRandomGenerator;
import org.jwildfire.create.tina.variation.FlameTransformationContext;

public class RenderIterationState implements Serializable {
  private static final long serialVersionUID = 2L;
  public static RGBPalette density_default_palette;
  // public RGBPalette density_default_palette;
  private boolean density_palette_missing = false;
  
  protected final AbstractRenderThread renderThread;
  protected final FlameRenderer renderer;
  protected final FlameRendererView view;
  protected final RenderPacket packet;
  protected final Flame flame;
  protected final Layer layer;
  protected AbstractRaster raster;
  protected final FlameTransformationContext ctx;
  protected final AbstractRandomGenerator randGen;
  protected final List<IterationObserver> observers;
  protected final RenderColor[] colorMap;
  protected final double paletteIdxScl;

  public RenderIterationState(AbstractRenderThread pRenderThread, FlameRenderer pRenderer, RenderPacket pPacket, Layer pLayer, FlameTransformationContext pCtx, AbstractRandomGenerator pRandGen) {
    renderThread = pRenderThread;
    renderer = pRenderer;
    packet = pPacket;
    view = pPacket.getView();
    flame = pPacket.getFlame();
    layer = pLayer;
    ctx = pCtx;
    randGen = pRandGen;
    observers = renderer.getIterationObservers();

    if (pRenderer.prefs.isDensityPostProcess()) {
      if (density_default_palette == null && !density_palette_missing) {
        createDensityDefaultPalette();
      }
      if (density_palette_missing) {
        colorMap = pLayer.getPalette().createRenderPalette(flame.getWhiteLevel());
      }
      else {
        colorMap = density_default_palette.createRenderPalette(flame.getWhiteLevel());
      }
    }
    else {
      colorMap = pLayer.getPalette().createRenderPalette(flame.getWhiteLevel());
    }
    paletteIdxScl = colorMap.length - 2;
    raster = pRenderer.getRaster();
  }
  
 //  public static void createDensityDefaultPalette() {
  public void createDensityDefaultPalette() {
    // String palette_file = "/Users/gregg";
    String palette_file = renderer.prefs.getDefaultDensityColormap();
    System.out.println("getting default density colormap: " + palette_file);
    MapGradientReader greader = new MapGradientReader();
    try {
      List<RGBPalette> palettes = greader.readPalettes(palette_file);
      density_default_palette = palettes.get(0);
    }
    catch (Exception ex) {
      ex.printStackTrace();
      density_default_palette = null;
      density_palette_missing = true;
    }
    System.out.println("default density palette: " + density_default_palette);

  }

}
