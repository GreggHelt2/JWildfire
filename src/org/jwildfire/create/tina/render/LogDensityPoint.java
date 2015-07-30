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

import org.jwildfire.create.tina.base.raster.RasterPoint;

public class LogDensityPoint {
  public double red;
  public double green;
  public double blue;
  public double intensity;

  public final RasterPoint rp = new RasterPoint();
  public final RasterPoint lu = new RasterPoint();
  public final RasterPoint ru = new RasterPoint();
  public final RasterPoint lb = new RasterPoint();
  public final RasterPoint rb = new RasterPoint();

  public void clear() {
    red = green = blue = intensity = 0.0;
  }

  public void clip() {
    if (red < 0.0) {
      red = 0.0;
    }
    if (green < 0.0) {
      green = 0.0;
    }
    if (blue < 0.0) {
      blue = 0.0;
    }
    if (intensity < 0.0) {
      intensity = 0.0;
    }
  }

}
