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

import org.jwildfire.create.tina.base.Flame;
import org.jwildfire.create.tina.random.AbstractRandomGenerator;
import org.jwildfire.envelope.Envelope;

public class IntensityColorFunc implements ColorFunc {

  /*  public IntensityColorFunc() {
    super();
    System.out.println("called IntensityColorFunc() constructor");
  }
  */
  
  private Envelope envelope;

  public double mapIntensity(double intensity) {
    // return envelope.evaluate(intensity) / intensity;
    return envelope.evaluate(intensity);
    /* double new_intensity;
    if (intensity > 3) { 
      // System.out.println("intensity > 3");
      new_intensity = 0.2;
    }
    else if (intensity < 0.5) { 
      // System.out.println("intensity < 0.5");
      new_intensity = 5.0; 
    }
    else { new_intensity = intensity; }
    return new_intensity;
            */
  }
  
  @Override
  public double mapRGBToR(double pR, double pG, double pB) {
    // double brightness = (pR + pG + pB)/3.0;
    //    double brightness = 0.2990 * pR + 0.5880 * pG + 0.1130 * pB;
    // double brightness = pR;
    // return envelope.evaluate(brightness) / brightness * pR;
    return pR;
  }

  @Override
  public double mapRGBToG(double pR, double pG, double pB) {
    // double brightness = (pR + pG + pB)/3.0;
    // double brightness = 0.2990 * pR + 0.5880 * pG + 0.1130 * pB;
    // double brightness = pG;
    // return envelope.evaluate(brightness) / brightness * pG;
    return pG;
  }

  @Override
  public double mapRGBToB(double pR, double pG, double pB) {
    // double brightness = (pR + pG + pB)/3.0;
    // double brightness = 0.2990 * pR + 0.5880 * pG + 0.1130 * pB;
    //    double brightness = pB;
    //return envelope.evaluate(brightness) / brightness * pB;
    return pB;
  }

  @Override
  public void prepare(Flame pFlame, AbstractRandomGenerator pRandGen) {
    envelope = pFlame.getMixerRRCurve().toEnvelope();
    envelope.setUseBisection(true);
  }

}
