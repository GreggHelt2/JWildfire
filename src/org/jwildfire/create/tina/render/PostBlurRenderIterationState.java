package org.jwildfire.create.tina.render;

import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.palette.RenderColor;
import org.jwildfire.create.tina.random.AbstractRandomGenerator;
import org.jwildfire.create.tina.variation.FlameTransformationContext;

public class PostBlurRenderIterationState extends DefaultRenderIterationState {
  private static final long serialVersionUID = 1L;
  private int blurRadius;
  private double blurFade;
  private double blurKernel[][];
  private int rasterWidth, rasterHeight;

  public PostBlurRenderIterationState(AbstractRenderThread pRenderThread, FlameRenderer pRenderer, RenderPacket pPacket, Layer pLayer, FlameTransformationContext pCtx, AbstractRandomGenerator pRandGen) {
    super(pRenderThread, pRenderer, pPacket, pLayer, pCtx, pRandGen);
  }

  @Override
  public void init() {
    blurRadius = flame.getPostBlurRadius();
    blurKernel = createBlurKernel(blurRadius, flame.getPostBlurFallOff());

    blurFade = flame.getPostBlurFade();
    if (blurFade < 0.0) {
      blurFade = 0.0;
    }
    else if (blurFade > 1.0) {
      blurFade = 1.0;
    }
    rasterWidth = renderer.rasterWidth;
    rasterHeight = renderer.rasterHeight;
  }

  private double[][] createBlurKernel(int radius, double fallOff) {
    double kernel[][];
    if (radius <= 0) {
      kernel = new double[1][1];
      kernel[0][0] = 1.0;
    }
    else {
      kernel = new double[2 * radius + 1][2 * radius + 1];
      kernel[0][0] = 1.0;
      for (int y = -radius, i = 0; y <= radius; y++, i++) {
        for (int x = -radius, j = 0; x <= radius; x++, j++) {
          double r = Math.sqrt(x * x + y * y) / (double) (radius);
          kernel[i][j] = Math.pow(2, -r * fallOff);
          //          System.out.println(i + " " + j + " " + " r=" + r + ", k=" + kernel[i][j]);
        }
      }
    }
    return kernel;
  }

  @Override
  protected void plotPoint(int screenX, int screenY, double rawX, double rawY, double intensity, XYZPoint origin) {
    if (p.rgbColor) {
      plotRed = p.redColor;
      plotGreen = p.greenColor;
      plotBlue = p.blueColor;
    }
    else {
      RenderColor color = colorProvider.getColor(p, q);
      plotRed = color.red;
      plotGreen = color.green;
      plotBlue = color.blue;
    }
    transformPlotColor(p);

    if (ctx.random() > blurFade) {
      for (int k = screenY - blurRadius, yk = 0; k <= screenY + blurRadius; k++, yk++) {
        if (k >= 0 && k < rasterHeight) {
          for (int l = screenX - blurRadius, xk = 0; l <= screenX + blurRadius; l++, xk++) {
            if (l >= 0 && l < rasterWidth) {
              // y, x

              double scl = blurKernel[yk][xk];

              double finalRed = plotRed * scl * prj.intensity;
              double finalGreen = plotGreen * scl * prj.intensity;
              double finalBlue = plotBlue * scl * prj.intensity;
              plotBuffer[plotBufferIdx++].set(l, k, finalRed, finalGreen, finalBlue, rawX, rawY, prj.z * view.bws, p.material, prj.dofDist, origin.x, origin.y, origin.z, p.receiveOnlyShadows);
              if (plotBufferIdx >= plotBuffer.length) {
                applySamplesToRaster();
              }

              if (observers != null && observers.size() > 0) {
                for (IterationObserver observer : observers) {
                  observer.notifyIterationFinished(renderThread, l, k, prj, q.x, q.y, q.z, finalRed, finalGreen, finalBlue);
                }
              }
            }
          }
        }
      }
    }
    else {
      double finalRed = plotRed * prj.intensity;
      double finalGreen = plotGreen * prj.intensity;
      double finalBlue = plotBlue * prj.intensity;
      plotBuffer[plotBufferIdx++].set(screenX, screenY, finalRed, finalGreen, finalBlue, rawX, rawY, prj.z * view.bws, p.material, prj.dofDist, origin.x, origin.y, origin.z, p.receiveOnlyShadows);
      if (plotBufferIdx >= plotBuffer.length) {
        applySamplesToRaster();
      }
      if (observers != null && observers.size() > 0) {
        for (IterationObserver observer : observers) {
          observer.notifyIterationFinished(renderThread, screenX, screenY, prj, q.x, q.y, q.z, finalRed, finalGreen, finalBlue);
        }
      }
    }
  }
}
