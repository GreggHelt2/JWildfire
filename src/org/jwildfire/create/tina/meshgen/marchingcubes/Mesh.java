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
package org.jwildfire.create.tina.meshgen.marchingcubes;

import java.util.Collections;
import java.util.List;

public class Mesh {
  private List<Point3f> vertices;
  private List<Point3f> vertexNormals;
  private List<Face> faces;
  private List<Point2f> textureCoords;
  private boolean dirty;
  private Point3f[] bounds = null;

  public Mesh() {
    setVertices(null);
    setVertexNormals(null);
    setFaces(null);
    setTextureCoords(null);
    dirty = true;
  }

  public Mesh(List<Point3f> pVertices, List<Face> pFaces) {
    setVertices(pVertices);
    setVertexNormals(null);
    setFaces(pFaces);
    setTextureCoords(null);
    dirty = true;
  }

  public Mesh(List<Point3f> pVertices, List<Point3f> pVertexNormals, List<Face> pFaces) {
    setVertices(pVertices);
    setVertexNormals(pVertexNormals);
    setFaces(pFaces);
    setTextureCoords(null);
    dirty = true;
  }

  public List<Point3f> getVertices() {
    return vertices;
  }

  public void setVertices(List<Point3f> pVertices) {
    if (pVertices != null) {
      vertices = pVertices;
    }
    else {
      vertices = Collections.emptyList();
    }
    dirty = true;
  }

  public List<Face> getFaces() {
    return faces;
  }

  public void setFaces(List<Face> pFaces) {
    if (pFaces != null) {
      faces = pFaces;
    }
    else {
      faces = Collections.emptyList();
    }
    dirty = true;
  }

  private Point3f[] getBounds() {
    if (dirty) {
      Point3f min, max;
      bounds = new Point3f[2];
      bounds[0] = min = new Point3f();
      bounds[1] = max = new Point3f();
      if (vertices.size() > 0) {
        Point3f first = vertices.get(0);
        bounds[0].x = bounds[1].x = first.x;
        bounds[0].y = bounds[1].y = first.y;
        bounds[0].z = bounds[1].z = first.z;
        for (int i = 1; i < vertices.size(); i++) {
          Point3f vertex = vertices.get(i);
          if (vertex.x < min.x)
            min.x = vertex.x;
          else if (vertex.x > max.x)
            max.x = vertex.x;
          if (vertex.y < min.y)
            min.y = vertex.y;
          else if (vertex.y > max.y)
            max.y = vertex.y;
          if (vertex.z < min.z)
            min.z = vertex.z;
          else if (vertex.z > max.z)
            max.z = vertex.z;
        }
      }
      dirty = false;
    }
    return bounds;
  }

  public Point3f getPMin() {
    return getBounds()[0];
  }

  public Point3f getPMax() {
    return getBounds()[1];
  }

  public List<Point3f> getVertexNormals() {
    return vertexNormals;
  }

  public void setVertexNormals(List<Point3f> pVertexNormals) {
    if (pVertexNormals != null) {
      vertexNormals = pVertexNormals;
    }
    else {
      vertexNormals = Collections.emptyList();
    }
  }

  public List<Point2f> getTextureCoords() {
    return textureCoords;
  }

  public void setTextureCoords(List<Point2f> pTextureCoords) {
    if (pTextureCoords != null) {
      textureCoords = pTextureCoords;
    }
    else {
      textureCoords = Collections.emptyList();
    }
    dirty = true;
  }
}
