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
package org.jwildfire.create.tina.script.ui;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.DefaultListModel;
import javax.swing.ImageIcon;

import javax.swing.JPanel;
import javax.swing.JList;
import javax.swing.JScrollPane;

public class ImageListBuilder extends FieldBuilder {
  private String initialValue = "";
  private String imageDir;

  protected ImageListBuilder(ContainerBuilder pParent) {
    super(pParent);
  }

  @Override
  public ImageListBuilder withCaption(String pCaption) {
    super.withCaption(pCaption);
    return this;
  }

  @Override
  public ImageListBuilder withPropertyName(String pPropertyName) {
    super.withPropertyName(pPropertyName);
    return this;
  }

  public ImageListBuilder withImageDirectory(String pImageDir) {
    imageDir = pImageDir;
    return this;
  }

  public ImageListBuilder withInitialValue(String pInitialValue) {
    initialValue = pInitialValue;
    return this;
  }

  public ContainerBuilder closeImageList() {
    return parent;
  }

  public JList buildPart(ScriptParamsForm pForm, JPanel pPanel, int xOff, int yOff) {
    createLabel(pPanel, xOff, yOff);
    JList imageList = new JList();
    int xydim = 500;
    imageList.setName(propertyName);
    imageList.removeAll();

    DefaultListModel listModel = new DefaultListModel();
    File folder = new File(imageDir);
    File[] listOfFiles = folder.listFiles();
    for (int i = 0; i < listOfFiles.length; i++) {
      File fil = listOfFiles[i];
      String name = fil.getName();
      // load only JPEGs and PNGs
      if ( name.endsWith(".jpg") || name.endsWith(".png")) {
        ImageIcon imicon;
        try {
          imicon = new ImageIcon(ImageIO.read(fil));
          name = name.substring(0, name.lastIndexOf("."));          
          imicon.setDescription(name); // probably want to strip off suffix?
          listModel.add(i, imicon);
        } catch (IOException ex) {
          Logger.getLogger(ImageListBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }

      }
    }

    imageList.setModel(listModel);
    if (initialValue != null) {
      imageList.setSelectedValue(initialValue, true);
      //imageList.setSelectedItem(initialValue);
    }
    else {
      imageList.setSelectedIndex(0);
    }
    imageList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
    imageList.setVisibleRowCount(-1);

    JScrollPane listScroller = new JScrollPane(imageList);
    listScroller.setPreferredSize(new Dimension(xydim, xydim));
    listScroller.setSize(xydim, xydim);
    listScroller.setMinimumSize(new Dimension(xydim, xydim));
    listScroller.setLocation(new Point(xOff + LABEL_WIDTH + H_BORDER, yOff));
    listScroller.setFont(new Font("Dialog", Font.PLAIN, 10));
    listScroller.setBounds(xOff + LABEL_WIDTH + H_BORDER, yOff, xydim, xydim);
    pPanel.add(listScroller);
    // pPanel.add(imageList);
    return imageList;
  }

}
