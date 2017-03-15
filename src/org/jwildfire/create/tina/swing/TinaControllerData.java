/*
  JWildfire - an image and animation processor written in Java 
  Copyright (C) 1995-2016 Andreas Maschke

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
package org.jwildfire.create.tina.swing;

import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.JToggleButton;
import javax.swing.JTree;

import org.jwildfire.create.tina.palette.RGBColor;
import org.jwildfire.swing.ImagePanel;

public class TinaControllerData {
  public JComboBox qualityProfileCmb;
  public JComboBox resolutionProfileCmb;
  public JComboBox gpuQualityProfileCmb;
  public JComboBox gpuResolutionProfileCmb;
  public JComboBox batchQualityProfileCmb;
  public JComboBox batchResolutionProfileCmb;
  public JComboBox interactiveResolutionProfileCmb;
  public JComboBox swfAnimatorResolutionProfileCmb;
  public JComboBox swfAnimatorQualityProfileCmb;
  public JWFNumberField cameraRollREd;
  public JSlider cameraRollSlider;
  public JWFNumberField cameraPitchREd;
  public JSlider cameraPitchSlider;
  public JWFNumberField cameraYawREd;
  public JSlider cameraYawSlider;
  public JWFNumberField cameraPerspectiveREd;
  public JSlider cameraPerspectiveSlider;
  public JWFNumberField camPosXREd;
  public JSlider camPosXSlider;
  public JWFNumberField camPosYREd;
  public JSlider camPosYSlider;
  public JWFNumberField camPosZREd;
  public JSlider camPosZSlider;
  public JWFNumberField cameraCentreXREd;
  public JSlider cameraCentreXSlider;
  public JWFNumberField cameraCentreYREd;
  public JSlider cameraCentreYSlider;
  public JWFNumberField cameraZoomREd;
  public JSlider cameraZoomSlider;
  public JWFNumberField focusXREd;
  public JSlider focusXSlider;
  public JWFNumberField focusYREd;
  public JSlider focusYSlider;
  public JWFNumberField focusZREd;
  public JSlider focusZSlider;
  public JWFNumberField dimishZREd;
  public JSlider dimishZSlider;
  public JWFNumberField camZREd;
  public JSlider camZSlider;
  public JWFNumberField cameraDOFREd;
  public JSlider cameraDOFSlider;
  public JWFNumberField cameraDOFAreaREd;
  public JSlider cameraDOFAreaSlider;
  public JWFNumberField cameraDOFExponentREd;
  public JSlider cameraDOFExponentSlider;
  public JCheckBox newDOFCBx;
  public JWFNumberField pixelsPerUnitREd;
  public JSlider pixelsPerUnitSlider;
  public JWFNumberField brightnessREd;
  public JSlider brightnessSlider;
  public JWFNumberField contrastREd;
  public JSlider contrastSlider;
  public JWFNumberField whiteLevelREd;
  public JSlider whiteLevelSlider;
  public JWFNumberField gammaREd;
  public JSlider gammaSlider;
  public JWFNumberField vibrancyREd;
  public JSlider vibrancySlider;
  public JWFNumberField saturationREd;
  public JSlider saturationSlider;
  public JWFNumberField filterRadiusREd;
  public JSlider filterRadiusSlider;
  public JComboBox filterKernelCmb;
  public JWFNumberField gammaThresholdREd;
  public JSlider gammaThresholdSlider;
  public JCheckBox bgTransparencyCBx;
  public JWFNumberField postBlurRadiusREd;
  public JSlider postBlurRadiusSlider;
  public JWFNumberField postBlurFadeREd;
  public JSlider postBlurFadeSlider;
  public JWFNumberField postBlurFallOffREd;
  public JSlider postBlurFallOffSlider;
  public JWFNumberField tinaZBufferScaleREd;
  public JSlider tinaZBufferScaleSlider;
  public JTextField paletteRandomPointsREd;
  public JComboBox paletteRandomGeneratorCmb;
  public JCheckBox paletteFadeColorsCBx;
  public JPanel paletteImgPanel;
  public JPanel colorChooserPaletteImgPanel;
  public ImagePanel palettePanel;
  public ImagePanel colorChooserPalettePanel;
  public ImagePanel filterKernelPreviewPanel;
  public JWFNumberField paletteShiftREd;
  public JSlider paletteShiftSlider;
  public JWFNumberField paletteRedREd;
  public JSlider paletteRedSlider;
  public JWFNumberField paletteGreenREd;
  public JSlider paletteGreenSlider;
  public JWFNumberField paletteBlueREd;
  public JSlider paletteBlueSlider;
  public JWFNumberField paletteHueREd;
  public JSlider paletteHueSlider;
  public JWFNumberField paletteSaturationREd;
  public JSlider paletteSaturationSlider;
  public JWFNumberField paletteContrastREd;
  public JSlider paletteContrastSlider;
  public JWFNumberField paletteGammaREd;
  public JSlider paletteGammaSlider;
  public JWFNumberField paletteBrightnessREd;
  public JSlider paletteBrightnessSlider;
  public JWFNumberField paletteSwapRGBREd;
  public JSlider paletteSwapRGBSlider;
  public JWFNumberField paletteFrequencyREd;
  public JSlider paletteFrequencySlider;
  public JWFNumberField paletteBlurREd;
  public JSlider paletteBlurSlider;
  public JButton paletteInvertBtn;
  public JButton paletteReverseBtn;
  public JTable renderBatchJobsTable;
  public JProgressBar batchRenderJobProgressBar;
  public JProgressBar batchRenderTotalProgressBar;
  public JPanel batchPreviewRootPanel;
  public JToggleButton affineScaleXButton;
  public JToggleButton affineScaleYButton;
  public JTable transformationsTable;
  public JButton affineResetTransformButton;
  public JWFNumberField affineC00REd;
  public JWFNumberField affineC01REd;
  public JWFNumberField affineC10REd;
  public JWFNumberField affineC11REd;
  public JWFNumberField affineC20REd;
  public JWFNumberField affineC21REd;
  public JWFNumberField affineRotateAmountREd;
  public JWFNumberField affineScaleAmountREd;
  public JWFNumberField affineMoveHorizAmountREd;
  public JWFNumberField affineMoveVertAmountREd;
  public JButton affineRotateLeftButton;
  public JButton affineRotateRightButton;
  public JButton affineEnlargeButton;
  public JButton affineShrinkButton;
  public JButton affineMoveUpButton;
  public JButton affineMoveLeftButton;
  public JButton affineMoveRightButton;
  public JButton affineMoveDownButton;
  public JButton affineFlipHorizontalButton;
  public JButton affineFlipVerticalButton;
  public JButton addTransformationButton;
  public JButton addLinkedTransformationButton;
  public JButton duplicateTransformationButton;
  public JButton deleteTransformationButton;
  public JButton addFinalTransformationButton;
  public JWFNumberField transformationWeightREd;
  public JButton affineRotateEditMotionCurveBtn;
  public JButton affineScaleEditMotionCurveBtn;
  public JToggleButton affineEditPostTransformButton;
  public JToggleButton affineEditPostTransformSmallButton;
  public JToggleButton affinePreserveZButton;
  public JToggleButton affineMirrorPrePostTranslationsButton;
  public JButton randomizeButton;
  public JToggleButton toggleVariationsButton;
  public JToggleButton toggleTransparencyButton;
  public JPanel gradientLibraryPanel;
  public JPanel randomBatchPanel;
  public JScrollPane randomBatchScrollPane;
  public TinaNonlinearControlsRow[] TinaNonlinearControlsRows;
  public VariationControlsDelegate[] variationControlsDelegates;
  public JWFNumberField xFormColorREd;
  public JSlider xFormColorSlider;
  public JWFNumberField xFormSymmetryREd;
  public JSlider xFormSymmetrySlider;
  public JWFNumberField xFormMaterialREd;
  public JSlider xFormMaterialSlider;
  public JWFNumberField xFormMaterialSpeedREd;
  public JSlider xFormMaterialSpeedSlider;
  public JWFNumberField xFormModGammaREd;
  public JSlider xFormModGammaSlider;
  public JWFNumberField xFormModGammaSpeedREd;
  public JSlider xFormModGammaSpeedSlider;
  public JWFNumberField xFormModContrastREd;
  public JSlider xFormModContrastSlider;
  public JWFNumberField xFormModContrastSpeedREd;
  public JSlider xFormModContrastSpeedSlider;
  public JWFNumberField xFormModSaturationREd;
  public JSlider xFormModSaturationSlider;
  public JWFNumberField xFormModSaturationSpeedREd;
  public JSlider xFormModSaturationSpeedSlider;
  public JWFNumberField xFormOpacityREd;
  public JSlider xFormOpacitySlider;
  public JComboBox xFormDrawModeCmb;
  public JWFNumberField xFormAntialiasAmountREd;
  public JSlider xFormAntialiasAmountSlider;
  public JWFNumberField xFormAntialiasRadiusREd;
  public JSlider xFormAntialiasRadiusSlider;
  public JTable relWeightsTable;
  public JButton relWeightsZeroButton;
  public JButton relWeightsOneButton;
  public JWFNumberField relWeightREd;
  public JTable createPaletteColorsTable;
  public List<RGBColor> paletteKeyFrames;
  public JButton renderFlameButton;
  public JButton renderMainButton;
  public JButton appendToMovieButton;
  public JToggleButton mouseTransformMoveTrianglesButton;
  public JToggleButton mouseTransformRotateTrianglesButton;
  public JToggleButton mouseTransformScaleTrianglesButton;
  public JToggleButton mouseTransformEditFocusPointButton;
  public JToggleButton mouseTransformEditPointsButton;
  public JToggleButton mouseTransformEditGradientButton;
  public JToggleButton mouseTransformEditTriangleViewButton;
  public JToggleButton mouseTransformEditViewButton;
  public JToggleButton mouseTransformSlowButton;
  public JToggleButton toggleTriangleWithColorsButton;
  public JToggleButton realtimePreviewToggleButton;
  public JButton batchRenderAddFilesButton;
  public JButton batchRenderFilesMoveDownButton;
  public JButton batchRenderFilesMoveUpButton;
  public JButton batchRenderFilesRemoveButton;
  public JButton batchRenderFilesRemoveAllButton;
  public JButton batchRenderStartButton;
  public JTextPane helpPane;
  public JTextPane apophysisHintsPane;
  public JButton undoButton;
  public JButton redoButton;
  public JButton editTransformCaptionButton;
  public JButton editFlameTileButton;
  public JButton snapShotButton;
  public JButton qSaveButton;
  public JButton sendToIRButton;
  public JButton bokehButton;
  public JButton movieButton;
  public JToggleButton transformSlowButton;
  public JToggleButton transparencyButton;
  public JTree scriptTree;
  public JTextArea scriptDescriptionTextArea;
  public JTextArea scriptTextArea;
  public JButton rescanScriptsBtn;
  public JButton newScriptBtn;
  public JButton newScriptFromFlameBtn;
  public JButton deleteScriptBtn;
  public JButton scriptRenameBtn;
  public JButton scriptDuplicateBtn;
  public JButton scriptRunBtn;
  public JButton scriptEditBtn;
  public JTree gradientLibTree;
  public JButton backgroundColorIndicatorBtn;
  public JWFNumberField layerWeightEd;
  public JButton layerAddBtn;
  public JButton layerDuplicateBtn;
  public JButton layerDeleteBtn;
  public JTable layersTable;
  public JToggleButton layerVisibleBtn;
  public JToggleButton layerAppendBtn;
  public JToggleButton layerPreviewBtn;
  public JButton layerHideOthersBtn;
  public JButton layerShowAllBtn;
  public JWFNumberField motionBlurLengthField;
  public JSlider motionBlurLengthSlider;
  public JWFNumberField motionBlurTimeStepField;
  public JSlider motionBlurTimeStepSlider;
  public JWFNumberField motionBlurDecayField;
  public JSlider motionBlurDecaySlider;
  public JComboBox postSymmetryTypeCmb;
  public JWFNumberField postSymmetryDistanceREd;
  public JSlider postSymmetryDistanceSlider;
  public JWFNumberField postSymmetryRotationREd;
  public JSlider postSymmetryRotationSlider;
  public JWFNumberField postSymmetryOrderREd;
  public JSlider postSymmetryOrderSlider;
  public JWFNumberField postSymmetryCentreXREd;
  public JSlider postSymmetryCentreXSlider;
  public JWFNumberField postSymmetryCentreYREd;
  public JSlider postSymmetryCentreYSlider;
  public JComboBox stereo3dModeCmb;
  public JWFNumberField stereo3dAngleREd;
  public JSlider stereo3dAngleSlider;
  public JWFNumberField stereo3dEyeDistREd;
  public JSlider stereo3dEyeDistSlider;
  public JComboBox stereo3dLeftEyeColorCmb;
  public JComboBox stereo3dRightEyeColorCmb;
  public JWFNumberField stereo3dInterpolatedImageCountREd;
  public JSlider stereo3dInterpolatedImageCountSlider;
  public JComboBox stereo3dPreviewCmb;
  public JWFNumberField stereo3dFocalOffsetREd;
  public JSlider stereo3dFocalOffsetSlider;
  public JCheckBox stereo3dSwapSidesCBx;
  public JToggleButton toggleDrawGridButton;
  public JToggleButton toggleDrawGuidesButton;
  public JComboBox triangleStyleCmb;
  public JButton channelMixerResetBtn;
  public JComboBox channelMixerModeCmb;
  public JPanel channelMixerRRRootPanel;
  public JPanel channelMixerRGRootPanel;
  public JPanel channelMixerRBRootPanel;
  public JPanel channelMixerGRRootPanel;
  public JPanel channelMixerGGRootPanel;
  public JPanel channelMixerGBRootPanel;
  public JPanel channelMixerBRRootPanel;
  public JPanel channelMixerBGRootPanel;
  public JPanel channelMixerBBRootPanel;
  public JComboBox dofDOFShapeCmb;
  public JWFNumberField dofDOFScaleREd;
  public JSlider dofDOFScaleSlider;
  public JWFNumberField dofDOFAngleREd;
  public JSlider dofDOFAngleSlider;
  public JWFNumberField dofDOFFadeREd;
  public JSlider dofDOFFadeSlider;
  public JWFNumberField dofDOFParam1REd;
  public JSlider dofDOFParam1Slider;
  public JLabel dofDOFParam1Lbl;
  public JWFNumberField dofDOFParam2REd;
  public JSlider dofDOFParam2Slider;
  public JLabel dofDOFParam2Lbl;
  public JWFNumberField dofDOFParam3REd;
  public JSlider dofDOFParam3Slider;
  public JLabel dofDOFParam3Lbl;
  public JWFNumberField dofDOFParam4REd;
  public JSlider dofDOFParam4Slider;
  public JLabel dofDOFParam4Lbl;
  public JWFNumberField dofDOFParam5REd;
  public JSlider dofDOFParam5Slider;
  public JLabel dofDOFParam5Lbl;
  public JWFNumberField dofDOFParam6REd;
  public JSlider dofDOFParam6Slider;
  public JLabel dofDOFParam6Lbl;
  public JButton resetCameraSettingsBtn;
  public JButton resetDOFSettingsButton;
  public JButton resetBokehOptionsButton;
  public JButton resetColoringOptionsButton;
  public JButton resetAntialiasOptionsButton;
  public JButton resetShadingSettingsBtn;
  public JButton resetStereo3DSettingsBtn;
  public JButton resetPostSymmetrySettingsBtn;
  public JButton resetMotionBlurSettingsBtn;
  public JRadioButton xaosViewAsToBtn;
  public JRadioButton xaosViewAsFromBtn;
  public JPanel previewEastMainPanel;
  public JPanel macroButtonVertPanel;
  public JPanel macroButtonHorizPanel;
  public JTable macroButtonsTable;
  public JButton macroButtonMoveUpBtn;
  public JButton macroButtonMoveDownBtn;
  public JButton macroButtonDeleteBtn;
  public JToggleButton toggleDetachedPreviewButton;
  public JButton gradientResetBtn;
  public JPanel macroButtonHorizRootPanel;
  public JToggleButton affineXYEditPlaneToggleBtn;
  public JToggleButton affineYZEditPlaneToggleBtn;
  public JToggleButton affineZXEditPlaneToggleBtn;
  public JWFNumberField gradientColorMapHorizOffsetREd;
  public JSlider gradientColorMapHorizOffsetSlider;
  public JWFNumberField gradientColorMapHorizScaleREd;
  public JSlider gradientColorMapHorizScaleSlider;
  public JWFNumberField gradientColorMapVertOffsetREd;
  public JSlider gradientColorMapVertOffsetSlider;
  public JWFNumberField gradientColorMapVertScaleREd;
  public JSlider gradientColorMapVertScaleSlider;
  public JWFNumberField gradientColorMapLocalColorAddREd;
  public JSlider gradientColorMapLocalColorAddSlider;
  public JWFNumberField gradientColorMapLocalColorScaleREd;
  public JSlider gradientColorMapLocalColorScaleSlider;
  public JWFNumberField flameFPSField;
  public JPanel filterKernelPreviewRootPnl;
  public JWFNumberField tinaSpatialOversamplingREd;
  public JSlider tinaSpatialOversamplingSlider;
  public JWFNumberField tinaColorOversamplingREd;
  public JSlider tinaColorOversamplingSlider;
  public JCheckBox tinaSampleJitteringCheckBox;
  public JToggleButton filterKernelFlatPreviewBtn;
  public JCheckBox tinaPostNoiseFilterCheckBox;
  public JWFNumberField tinaPostNoiseThresholdField;
  public JSlider tinaPostNoiseThresholdSlider;
  public JWFNumberField foregroundOpacityField;
  public JSlider foregroundOpacitySlider;

  public JToggleButton solidRenderingToggleBtn;
  public JCheckBox tinaSolidRenderingEnableAOCBx;
  public JWFNumberField tinaSolidRenderingAOIntensityREd;
  public JSlider tinaSolidRenderingAOIntensitySlider;
  public JWFNumberField tinaSolidRenderingAOSearchRadiusREd;
  public JSlider tinaSolidRenderingAOSearchRadiusSlider;
  public JWFNumberField tinaSolidRenderingAOBlurRadiusREd;
  public JSlider tinaSolidRenderingAOBlurRadiusSlider;
  public JWFNumberField tinaSolidRenderingAOFalloffREd;
  public JSlider tinaSolidRenderingAOFalloffSlider;
  public JWFNumberField tinaSolidRenderingAORadiusSamplesREd;
  public JSlider tinaSolidRenderingAORadiusSamplesSlider;
  public JWFNumberField tinaSolidRenderingAOAzimuthSamplesREd;
  public JSlider tinaSolidRenderingAOAzimuthSamplesSlider;
  public JWFNumberField tinaSolidRenderingAOAffectDiffuseREd;
  public JSlider tinaSolidRenderingAOAffectDiffuseSlider;
  public JComboBox tinaSolidRenderingShadowTypeCmb;
  public JComboBox tinaSolidRenderingShadowmapSizeCmb;
  public JWFNumberField tinaSolidRenderingShadowSmoothRadiusREd;
  public JSlider tinaSolidRenderingShadowSmoothRadiusSlider;
  public JWFNumberField tinaSolidRenderingShadowmapBiasREd;
  public JSlider tinaSolidRenderingShadowmapBiasSlider;
  public JButton resetSolidRenderingMaterialsBtn;
  public JButton resetSolidRenderingLightsBtn;
  public JButton resetSolidRenderingHardShadowOptionsBtn;
  public JButton resetSolidRenderingAmbientShadowOptionsBtn;
  public JComboBox tinaSolidRenderingSelectedLightCmb;
  public JButton tinaSolidRenderingAddLightBtn;
  public JButton tinaSolidRenderingDeleteLightBtn;
  public JWFNumberField tinaSolidRenderingLightAltitudeREd;
  public JWFNumberField tinaSolidRenderingLightAzimuthREd;
  public JSlider tinaSolidRenderingLightAltitudeSlider;
  public JSlider tinaSolidRenderingLightAzimuthSlider;
  public JButton tinaSolidRenderingLightColorBtn;
  public JCheckBox tinaSolidRenderingLightCastShadowsCBx;
  public JWFNumberField tinaSolidRenderingLightIntensityREd;
  public JSlider tinaSolidRenderingLightIntensitySlider;
  public JWFNumberField tinaSolidRenderingShadowIntensityREd;
  public JSlider tinaSolidRenderingShadowIntensitySlider;
  public JComboBox tinaSolidRenderingSelectedMaterialCmb;
  public JButton tinaSolidRenderingAddMaterialBtn;
  public JButton tinaSolidRenderingDeleteMaterialBtn;
  public JWFNumberField tinaSolidRenderingMaterialDiffuseREd;
  public JSlider tinaSolidRenderingMaterialDiffuseSlider;
  public JWFNumberField tinaSolidRenderingMaterialAmbientREd;
  public JSlider tinaSolidRenderingMaterialAmbientSlider;
  public JWFNumberField tinaSolidRenderingMaterialSpecularREd;
  public JSlider tinaSolidRenderingMaterialSpecularSlider;
  public JWFNumberField tinaSolidRenderingMaterialSpecularSharpnessREd;
  public JSlider tinaSolidRenderingMaterialSpecularSharpnessSlider;
  public JButton tinaSolidRenderingMaterialSpecularColorBtn;
  public JComboBox tinaSolidRenderingMaterialDiffuseResponseCmb;
  public JComboBox tinaSolidRenderingMaterialReflectionMappingCmb;
  public JWFNumberField tinaSolidRenderingMaterialReflectionMapIntensityREd;
  public JSlider tinaSolidRenderingMaterialReflectionMapIntensitySlider;
  public JButton tinaSolidRenderingMaterialReflMapBtn;
  public JButton tinaSolidRenderingMaterialSelectReflMapBtn;
  public JButton tinaSolidRenderingMaterialRemoveReflMapBtn;
  public JWFNumberField xFormModHueREd;
  public JSlider xFormModHueSlider;
  public JWFNumberField xFormModHueSpeedREd;
  public JSlider xFormModHueSpeedSlider;
  public JPanel bokehSettingsPnl;
  public JPanel postBokehSettingsPnl;
  public JButton resetPostBokehSettingsBtn;
  public JWFNumberField postBokehIntensityREd;
  public JSlider postBokehIntensitySlider;
  public JWFNumberField postBokehBrightnessREd;
  public JSlider postBokehBrightnessSlider;
  public JWFNumberField postBokehSizeREd;
  public JSlider postBokehSizeSlider;
  public JWFNumberField postBokehActivationREd;
  public JSlider postBokehActivationSlider;
  public JComboBox postBokehFilterKernelCmb;

}