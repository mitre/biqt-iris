> NOTICE
> 
> This software (or technical data) was produced for the U. S. Government under contract, and is subject to the Rights in Data-General Clause 52.227-14, Alt. IV (DEC 2007)
> 
> (C) 2019 The MITRE Corporation. All Rights Reserved.
> Approved for Public Release; Distribution Unlimited. Public Release Case Number 18-0812.

## Summary ##

BIQTIris is a reference library for computing iris image quality metrics. It is
part of the open-source BIQT Framework.

Support is currently limited to mage files of size 480x480 and 640x480 pixels.
Raw image data is assumed to have a color depth of 8 bpp.

### Quality Attributes ###

The following quality attributes are reported by this provider:

  * `quality` - An overall quality ranging from 0 (low) to 100 (high)
  * `contrast` - Raw score quantifying overall image contrast.
  * `sharpness` - Raw score quantifying the "sharpness" of the image.
  * `iris_diameter` - Raw diameter of the iris measured in pixels.
  * `iris_sclera_gs` - Raw measure quantifying how distinguishable the boundary is between the iris and the sclera.
  * `iris_pupil_gs` - Raw measure quantifying how distinguishable the boundary is between the pupil and the iris.
  * `percent_visible_iris` - Percent of the iris texture that is visible within the boundary of the iris.

The following metrics are also reported by this provider. Each metric reports
a value between 0 (low) and 1 (high).

  * `normalized_contrast`
  * `normalized_sharpness`
  * `normalized_iris_diameter`
  * `normalized_iris_pupil_gs`
  * `normalized_iris_sclera_gs`  
  * `normalized_percent_visible_iris`

### Features ###

The provider reports the following features:

  * `image_width` - Image width measured in pixels.
  * `image_height` - Image height measured in pixels.
  * `iris_center_x` - The x-coordinate of the iris center in the image.
  * `iris_center_y` - The y-coordinate of the iris center in the image.
