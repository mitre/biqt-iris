> NOTICE
> 
> This software (or technical data) was produced for the U. S. Government under contract, and is subject to the Rights in Data-General Clause 52.227-14, Alt. IV (DEC 2007)
> (C) 2023 The MITRE Corporation. All Rights Reserved.
> Approved for Public Release; Distribution Unlimited. Public Release Case Number 18-0812.

## Summary ##

BIQTIris is a reference library for computing iris image statistics using various quality attributes, features, and ISO metrics. It is
part of the open-source BIQT Framework.

Images smaller than 256x256 or larger than 1000x680 are not supported. Raw image data is assumed to have a color depth of 8 bpp.

### Features ###

The following features are reported by this provider:

* `image_height` _value range: [256, +inf)_ - Image height measured in pixels.
* `image_width` _value range: [256, +inf)_ - Image width measured in pixels.
* `iris_center_x` _value range: [1, +inf)_ - The x-coordinate of the iris center in the image.
* `iris_center_y` _value range: [1, +inf)_ - The y-coordinate of the iris center in the image.
* `iris_diameter` _value range: [224, 432]_ - Raw diameter of the iris measured in pixels.
* `pupil_center_x` _value range: [1, +inf)_ - The x-coordinate of the pupil center in the image.
* `pupil_center_y` _value range: [1, +inf)_ - The y-coordinate of the pupil center in the image.
* `pupil_diameter` _value range: [134, 202]_ - Raw diameter of the pupil measured in pixels.
* `pupil_radius` _value range: [32, 101]_ - Raw radius of the pupil measured in pixels.

### Quality Metrics ###

The following quality metrics are reported by this provider:

  * `contrast` _value range: [0, +inf)_ - Raw score quantifying overall image contrast.
  * `iris_pupil_gs` _value range: [0, +inf)_ - Raw measure quantifying how distinguishable the boundary is between the pupil and the iris.
  * `iris_sclera_gs` _value range: [0, +inf)_ - Raw measure quantifying how distinguishable the boundary is between the iris and the sclera.
  * `pupil_circularity_avg_deviation` _value range: [0, +inf)_ - Average deviation for the dimensionless pupil radii at corresponding angle segments.
  * `quality` _value range: [0, 100]_ - An overall quality score that leverages several statistics together.
  * `sharpness` _value range: [0, +inf)_ - Raw score quantifying the sharpness of the image.

NEW! The following ISO metrics are reported by this provider (please refer to ISO/IEC 29794-6:2015 [Information technology — Biometric sample quality — Part 6: Iris image data] documentation for more information):

  * `iso_overall_quality` _value range: [0, 100)_ - The overall ISO quality score based on the product of normalized individual iso metrics.
  * `iso_greyscale_utilization` _value range: [0, +inf), recommended value: 6 or greater_ - The spread of intensity values regarding the pixel values within the iris portion of the image.
  * `iso_iris_pupil_concentricity` _value range: [0, 100], recommended value: 90 or greater_ - The degree to which the pupil centre and the iris centre are in the same location.
  * `iso_iris_pupil_contrast` _value range: [0, 100), recommended value: 30 or greater_ - The image characteristics at the boundary between the iris region and the pupil.
  * `iso_iris_pupil_ratio` _value range: (9.58, 121.30), recommended value: between 20 and 70_ - The degree to which the pupil is dilated or constricted.
  * `iso_iris_sclera_contrast` _value range: [0, 100), recommended value: greater than 5_ - The image characteristics at the boundary between the iris region and the sclera.
  * `iso_margin_adequacy` _value range: [0, 100], recommended value: greater than 80_ - The degree to which the iris portion of the image is centred relative to the edges of the entire image.
  * `iso_pupil_boundary_circularity` _value range: [0, 100]_ - The circularity of the iris-pupil boundary.    
  * `iso_sharpness` _value range: [0, 100)_ - The degree of focus present in the image.
  * `iso_usable_iris_area` _value range: [0, 100]_ - The fraction of the iris portion of the image that is not occluded by eyelids, eyelashes, or specular reflections.

The following normalized statistics are also reported by this provider:

  * `normalized_contrast` _value range: [0, 1]_ - Normalized value for `contrast`.
  * `normalized_iris_diameter` _value range: [0, 1]_ - Normalized value for `iris_diameter`.
  * `normalized_iris_pupil_gs` _value range: [0, 1]_ - Normalized value for `iris_pupil_gs`.
  * `normalized_iris_sclera_gs` _value range: [0, 1]_ - Normalized value for `iris_sclera_gs`.
  * `normalized_sharpness` _value range: [0, 1]_ - Normalized value for `sharpness`.

NEW! The following ISO normalized statistics are also reported by this provider 
  * `normalized_iso_greyscale_utilization` _value range: [0, 1]_ - Normalized value for `iso_greyscale_utilization`.
  * `normalized_iso_iris_diameter` _value range: [0, 1]_ - Normalized value for `iris_diameter`.
  * `normalized_iso_iris_pupil_concentricity` _value range: [0, 1]_ - Normalized value for `iso_iris_pupil_concentricity`.
  * `normalized_iso_iris_pupil_contrast` _value range: [0, 1]_ - Normalized value for `iso_iris_pupil_contrast`.
  * `normalized_iso_iris_pupil_ratio` _value range: [0, 1]_ - Normalized value for `iso_iris_pupil_ratio`.
  * `normalized_iso_iris_sclera_contrast` _value range: [0, 1]_ - Normalized value for `iso_iris_sclera_contrast`.
  * `normalized_iso_margin_adequacy` _value range: [0, 1]_ - Normalized value for `iso_margin_adequacy`.  
  * `normalized_iso_sharpness` _value range: [0, 1]_ - Normalized value for `iso_sharpness`.
  * `normalized_iso_usable_iris_area` _value range: [0, 1]_ - Normalized value for `iso_usable_iris_area`.
