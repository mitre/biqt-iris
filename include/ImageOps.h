// #######################################################################
// NOTICE
//
// This software (or technical data) was produced for the U.S. Government
// under contract, and is subject to the Rights in Data-General Clause
// 52.227-14, Alt. IV (DEC 2007).
//
// Copyright 2022 The MITRE Corporation. All Rights Reserved.
// #######################################################################

#ifndef IMAGE_OPS_H
#define IMAGE_OPS_H

#ifdef _WIN32
#define IMAGEOPS_EXPORT __declspec(dllexport)
#include <windows.h>
#else
#define IMAGEOPS_EXPORT
#endif

#define FOCUS_KERNEL_DIM 8

#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

struct Point {
  int X;
  int Y;
};

// Class for various image operations
class IMAGEOPS_EXPORT MFilter {
 public:
  MFilter();
  ~MFilter();
  void Initialize();

  // Main Functions
  int GetQualityFromImageFrame(const uint8_t *frame_bytes, int width, int height);

  // Quality Attributes
  int GetContrastScore();
  int GetDefocusScore();

  int GetIrisRadius();
  double GetISGSDiffMeanAvg();
  double GetIrisPupilGSDiff();
  int GetUsableIrisAreaPercent();

  int GetPupilRadius();
  double GetPupilCircularityDeviationAvg();

  // (Normalized) Quality Attributes
  double GetNContrast();
  double GetNDefocus();

  double GetNIrisID();
  double GetNISGSMean();
  double GetNIPGSDiff();
  double GetNIrisVis();

  // ISO Metrics
  double GetISOIrisScleraContrast();
  double GetISOIrisPupilContrast();
  double GetISOPupilBoundaryCircularity();
  double GetISOGreyscaleUtilization();
  double GetISOPIRatio();
  double GetISOIPConcentricity();
  double GetISOMarginAdequacy();
  double GetISOSharpness();
  double GetNormalizedISOIrisScleraContrast();
  double GetNormalizedISOIrisPupilContrast();
  double GetNormalizedISOGreyscaleUtilization();
  double GetNormalizedISOPIRatio();
  double GetNormalizedISOIPConcentricity();
  double GetNormalizedISOIrisDiameter();
  double GetNormalizedISOMarginAdequacy();
  double GetNormalizedISOSharpness();
  int GetIsoOverallQuality();

  // Features
  int GetIrisCenterX();
  int GetIrisCenterY();
  int GetPupilCenterX();
  int GetPupilCenterY();
  int min_width_ = 256;
  int min_height_ = 256;
  int max_width_ = 1000;
  int max_height_ = 680;

 private:
  void CreateImagePointers(int width, int height);
  void DeleteImagePointers(int width, int height);
  void InitMathTables();

  int CalcOverallQuality(int contrast, int defocus,
                         int isgs_mean, double margin, int iris_d,
                         double iris_vis, double ipgs_diff,
                         double iris_pupil_ratio);

  int Calc_ISO_Overall_Quality(double n_iso_sharpness_value_, double n_iso_greyscale_value_,
                                 double n_iso_ip_concentricity_value_, double n_iso_iris_sclera_contrast_value_,
                                 double n_iso_margin_adequacy_value_, double n_iso_iris_pupil_contrast_value_,
                                 double n_iso_iris_pupil_ratio_value_);
  int GetQuality(double n_contrast, double n_defocus, double n_iris_d,
                 double n_isgs_mean, double n_iris_vis, double n_margin,
                 double n_iris_pupil_ratio, double n_ipgs_diff);
  int Defocus(uint8_t **raw_img, int width, int height);
  int Contrast(uint8_t **raw_img, int width, int height);

  double NormalizeDefocus(int defocus);
  double NormalizeContrast(int contrast);
  double NormalizeIrisDiameter(int iris_d);
  double NormalizeISGS(int isgs_mean);
  double NormalizeIPGSDiff(double ipgs_diff);
  double NormalizeIrisPupilRatio(double iris_pupil_ratio);
  double NormalizeIrisVisibility(double iris_vis);

  void ISOContrast(const uint8_t *raw_img, const int width, const int height);
  void ISOPupilBoundaryCircularity(const double *pupil_radii, const int m, const int n);
  void ISOGreyscaleUtilization(const uint8_t *raw_img, const int width, const int height);
  void ISOIrisPupilConcentricity();
  void ISOMarginAdequacy(int width, int height);
  void ISOSharpness(const uint8_t *raw_img, const int width, const int height);

  void FindIris(const uint8_t *frame_bytes, int width, int height);
  void FindIris(uint8_t **raw_bytes, int width, int height);
  void FindIrisCenter(int **l_edge_vals, int **r_edge_vals, int width,
                      int height, uint8_t **raw_bytes);
  void FindFineIris(uint8_t **raw_bytes, int width, int height, int rough_iris_center_x,
                    int rough_iris_center_y, int rough_iris_diameter);
  void FindPupilCenter(int **edge_v, int width, int height, int iris_center_x,
                       int iris_center_y, int iris_diameter);
  void FindFinePupil(int rough_pupil_center_x, int rough_pupil_center_y, int width, int height);
  void FindOcclusions(uint8_t **raw_img, int width, int height, int iris_center_x, int iris_center_y,
                      int iris_radius, int pupil_center_x, int pupil_center_y, int pupil_radius);

  void RawByteFrameToTwoDimArray(const uint8_t *frame_bytes, uint8_t **raw_img, int width, int height);
  void CheckMargins(int width, int height, int iris_center_x, int iris_center_y, int iris_diameter);
  void DownSize(uint8_t **raw_img, int width, int height, int downsize_scale, uint8_t **ds_img);
  void EdgeMap(uint8_t **raw_img, int width, int height);
  void GetVertEdges(uint8_t **raw_img, int width, int height);

  int max_segs_ = 256;
  int max_radius_ = 250;
  int lower_seg_trial_ = 47;
  int upper_seg_trial_ = 35;
  int p_lower_seg_trial_ = 63;
  int p_upper_seg_trial_ = 35;
  int focus_kernel_[FOCUS_KERNEL_DIM][FOCUS_KERNEL_DIM] = {
      {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1},
      {-1, -1, 3, 3, 3, 3, -1, -1}, {-1, -1, 3, 3, 3, 3, -1, -1},
      {-1, -1, 3, 3, 3, 3, -1, -1}, {-1, -1, 3, 3, 3, 3, -1, -1},
      {-1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1}};

  bool math_tables_defined_;

  int max_point_response_;
  int min_r_trial_;
  int max_r_trial_;
  int min_p_r_trial_;
  int max_p_r_trial_;
  int pupil_cx_;
  int pupil_cy_;
  int pupil_rad_;
  int pupil_d_;
  int iris_center_x_;
  int iris_center_y_;
  int iris_radius_;
  int pupil_center_x_;
  int pupil_center_y_;
  int pupil_radius_;
  int overall_quality_;
  int iso_overall_quality_;
  int defocus_score_;
  int contrast_score_;
  int usable_iris_area_percent_;

  double combined_quality_;
  double overall_margin_;
  double isgs_diff_mean_avg_;
  double iris_pupil_gs_diff_;
  double iris_pupil_diameter_ratio_;
  double pupil_circularity_avg_deviation_;
  double iso_iris_sclera_contrast_value_;
  double iso_iris_pupil_contrast_value_;
  double iso_pupil_boundary_circularity_value_;
  double iso_greyscale_value_;
  double iso_ip_concentricity_value_;
  double iso_margin_adequacy_value_;
  double iso_sharpness_value_;
  double n_iso_iris_sclera_contrast_value_;
  double n_iso_iris_pupil_contrast_value_;
  double n_iso_iris_pupil_ratio_value_;
  double n_iso_pupil_boundary_circularity_value_;
  double n_iso_greyscale_value_;
  double n_iso_pupil_iris_contrast_value_;
  double n_iso_ip_concentricity_value_;
  double n_iso_margin_adequacy_value_;
  double n_iso_sharpness_value_;
  double n_iso_usable_iris_area_value_;
  double n_iso_iris_diameter_value_;
  double n_contrast_;
  double n_defocus_;
  double n_iris_id_;
  double n_isgs_mean_;
  double n_iris_vis_;
  double n_ipgs_diff_;
  double n_iris_pupil_ratio_;
  double n_defocus_factor_;
  double n_contrast_factor_;
  double n_i_diam_factor_;
  double n_isgs_factor_;
  double n_ipgs_diff_factor_;
  double n_iris_pupil_ratio_factor_;
  double n_iris_vis_factor_;
  double defocus_low_limit_;
  double defocus_upper_limit_;
  double contrast_low_limit_;
  double contrast_upper_limit_;
  double i_d_low_limit_;
  double i_d_med_limit_1_;
  double i_d_med_limit_2_;
  double i_d_upper_limit_;
  double isgs_mean_low_limit_;
  double isgs_mean_upper_limit_;
  double iris_vis_low_limit_;
  double iris_vis_upper_limit_;
  double ipgs_diff_low_limit_;
  double ipgs_diff_upper_limit_;
  double iris_pupil_ratio_low_limit_;
  double iris_pupil_ratio_high_limit_;
  double combined_quality_low_limit_;
  double combined_quality_upper_limit_;
  double dimless_r_[256];

  /** The angle in the range 0-255. 0 is 0 degrees. 64 is 90 deg. 128 is 180
   * deg 192 is 270 deg */
  int *isgs_vals_;
  int *cos;
  int **val_f16_edge_;
  int **ds_val_bin_;
  int **neg_ds_val_bin_;
  int **edge_val_bin_;
  int **pupil_edge_bin_;
  int **transposed_edge_bin_;
  int **vert_edge_;
  int ***trig_vals_;

  unsigned char **raw_bin_;
  unsigned char **transposed_raw_bin_;
  unsigned char **masked_raw_bin_;
  unsigned char **ds_raw_bin_;
  unsigned char **neg_ds_raw_bin_;
  unsigned char **edge_raw_bin_;
  unsigned char **vert_edge_bin_;
};

inline int MFilter::GetContrastScore() { return contrast_score_; }

inline int MFilter::GetDefocusScore() { return defocus_score_; }

inline int MFilter::GetIrisRadius() { return iris_radius_; }

inline double MFilter::GetISGSDiffMeanAvg() { return isgs_diff_mean_avg_; }

inline double MFilter::GetIrisPupilGSDiff() { return iris_pupil_gs_diff_; }

inline int MFilter::GetUsableIrisAreaPercent() { return usable_iris_area_percent_; }

inline int MFilter::GetPupilRadius() { return pupil_radius_; }

inline double MFilter::GetPupilCircularityDeviationAvg() { return pupil_circularity_avg_deviation_; }

inline double MFilter::GetNContrast() { return n_contrast_; }

inline double MFilter::GetNDefocus() { return n_defocus_; }

inline double MFilter::GetNIrisID() { return n_iris_id_; }

inline double MFilter::GetNISGSMean() { return n_isgs_mean_; }

inline double MFilter::GetNIPGSDiff() { return n_ipgs_diff_; }

inline double MFilter::GetNIrisVis() { return n_iris_vis_; }

inline double MFilter::GetISOIrisScleraContrast() { return iso_iris_sclera_contrast_value_; }

inline double MFilter::GetISOIrisPupilContrast() { return iso_iris_pupil_contrast_value_; }

inline double MFilter::GetISOPupilBoundaryCircularity() { return iso_pupil_boundary_circularity_value_; }

inline double MFilter::GetISOGreyscaleUtilization() { return iso_greyscale_value_; }

inline double MFilter::GetISOPIRatio() { return ((double) pupil_rad_ / (double) iris_radius_) * 100.0; }

inline double MFilter::GetISOIPConcentricity() { return iso_ip_concentricity_value_; }

inline double MFilter::GetISOMarginAdequacy() { return iso_margin_adequacy_value_; }

inline double MFilter::GetISOSharpness() { return iso_sharpness_value_; }

inline double MFilter::GetNormalizedISOSharpness() { return n_iso_sharpness_value_; }

inline double MFilter::GetNormalizedISOGreyscaleUtilization() { return n_iso_greyscale_value_; }

inline double MFilter::GetNormalizedISOIPConcentricity() { return n_iso_ip_concentricity_value_; }

inline double MFilter::GetNormalizedISOIrisDiameter() { return n_iso_iris_diameter_value_; }

inline double MFilter::GetNormalizedISOIrisScleraContrast() { return n_iso_iris_sclera_contrast_value_; }

inline double MFilter::GetNormalizedISOMarginAdequacy() { return n_iso_margin_adequacy_value_; }

inline double MFilter::GetNormalizedISOIrisPupilContrast() { return n_iso_iris_pupil_contrast_value_; }

inline double MFilter::GetNormalizedISOPIRatio() { return n_iso_iris_pupil_ratio_value_; }

inline int MFilter::GetIsoOverallQuality() { return iso_overall_quality_; }

inline int MFilter::GetIrisCenterX() { return iris_center_x_; }

inline int MFilter::GetIrisCenterY() { return iris_center_y_; }

inline int MFilter::GetPupilCenterX() { return pupil_center_x_; }

inline int MFilter::GetPupilCenterY() { return pupil_center_y_; }

/**
 *
 * @param a The first item to compare
 * @param b The second item to compare
 * @return An integer indicating the relationship between the 2 numbers. 0
 *      indicates equality, a negative number indicates b > a, and a positive
 *      number indicates b < a
 */
inline int compare(const void *a, const void *b) { return *(uint8_t *) a - *(uint8_t *) b; }

#endif
