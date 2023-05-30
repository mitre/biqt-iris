// #######################################################################
// NOTICE
//
// This software (or technical data) was produced for the U.S. Government
// under contract, and is subject to the Rights in Data-General Clause
// 52.227-14, Alt. IV (DEC 2007).
//
// Copyright 2022 The MITRE Corporation. All Rights Reserved.
// #######################################################################

#include "ImageOps.h"

MFilter::MFilter()
    : math_tables_defined_(false), defocus_score_(0),
      iris_center_x_(0), iris_center_y_(0), iris_radius_(0), usable_iris_area_percent_(0),
      contrast_low_limit_(30.0), contrast_upper_limit_(50.0),
      n_contrast_factor_(1), defocus_low_limit_(60), defocus_upper_limit_(80),
      n_defocus_factor_(1), i_d_low_limit_(160), i_d_med_limit_1_(200), i_d_med_limit_2_(350),
      i_d_upper_limit_(400), n_i_diam_factor_(1), isgs_mean_low_limit_(20),
      isgs_mean_upper_limit_(30), n_isgs_factor_(1), iris_vis_low_limit_(0.4),
      iris_vis_upper_limit_(0.7), ipgs_diff_low_limit_(15), ipgs_diff_upper_limit_(30),
      n_ipgs_diff_factor_(1), n_iris_pupil_ratio_(1),
      iris_pupil_ratio_low_limit_(0.20), iris_pupil_ratio_high_limit_(0.60),
      n_iris_pupil_ratio_factor_(1), n_iris_vis_factor_(1), combined_quality_low_limit_(0),
      combined_quality_upper_limit_(0.80), overall_quality_(0), min_r_trial_(20),
      max_r_trial_(54), min_p_r_trial_(8), max_p_r_trial_(32), pupil_cx_(0), pupil_cy_(0),
      pupil_rad_(5), pupil_d_(10), contrast_score_(0),
      overall_margin_(1.0), isgs_diff_mean_avg_(0.0),
      iris_pupil_gs_diff_(0.0), max_point_response_(1000), pupil_circularity_avg_deviation_(0.0),
      n_iso_sharpness_value_(1.0), n_iso_greyscale_value_(1.0),
      n_iso_margin_adequacy_value_(1.0), n_iso_iris_pupil_contrast_value_(1.0), n_iso_iris_pupil_ratio_value_(1.0) ,
      n_iso_ip_concentricity_value_(1.0), n_iso_iris_sclera_contrast_value_(1.0), iso_overall_quality_(100){

//  Updates for 01 March 2023 minor version
//  Increased max_r_trial from 42 to 54 to allow finding rough iris diameters up to 432 pixels
//  Decrease min_r_trial from 28 to 20 to allow finding rough iris diameters down to 80 pixels
  
  cos = new int[101];

  isgs_vals_ = new int[256];
  memset(isgs_vals_, 0, sizeof(int) * max_segs_);

  trig_vals_ = new int **[max_segs_ + 1];
  for (int i = 0; i < max_segs_ + 1; i++) {
    trig_vals_[i] = new int *[max_radius_ + 1];
    for (int j = 0; j < max_radius_ + 1; j++) {
      trig_vals_[i][j] = new int[2];

      memset(trig_vals_[i][j], 0, sizeof(int) * 2);
    }
  }
}

MFilter::~MFilter() {
  delete[] isgs_vals_;
  delete[] cos;

  for (int i = 0; i < max_segs_ + 1; i++) {
    for (int j = 0; j < max_radius_ + 1; j++) {
      delete[] trig_vals_[i][j];
    }
    delete[] trig_vals_[i];
  }
  delete[] trig_vals_;
}

void MFilter::Initialize() {
  if (!math_tables_defined_) {
    InitMathTables();
  }
}

int MFilter::GetQualityFromImageFrame(const uint8_t *frame_bytes, int width, int height) {
  int overall_quality = 0;

  try {
    // Create all the image-based pointers (using width / height information)
    CreateImagePointers(width, height);

    FindIris(frame_bytes, width, height);
    Defocus(raw_bin_, width, height);
    Contrast(raw_bin_, width, height);

    CheckMargins(width, height, iris_center_x_, iris_center_y_, iris_radius_ * 2);
    FindPupilCenter(edge_val_bin_, width / 4, height / 4, iris_center_x_ / 4, iris_center_y_ / 4, iris_radius_ * 2 / 4);
    FindOcclusions(raw_bin_,
                   width,
                   height,
                   iris_center_x_,
                   iris_center_y_,
                   iris_radius_,
                   pupil_cx_,
                   pupil_cy_,
                   pupil_rad_);

    ISOContrast(frame_bytes, width, height);
    ISOPupilBoundaryCircularity(dimless_r_, 17, sizeof(dimless_r_) / sizeof(dimless_r_[0]));
    ISOGreyscaleUtilization(frame_bytes, width, height);
    ISOIrisPupilConcentricity();
    ISOMarginAdequacy(width, height);
    ISOSharpness(frame_bytes, width, height);
    double iris_vis = usable_iris_area_percent_ / 100.0;
    overall_quality = CalcOverallQuality(contrast_score_, defocus_score_, (int) isgs_diff_mean_avg_,
                                         overall_margin_, iris_radius_ * 2, iris_vis,
                                         iris_pupil_gs_diff_, iris_pupil_diameter_ratio_);

    iso_overall_quality_ = Calc_ISO_Overall_Quality(n_iso_sharpness_value_, n_iso_greyscale_value_,
                                                    n_iso_ip_concentricity_value_, n_iso_iris_sclera_contrast_value_,
                                                    n_iso_margin_adequacy_value_, n_iso_iris_pupil_contrast_value_,
                                                    n_iso_iris_pupil_ratio_value_);
  } catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
  }

  // Delete all the image-based pointers (using width / height information)
  DeleteImagePointers(width, height);

  return overall_quality;
}

void MFilter::CreateImagePointers(int width, int height) {
  raw_bin_ = new uint8_t *[width];
  val_f16_edge_ = new int *[width];
  masked_raw_bin_ = new uint8_t *[width];
  for (int i = 0; i < width; i++) {
    raw_bin_[i] = new uint8_t[height];
    val_f16_edge_[i] = new int[height];
    masked_raw_bin_[i] = new uint8_t[height];

    memset(raw_bin_[i], 0, sizeof(uint8_t) * height);
    memset(val_f16_edge_[i], 0, sizeof(int) * height);
    memset(masked_raw_bin_[i], 0, sizeof(uint8_t) * height);
  }

  ds_raw_bin_ = new uint8_t *[(int) (width / 4)];
  edge_val_bin_ = new int *[(int) (width / 4)];
  pupil_edge_bin_ = new int *[(int) (width / 4)];
  edge_raw_bin_ = new uint8_t *[(int) (width / 4)];
  vert_edge_ = new int *[(int) (width / 4)];
  vert_edge_bin_ = new uint8_t *[(int) (width / 4)];
  for (int i = 0; i < (int) (width / 4); i++) {
    ds_raw_bin_[i] = new uint8_t[(int) (height / 4)];
    edge_val_bin_[i] = new int[(int) (height / 4)];
    pupil_edge_bin_[i] = new int[(int) (height / 4)];
    edge_raw_bin_[i] = new uint8_t[(int) (height / 4)];
    vert_edge_[i] = new int[(int) (height / 4)];
    vert_edge_bin_[i] = new uint8_t[(int) (height / 4)];

    memset(ds_raw_bin_[i], 0, sizeof(uint8_t) * height / 4);
    memset(edge_val_bin_[i], 0, sizeof(int) * height / 4);
    memset(pupil_edge_bin_[i], 0, sizeof(int) * height / 4);
    memset(edge_raw_bin_[i], 0, sizeof(uint8_t) * height / 4);
    memset(vert_edge_[i], 0, sizeof(int) * height / 4);
    memset(vert_edge_bin_[i], 0, sizeof(uint8_t) * height / 4);
  }
}

void MFilter::DeleteImagePointers(int width, int height) {
  for (int i = 0; i < width; i++) {
    delete[] raw_bin_[i];
    delete[] val_f16_edge_[i];
    delete[] masked_raw_bin_[i];
  }
  delete[] raw_bin_;
  delete[] val_f16_edge_;
  delete[] masked_raw_bin_;

  for (int i = 0; i < (int) (width / 4); i++) {
    delete[] ds_raw_bin_[i];
    delete[] edge_val_bin_[i];
    delete[] pupil_edge_bin_[i];
    delete[] edge_raw_bin_[i];
    delete[] vert_edge_[i];
    delete[] vert_edge_bin_[i];
  }
  delete[] ds_raw_bin_;
  delete[] edge_val_bin_;
  delete[] pupil_edge_bin_;
  delete[] edge_raw_bin_;
  delete[] vert_edge_;
  delete[] vert_edge_bin_;
}

void MFilter::InitMathTables() {
  int half_max_segs = max_segs_ / 2;
  int seg_90 = (int) ((double) max_segs_ / 4.0);
  int seg_180 = (int) ((double) max_segs_ / 2.0);
  double single_segment = (6.28318530718 / (double) max_segs_);

  for (int angle_seg_num = 0; angle_seg_num < seg_90 + 1; angle_seg_num++) {
    double d_sin = sin(single_segment * (double) (angle_seg_num));
    double d_cos = ::cos(single_segment * (double) (angle_seg_num));

    for (int radius_length = 0; radius_length < max_radius_ + 1; radius_length++) {
      if (d_cos >= 0) {
        trig_vals_[angle_seg_num][radius_length][0] = (int) floor(d_cos * (double) radius_length + 0.5);
      } else {
        trig_vals_[angle_seg_num][radius_length][0] = (int) ceil(d_cos * (double) radius_length - 0.5);
      }
      trig_vals_[half_max_segs - angle_seg_num][radius_length][0] = -trig_vals_[angle_seg_num][radius_length][0] - 1;
      trig_vals_[half_max_segs + angle_seg_num][radius_length][0] = -trig_vals_[angle_seg_num][radius_length][0] - 1;
      trig_vals_[max_segs_ - angle_seg_num][radius_length][0] = trig_vals_[angle_seg_num][radius_length][0];

      if (d_sin >= 0) {
        trig_vals_[angle_seg_num][radius_length][1] = (int) floor(d_sin * (double) radius_length + 0.5);
      } else {
        trig_vals_[angle_seg_num][radius_length][1] = (int) ceil(d_sin * (double) radius_length - 0.5);
      }
      trig_vals_[half_max_segs - angle_seg_num][radius_length][1] = trig_vals_[angle_seg_num][radius_length][1];
      trig_vals_[half_max_segs + angle_seg_num][radius_length][1] = -trig_vals_[angle_seg_num][radius_length][1] - 1;
      trig_vals_[max_segs_ - angle_seg_num][radius_length][1] = -trig_vals_[angle_seg_num][radius_length][1];
    }
  }
  trig_vals_[0][max_radius_][0] = max_radius_ - 1;
  trig_vals_[seg_90][max_radius_][1] = max_radius_ - 1;
  trig_vals_[half_max_segs][max_radius_][0] = -max_radius_;
  trig_vals_[max_segs_ - seg_90][max_radius_][1] = -max_radius_;
  trig_vals_[max_segs_][max_radius_][0] = max_radius_ - 1;

  math_tables_defined_ = true;
}

int MFilter::CalcOverallQuality(int contrast, int defocus,
                                int isgs_mean, double margin, int iris_d,
                                double iris_vis, double ipgs_diff,
                                double iris_pupil_ratio) {
  overall_quality_ = 0;
  // Normalize Contrast value
  NormalizeContrast(contrast);
  // Normalize Defocus score
  NormalizeDefocus(defocus);
  // Normalize ISGS score
  NormalizeISGS(isgs_mean);
  // Normalize IrisDiamter
  n_iris_id_ = NormalizeIrisDiameter(iris_d);
  // Normalize IPGSDifference
  n_ipgs_diff_ = NormalizeIPGSDiff(ipgs_diff);
  // Normalize IPDiamRatio
  n_iris_pupil_ratio_ = NormalizeIrisPupilRatio(iris_pupil_ratio);
  // Normalize IrisVisibility
  n_iris_vis_ = NormalizeIrisVisibility(iris_vis);
  // Calculate combined_quality_ Score
  overall_quality_ =
      GetQuality(n_contrast_, n_defocus_, n_iris_id_, n_isgs_mean_, n_iris_vis_, margin,
                 n_iris_pupil_ratio_, n_ipgs_diff_);
  // Finally Scale the combined_quality_ output metric
  return overall_quality_;
}

int MFilter::Calc_ISO_Overall_Quality(double n_iso_sharpness_value, double n_iso_greyscale_value, double n_iso_ip_concentricity_value, double  n_iso_iris_sclera_contrast_value,
                                      double n_iso_margin_adequacy_value, double n_iso_iris_pupil_contrast_value, double n_iso_iris_pupil_ratio_value) {
    double iso_overall_quality = 0.0;
    int iso_overall_quality_100 = 0;


    iso_overall_quality = n_iso_sharpness_value * n_iso_greyscale_value * n_iso_iris_sclera_contrast_value * n_iso_margin_adequacy_value *
                          n_iso_iris_pupil_contrast_value * n_iso_ip_concentricity_value;
    iso_overall_quality_100 = 100 * iso_overall_quality;
    return iso_overall_quality_100;
}


int MFilter::GetQuality(double n_contrast, double n_defocus, double n_iris_d,
                        double n_isgs_mean, double n_iris_vis, double n_margin,
                        double n_iris_pupil_ratio, double n_ipgs_diff) {
  combined_quality_ = (n_contrast * n_contrast_factor_ * n_defocus * n_defocus_factor_ *
      n_iris_d * n_isgs_mean * n_isgs_factor_ * n_iris_vis * n_iris_vis_factor_ *
      n_margin * n_iris_pupil_ratio * n_iris_pupil_ratio_factor_ *
      n_ipgs_diff * n_ipgs_diff_factor_);

  if (combined_quality_ < combined_quality_low_limit_) {
    overall_quality_ = 0;
  }
  if (combined_quality_ >= combined_quality_low_limit_ && combined_quality_ <= combined_quality_upper_limit_) {
    overall_quality_ = (int) (80.0 * (combined_quality_ - combined_quality_low_limit_)
        / (combined_quality_upper_limit_ - combined_quality_low_limit_));
  }
  if (combined_quality_ > combined_quality_upper_limit_) {
    overall_quality_ = (int) (80.0
        + 20.0 * (combined_quality_ - combined_quality_upper_limit_) / (1.0 - combined_quality_upper_limit_));
  }

  return overall_quality_;
}

int MFilter::Contrast(uint8_t **raw_img, int width, int height) {
  double total_c_square = 0;
  double total_c = 0;
  double avg_c = 0;

  try {
    double num_pix = (double) width * (double) height;

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        total_c += (double) raw_img[x][y];
      }
    }

    avg_c = total_c / num_pix;
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        double var_i = (double) raw_img[x][y] - avg_c;
        total_c_square += var_i * var_i;
      }
    }

    contrast_score_ = (int) sqrt(total_c_square / num_pix);
  } catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
    contrast_score_ = 0;
  }

  return contrast_score_;
}

int MFilter::Defocus(uint8_t **raw_img, int width, int height) {
  try {
    int k_width = (int) FOCUS_KERNEL_DIM;
    int k_height = (int) FOCUS_KERNEL_DIM;
    int mar_wid = k_width / 2;
    int mar_height = k_height / 2;
    int num_points = 0;
    int total_response = 0;
    for (int i_x = mar_wid * 4; i_x < width - mar_height * 4; i_x += 4) {
      for (int i_y = mar_height * 4; i_y < height - mar_height * 4; i_y += 4) {
        int pt_response = 0;
        for (int j_x = 0; j_x < k_width; j_x++) {
          for (int j_y = 0; j_y < k_height; j_y++) {
            int trial_pt = (int) raw_img[i_x - mar_wid + j_x][i_y - mar_height + j_y];
            pt_response = (trial_pt * focus_kernel_[j_x][j_y]) + pt_response;
          }
        }

        total_response += (int) abs(pt_response);
        num_points++;
      }
    }

    defocus_score_ = (int) (total_response / num_points);
  } catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
    defocus_score_ = 0;
  }

  return defocus_score_;
}

double MFilter::NormalizeDefocus(int defocus) {
  n_defocus_ = 0.0;
  if (defocus < defocus_low_limit_) {
    n_defocus_ = 0;
  }
  if (defocus >= defocus_low_limit_ && defocus <= defocus_upper_limit_) {
    n_defocus_ = ((double) (defocus) - defocus_low_limit_) / (defocus_upper_limit_ - defocus_low_limit_);
  }
  if (defocus >= defocus_upper_limit_) {
    n_defocus_ = 1;
  }

  return n_defocus_;
}

double MFilter::NormalizeContrast(int contrast) {
  n_contrast_ = 0.0;
  if (contrast < contrast_low_limit_) {
    n_contrast_ = 0;
  }
  if (contrast >= contrast_low_limit_ && contrast <= contrast_upper_limit_) {
    n_contrast_ = ((double) (contrast) - contrast_low_limit_) / (contrast_upper_limit_ - contrast_low_limit_);
  }
  if (contrast >= contrast_upper_limit_) {
    n_contrast_ = 1;
  }

  return n_contrast_;
}

double MFilter::NormalizeIrisDiameter(int iris_d) {
  n_iris_id_ = 0.0;
  if (iris_d < i_d_low_limit_) {
    n_iris_id_ = 0;
  }
  if (iris_d >= i_d_low_limit_ && iris_d <= i_d_med_limit_1_) {
    n_iris_id_ = (iris_d - i_d_low_limit_) / (i_d_med_limit_1_ - i_d_low_limit_);
  }
  if (iris_d >= i_d_med_limit_1_ && iris_d <= i_d_med_limit_2_) {
    n_iris_id_ = 1.0;
  }
  if (iris_d >= i_d_med_limit_2_ && iris_d <= i_d_upper_limit_) {
    n_iris_id_ = 1.0 - ((double) iris_d - i_d_med_limit_2_) / (i_d_upper_limit_ - i_d_med_limit_2_);
  }
  if (iris_d >= i_d_upper_limit_) {
    n_iris_id_ = 0;
  }

  return n_iris_id_;
}

double MFilter::NormalizeISGS(int isgs_mean) {
  n_isgs_mean_ = 0.0;
  if (isgs_mean < isgs_mean_low_limit_) {
    n_isgs_mean_ = 0.0; // Revision 2.2
  }
  if (isgs_mean >= isgs_mean_low_limit_ && isgs_mean <= isgs_mean_upper_limit_) {
    n_isgs_mean_ =
        0.8 + ((double) isgs_mean - isgs_mean_low_limit_) / (isgs_mean_upper_limit_ - isgs_mean_low_limit_) * 0.2;
  }
  if (isgs_mean >= isgs_mean_upper_limit_) {
    n_isgs_mean_ = 1;
  }

  return n_isgs_mean_;
}

double MFilter::NormalizeIPGSDiff(double ipgs_diff) {
  n_ipgs_diff_ = 0.0;
  if (ipgs_diff < ipgs_diff_low_limit_) {
    n_ipgs_diff_ = 0.85;
  }
  if (ipgs_diff >= ipgs_diff_low_limit_ && ipgs_diff <= ipgs_diff_upper_limit_) {
    double n_ipgs_diff_value_at_low_limit = 0.85;

    n_ipgs_diff_ = n_ipgs_diff_value_at_low_limit
        + (ipgs_diff - ipgs_diff_low_limit_) / (ipgs_diff_upper_limit_ - ipgs_diff_low_limit_)
            * (1.0 - n_ipgs_diff_value_at_low_limit);
  }
  if (ipgs_diff >= ipgs_diff_upper_limit_) {
    n_ipgs_diff_ = 1;
  }

  return n_ipgs_diff_;
}

double MFilter::NormalizeIrisPupilRatio(double iris_pupil_ratio) {
  n_iris_pupil_ratio_ = 0.0;
  if (iris_pupil_ratio < iris_pupil_ratio_low_limit_) {
    n_iris_pupil_ratio_ = 0.0;
  }
  if (iris_pupil_ratio >= iris_pupil_ratio_low_limit_ && iris_pupil_ratio <= iris_pupil_ratio_high_limit_) {
    n_iris_pupil_ratio_ = 1.0;
  }
  if (iris_pupil_ratio > iris_pupil_ratio_high_limit_) {
    n_iris_pupil_ratio_ = 0.0;
  }

  return n_iris_pupil_ratio_;
}

double MFilter::NormalizeIrisVisibility(double iris_vis) {
  n_iris_vis_ = 0.0;
  if (iris_vis < iris_vis_low_limit_) {
    n_iris_vis_ = 0;
  }
  if (iris_vis >= iris_vis_low_limit_ && iris_vis <= iris_vis_upper_limit_) {
    n_iris_vis_ = 0.7 + (iris_vis - iris_vis_low_limit_) / (iris_vis_upper_limit_ - iris_vis_low_limit_) * 0.3;
  }
  if (iris_vis >= iris_vis_upper_limit_) {
    n_iris_vis_ = 1;
  }

  return n_iris_vis_;
}

/**
 * Determines the ISO Iris/Sclera Contrast (6.2.2) and Pupil/Iris Contrast (6.2.3)
 * This must be called after FindFineIris and FindPupilCenter
 *
 * @param raw_img The raw image bytes as a 1 dimensional array
 * @param width The width of the image in pixels
 * @param height The height of the image in pixels
 */
void MFilter::ISOContrast(const uint8_t *raw_img, const int width, const int height) {
  if (iris_radius_ == 0 || pupil_rad_ == 0) {
    iso_iris_sclera_contrast_value_ = 0.0;
    iso_iris_pupil_contrast_value_ = 0.0;

    return;
  }

  const uint8_t pupil = 1;
  const uint8_t iris = 2;
  const uint8_t sclera = 3;
  const uint8_t iris_2 = 10;
  const uint8_t occluded = 255;

  // first we'll calculate the iris/sclera contrast
  // mark the pixels for use in median calculation
  uint8_t *pixel_marks = new uint8_t[width * height]();
  int **vals;

  const int inner_pupil_radius = pupil_rad_ * 0.8;
  for (int i = 0; i <= max_segs_; i++) {
    vals = trig_vals_[i];
    for (int j = 0; j <= inner_pupil_radius; j++) {
      int x = iris_center_x_ + vals[j][0];
      int y = iris_center_y_ - vals[j][1];
      if (x < 0 || y < 0 || x >= width || y >= height) {
        continue;
      }

      if (masked_raw_bin_[x][y] == 0xff) {
        pixel_marks[width * y + x] = occluded;
      } else {
        pixel_marks[width * y + x] = pupil;
      }
    }
  }

  const int inner_iris_radius = (pupil_rad_ + iris_radius_) / 2;
  const int outer_iris_radius = (int) (iris_radius_ * 0.9 + 0.5);
  for (int i = 0; i <= max_segs_; i++) {
    vals = trig_vals_[i];
    for (int j = inner_iris_radius; j <= outer_iris_radius; j++) {
      int x = iris_center_x_ + vals[j][0];
      int y = iris_center_y_ - vals[j][1];
      if (x < 0 || y < 0 || x >= width || y >= height) {
        continue;
      }

      if (masked_raw_bin_[x][y] == 0xff) {
        pixel_marks[width * y + x] = occluded;
      } else {
        pixel_marks[width * y + x] = iris;
      }
    }
  }
// Limit the
  const int inner_sclera_radius = (int) (iris_radius_ * 1.1 + 0.5);
  const int outer_sclera_radius = (int) (iris_radius_ * 1.2 + 0.5);
  for (int i = 0; i <= max_segs_; i++) {
    vals = trig_vals_[i];
    for (int j = inner_sclera_radius; j <= outer_sclera_radius; j++) {
      int x = iris_center_x_ + vals[j][0];
      int y = iris_center_y_ - vals[j][1];
      // Modified 2022-03-10 to prevent segment faults
      if (x < 0 || y < 0 || x >= width || y >= height) {
        continue;
      }

      if (masked_raw_bin_[x][y] == 0xff) {
        pixel_marks[width * y + x] = occluded;
      } else {
        pixel_marks[width * y + x] = sclera;
      }
    }
  }

  // put all the pupil, iris, and sclera pixels in a list for sorting
  uint8_t *pupil_pixels = new uint8_t[width * height]();
  uint8_t *iris_pixels = new uint8_t[width * height]();
  uint8_t *sclera_pixels = new uint8_t[width * height]();
  uint8_t *iris_pixel = iris_pixels;
  uint8_t *sclera_pixel = sclera_pixels;
  uint8_t *pupil_pixel = pupil_pixels;
  for (int i = 0; i < width * height; i++) {
    int x = pixel_marks[i];
    if (x == pupil) {
      *(pupil_pixel++) = raw_img[i];
    } else if (x == iris) {
      *(iris_pixel++) = raw_img[i];
    } else if (x == sclera) {
      *(sclera_pixel++) = raw_img[i];
    }
  }

  // sort the pixel values
  int pupil_count = pupil_pixel - pupil_pixels;
  int iris_count = iris_pixel - iris_pixels;
  int sclera_count = sclera_pixel - sclera_pixels;
  qsort((void *) pupil_pixels, pupil_count, sizeof(uint8_t), &compare);
  qsort((void *) iris_pixels, iris_count, sizeof(uint8_t), &compare);
  qsort((void *) sclera_pixels, sclera_count, sizeof(uint8_t), &compare);

  // pull the median from the sorted arrays;
  double pupil_median;
  if (pupil_count % 2) {
    pupil_median = (double) (pupil_pixels[pupil_count / 2] + pupil_pixels[pupil_count / 2 + 1]) / 2;
  } else {
    pupil_median = pupil_pixels[pupil_count / 2];
  }

  double iris_median;
  if (iris_count % 2) {
    iris_median = (double) (iris_pixels[iris_count / 2] + iris_pixels[iris_count / 2 + 1]) / 2;
  } else {
    iris_median = iris_pixels[iris_count / 2];
  }

  double sclera_median;
  if (sclera_count % 2) {
    sclera_median = (double) (sclera_pixels[sclera_count / 2] + sclera_pixels[sclera_count / 2 + 1]) / 2;
  } else {
    sclera_median = sclera_pixels[sclera_count / 2];
  }

  // now calculate the value
  double calc_divisor = sclera_median + iris_median - (2 * pupil_median);
  if (calc_divisor > 0) {
    iso_iris_sclera_contrast_value_ =
        100.0 * ((sclera_median - iris_median) / (sclera_median + iris_median - (2 * pupil_median)));
  } else {
    iso_iris_sclera_contrast_value_ = 0;
  }


  // take the absolute value
  if (iso_iris_sclera_contrast_value_ < 0) {
    // Limit the lower boundary of iso_iris_sclera_contrast_value_ to 0.
    iso_iris_sclera_contrast_value_ = 0.0;
  }
  
  n_iso_iris_sclera_contrast_value_ = 1.0;
  // calculate a normalized iso_iris_sclera_contrast value
  n_iso_iris_sclera_contrast_value_ = 1.0;
  // Default to 1.0;
  // Set limits for proportional value calculation
  // lower proportional limit is set to approximately lowest 2% of operational test set values
  // upper proportional limit is set to approximately lowest 5% of operational test set values
  double lpl_iso_iris_sclera_contrast = 0.0;
  // Substitute lower proportional limit below for stricter limits
  // double lpl_iso_iris_sclera_contrast = 50.0;
  double upl_iso_iris_sclera_contrast = 15.0;
  if (iso_iris_sclera_contrast_value_ < lpl_iso_iris_sclera_contrast){
      n_iso_iris_sclera_contrast_value_= 0.0;
  }
  if (iso_iris_sclera_contrast_value_ >= upl_iso_iris_sclera_contrast){
      n_iso_iris_sclera_contrast_value_= 1.0;
  }
  if (iso_iris_sclera_contrast_value_ >= lpl_iso_iris_sclera_contrast && iso_iris_sclera_contrast_value_ <= upl_iso_iris_sclera_contrast){
      n_iso_iris_sclera_contrast_value_= (iso_iris_sclera_contrast_value_ - lpl_iso_iris_sclera_contrast)/(upl_iso_iris_sclera_contrast - lpl_iso_iris_sclera_contrast);
  }

  // now calculate the pupil/iris contrast
  // we already have the pupil values, so just recalculate the iris ones.
  // first, mark the iris pixels for the calculation
  const int inner_iris_radius_p = (int) (iris_radius_ * 1.1); // inner and outer radii
  const int outer_iris_radius_p = (int) (iris_radius_ * 1.2); // for pupil-iris contrast
  for (int i = 0; i <= max_segs_; i++) {
    vals = trig_vals_[i];
    for (int j = inner_iris_radius_p; j <= outer_iris_radius_p; j++) {
      int x = iris_center_x_ + vals[j][0];
      int y = iris_center_y_ - vals[j][1];
      // Modified 2022-03-10 to prevent segfaults
      if (x < 0 || y < 0 || x >= width || y >= height) {
        continue;
      }

      if (masked_raw_bin_[x][y] != 0xff) {
        pixel_marks[width * y + x] = iris_2;
      }
    }
  }

  // put all the iris pixels in a list for sorting
  iris_pixel = iris_pixels;
  for (int i = 0; i < width * height; i++) {
    int x = pixel_marks[i];
    if (x == iris_2) {
      *(iris_pixel++) = raw_img[i];
    }
  }

  // sort the pixel values
  iris_count = iris_pixel - iris_pixels;
  qsort((void *) iris_pixels, iris_count, sizeof(uint8_t), &compare);

  // pull the median from the sorted array
  if (iris_count % 2) {
    iris_median = (double) (iris_pixels[iris_count / 2] + iris_pixels[iris_count / 2 + 1]) / 2;
  } else {
    iris_median = iris_pixels[iris_count / 2];
  }

  // now calculate the value
  double weber_ratio = (iris_median - pupil_median) / (20 + pupil_median);

  iso_iris_pupil_contrast_value_ = 100.0 * (weber_ratio / (0.75 + weber_ratio));
  // The spec doesn't mention an absolute value here, but a negative value would indicate something is wrong.
  // Limit lower boundary of iso_iris_pupil_contrast_value_ to 0.0.
  if (iso_iris_pupil_contrast_value_ < 0) {
    iso_iris_pupil_contrast_value_ = 0.0;
  }

  delete[] pupil_pixels;
  delete[] iris_pixels;
  delete[] sclera_pixels;
  delete[] pixel_marks;
}

/**
 * Determines the ISO Pupil Boundary Circularity (6.2.4)
 * 
 * @param pupil_radii The pupil radii as a 1 dimensional array
 * @param m The number of discrete Fourier coefficients to use
 * @param n The number of regularly spaced angular samples of radial gradient edge data
 */
void MFilter::ISOPupilBoundaryCircularity(const double *pupil_radii, const int m, const int n) {
  int num_valid_radii = 0;
  for (int theta = 0; theta < n; theta++) {
    if (pupil_radii[theta] > 0) {
      num_valid_radii++;
    }
  }

  if (num_valid_radii > 0) {
    // Pre-processing / validation step (i.e. only accept radii > 0)
    double valid_radii[num_valid_radii];
    memset(valid_radii, 0, num_valid_radii * sizeof(double));

    int j = 0, counter = 0;
    while (j < n) {
      if (pupil_radii[j] > 0) {
        valid_radii[counter++] = pupil_radii[j];
      }
      j++;
    }

    // Calculate ISO Pupil Boundary Circularity (using valid radii)
    std::complex<double> i(0, 1);

    double sum_of_squared_coefficients = 0;
    for (int k = 1; k < m; k++) {
      std::complex<double> coefficient(0, 0);
      for (int theta = 0; theta < num_valid_radii; theta++) {
        coefficient += valid_radii[theta] * std::exp((double(-2 * M_PI * k * theta) * i) / double(num_valid_radii));
      }
      sum_of_squared_coefficients += std::norm(coefficient);
    }

    iso_pupil_boundary_circularity_value_ = std::max(0.0, 100.0 - (sum_of_squared_coefficients / num_valid_radii));
  } else {
    iso_pupil_boundary_circularity_value_ = 0.0;
  }
}

/**
 * Determines the ISO Greyscale Utilization of the image (6.2.5)
 * 
 * @param raw_img The raw image bytes as a 1 dimensional array
 * @param width The width of the image in pixels
 * @param height The height of the image in pixels
 */
void MFilter::ISOGreyscaleUtilization(const uint8_t *raw_img, const int width, const int height) {
  iso_greyscale_value_ = 0.0;

  // Count the number of each greyscale value
  int probability[256];
  memset(probability, 0, 256 * sizeof(int));

  for (int j = 0; j < width * height; j++) {
    probability[raw_img[j]]++;
  }

  for (int j = 0; j < 256; j++) {
    if (probability[j]) {
      double i = (double) probability[j] / (double) (width * height);
      iso_greyscale_value_ -= i * log2(i);
    }
  }
}

/**
 * Determines ISO Pupil/Iris Concentricity (6.2.8)
 * This must be called after FindFineIris and FindPupilCenter
 */
void MFilter::ISOIrisPupilConcentricity() {
  if (iris_radius_ == 0 || pupil_rad_ == 0) {
    iso_ip_concentricity_value_ = 0.0;

    return;
  }

  int x_diff = pupil_cx_ - iris_center_x_;
  int y_diff = pupil_cy_ - iris_center_y_;

  iso_ip_concentricity_value_ = 1 - sqrt(x_diff * x_diff + y_diff * y_diff) / iris_radius_;
  if (iso_ip_concentricity_value_ < 0) {
    iso_ip_concentricity_value_ = 0.0;
  }

  iso_ip_concentricity_value_ *= 100.0;
}

/**
 * Determines the ISO margin adequacy of the image (6.2.9)
 * This must be called after FindFineIris and FindPupilCenter
 * 
 * @param width The width of the image in pixels
 * @param height The height of the image in pixels
 */
void MFilter::ISOMarginAdequacy(int width, int height) {
  if (iris_radius_ == 0 || pupil_rad_ == 0) {
    iso_margin_adequacy_value_ = 0.0;
    n_iso_margin_adequacy_value_ = 0.0;
    return;
  }

  double left_margin = (iris_center_x_ - iris_radius_) / (iris_radius_ * 0.6);
  if (left_margin < 0) {
    left_margin = 0;
  }

  double right_margin = (width - iris_center_x_ + iris_radius_) / (iris_radius_ * 0.6);
  if (right_margin < 0) {
    right_margin = 0;
  }

  double top_margin = (iris_center_y_ - iris_radius_) / (iris_radius_ * 0.2);
  if (top_margin < 0) {
    top_margin = 0;
  }

  double bottom_margin = (height - iris_center_y_ + iris_radius_) / (iris_radius_ * 0.2);
  if (bottom_margin < 0) {
    bottom_margin = 0;
  }
// March 2023 minor revisions
// Corrected bottom_margin check logic

  iso_margin_adequacy_value_ = left_margin;
  if (right_margin < iso_margin_adequacy_value_) {
    iso_margin_adequacy_value_ = right_margin;
  }
  if (top_margin < iso_margin_adequacy_value_) {
    iso_margin_adequacy_value_ = top_margin;
  }
  if (bottom_margin < iso_margin_adequacy_value_) {
    iso_margin_adequacy_value_ = bottom_margin;
  }

  iso_margin_adequacy_value_ *= 100.0;
  // Limit the upper boundary of reported iso_margin_adequacy_value to 100.
  if (iso_margin_adequacy_value_ > 100.0) {
    iso_margin_adequacy_value_ = 100.0;
  }
  n_iso_margin_adequacy_value_ = 1.0;
  
  if (iso_margin_adequacy_value_ < 80) {
    n_iso_margin_adequacy_value_ = 0.0;
  }
}

/**
 * Runs the ISO Standard sharpness algorithm (6.2.10)
 * 
 * @param raw_img The raw image bytes as a 1 dimensional array
 * @param width The width of the image in pixels
 * @param height The height of the image in pixels
 */
void MFilter::ISOSharpness(const uint8_t *raw_img, const int width, const int height) {
  int kernel[] = {0, 1, 1, 2, 2, 2, 1, 1, 0,
                  1, 2, 4, 5, 5, 5, 4, 2, 1,
                  1, 4, 5, 3, 0, 3, 5, 4, 1,
                  2, 5, 3, -12, -24, -12, 3, 5, 2,
                  2, 5, 0, -24, -40, -24, 0, 5, 2,
                  2, 5, 3, -12, -24, -12, 3, 5, 2,
                  1, 4, 5, 3, 0, 3, 5, 4, 1,
                  1, 2, 4, 5, 5, 5, 4, 2, 1,
                  0, 1, 1, 2, 2, 2, 1, 1, 0};

  int64_t ss = 0;
  for (int y = 4; y < height - 4; y += 4) {
    int64_t pixel_sum = 0;
    for (int x = 4; x < width - 4; x += 4) {
      int *kernel_value = kernel;
      for (int i = -4; i < 5; i++) {
        for (int j = -4; j < 5; j++) {
          int pixel_value = raw_img[(y * width + x) + (i * width + j)];
          pixel_sum += *(kernel_value++) * pixel_value;
        }
      }
    }
    ss += pixel_sum * pixel_sum;
  }

  // this constant is defined by the ISO standard (180000^2)
  const int64_t c = INT64_C(32400000000);

  double power = (double) ss / (double) (width * height);
  iso_sharpness_value_ = 100.0 * ((power * power) / (power * power + c));


  // calculate a normalized iso_sharpness
  // Default to 1.0;
  n_iso_sharpness_value_ = 1.0;
  // Set limits for proportional value calculation
  double lpl_iso_sharpness = 0.0;
  // Substitute lower proportional limit below for stricter limits
  double upl_iso_sharpness= 3.0;
  if (iso_sharpness_value_ < lpl_iso_sharpness){
    n_iso_sharpness_value_= 0.0;
  }
  if (iso_sharpness_value_ >= upl_iso_sharpness){
    n_iso_sharpness_value_= 1.0;
  }
  if (iso_sharpness_value_ >= lpl_iso_sharpness && iso_sharpness_value_ <= upl_iso_sharpness){
    n_iso_sharpness_value_= (iso_sharpness_value_ - lpl_iso_sharpness)/(upl_iso_sharpness - lpl_iso_sharpness);
  }
}

void MFilter::FindIris(const uint8_t *frame_bytes, int width, int height) {
  uint8_t **raw_img = new uint8_t *[width];
  for (int i = 0; i < width; i++) {
    raw_img[i] = new uint8_t[height];
  }

  RawByteFrameToTwoDimArray(frame_bytes, raw_img, width, height);
  FindIris(raw_img, width, height);

  for (int i = 0; i < width; i++) {
    delete[] raw_img[i];
  }

  delete[] raw_img;
}

void MFilter::FindIris(uint8_t **raw_bytes, int width, int height) {
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      raw_bin_[i][j] = raw_bytes[i][j];
    }
  }

  int ds_factor = 4;
  int ds_width = width / ds_factor;
  int ds_height = height / ds_factor;

  DownSize(raw_bytes, width, height, ds_factor, ds_raw_bin_);
  EdgeMap(ds_raw_bin_, ds_width, ds_height);
  GetVertEdges(ds_raw_bin_, ds_width, ds_height);
  FindIrisCenter(edge_val_bin_, edge_val_bin_, ds_width, ds_height, raw_bin_);
}

void MFilter::FindIrisCenter(int **l_edge_vals, int **r_edge_vals, int width,
                             int height, uint8_t **raw_bytes) {
  // Version 2.2.1 Change ... lower_seg_trial_ and upper_seg_trial_ are globals...
  //                          Reset these values due to conflict with Fine
  //                          Iris finding routine.
  lower_seg_trial_ = 47;
  upper_seg_trial_ = 40;

  int max_seg_trial = upper_seg_trial_;
  if (lower_seg_trial_ > max_seg_trial) {
    max_seg_trial = lower_seg_trial_;
  }

  Point iris_center_pt[224];
  int ctr = 0;
  for (int r = 2; r < 17; r += 2) {
    for (int j = 0; j < 61; j += 10) {
      iris_center_pt[ctr].X = trig_vals_[j][r][0];
      iris_center_pt[ctr].Y = trig_vals_[j][r][1];

      iris_center_pt[ctr + 56].X = trig_vals_[j][r][0];
      iris_center_pt[ctr + 56].Y = -trig_vals_[j][r][1];
      iris_center_pt[ctr + 112].X = -trig_vals_[j][r][0];
      iris_center_pt[ctr + 112].Y = trig_vals_[j][r][1];
      iris_center_pt[ctr + 168].X = -trig_vals_[j][r][0];
      iris_center_pt[ctr + 168].Y = -trig_vals_[j][r][1];
      ctr++;
    }
  }

  int best_response = 0;
  for (int r = min_r_trial_; r < max_r_trial_; r++, r++) {
    for (int i_x = width * 1 / 10; i_x < width * 9 / 10; i_x++, i_x++) {
      for (int i_y = height * 1 / 10; i_y < height * 9 / 10; i_y++, i_y++) {
// March 2023 minor revisions
// Increased search space for center of iris to within 10% of the left, right, top, and bottom of the image
        // Version 2.2.1 change initialize dark_iris_center to zero
        int dark_iris_center = 0;
        int trial_response = 0;
        int left_edge_response = 0;
        int right_edge_response = 0;
        int vert_edge_response = 0;
        int pt_dark_val = 0;
        int trial_x = 0;
        int trial_y = 0;

        for (int i_ctr = 0; i_ctr < 224; i_ctr++) {
          trial_x = i_x + iris_center_pt[i_ctr].X;
          trial_y = i_y - iris_center_pt[i_ctr].Y;

          if (trial_x > 0 && trial_x < width && trial_y < height && trial_y > 0) {
            pt_dark_val = (int) ds_raw_bin_[trial_x][trial_y];
            pt_dark_val = (55 - pt_dark_val);

            if (pt_dark_val < 0) {
              pt_dark_val = 0;
            }
            dark_iris_center += pt_dark_val;
          }
        }
        for (int j_x = 0; j_x < max_seg_trial; j_x++) {
          int r_x = i_x + trig_vals_[j_x][r][0];
          int u_y = i_y - trig_vals_[j_x][r][1];
          int l_x = i_x - trig_vals_[j_x][r][0];
          int l_y = i_y + trig_vals_[j_x][r][1];

          if (r_x + 4 < width - 1 && r_x - 2 >= 0) {
            if (u_y >= 4 && j_x < upper_seg_trial_) {
              right_edge_response += r_edge_vals[r_x][u_y];
              vert_edge_response += vert_edge_[r_x][u_y];
              right_edge_response += edge_val_bin_[r_x - 1][u_y];
              vert_edge_response += vert_edge_[r_x - 1][u_y];
              right_edge_response += edge_val_bin_[r_x + 1][u_y];
              vert_edge_response += vert_edge_[r_x + 1][u_y];
            }
            if (l_y < height - 3 && j_x < lower_seg_trial_) {
              right_edge_response += edge_val_bin_[r_x][l_y];
              vert_edge_response += vert_edge_[r_x][l_y];
              right_edge_response += edge_val_bin_[r_x - 1][l_y];
              vert_edge_response += vert_edge_[r_x - 1][l_y];
              right_edge_response += edge_val_bin_[r_x + 1][l_y];
              vert_edge_response += vert_edge_[r_x + 1][l_y];
            }
          }

          if (l_x - 2 >= 0 && r_x + 4 < width - 1) {
            if (u_y >= 4 && j_x < upper_seg_trial_) {
              left_edge_response += l_edge_vals[l_x][u_y];
              vert_edge_response += vert_edge_[l_x][u_y];
              left_edge_response += edge_val_bin_[l_x - 1][u_y];
              vert_edge_response += vert_edge_[l_x - 1][u_y];
              left_edge_response += edge_val_bin_[l_x + 1][u_y];
              vert_edge_response += vert_edge_[l_x + 1][u_y];
            }
            if (l_y < height - 1 && j_x < lower_seg_trial_) {
              left_edge_response += edge_val_bin_[l_x][l_y];
              vert_edge_response += vert_edge_[l_x][l_y];
              left_edge_response += edge_val_bin_[l_x - 1][l_y];
              vert_edge_response += vert_edge_[l_x - 1][l_y];
              left_edge_response += edge_val_bin_[l_x + 1][l_y];
              vert_edge_response += vert_edge_[l_x + 1][l_y];
            }
          }
        }

        dark_iris_center = dark_iris_center * 105 / 10;
        vert_edge_response = vert_edge_response * 3 / 10;
        trial_response =
            left_edge_response * 15 / 10 + right_edge_response * 15 / 10 + vert_edge_response + dark_iris_center;

        if (trial_response > best_response) {
          best_response = trial_response;
          iris_center_x_ = i_x * 4;
          iris_center_y_ = i_y * 4;
          iris_radius_ = r * 4;
        }
      }
    }
  }
  FindFineIris(raw_bin_, width * 4, height * 4, iris_center_x_, iris_center_y_, iris_radius_ * 2);
}

void MFilter::FindFineIris(uint8_t **raw_bytes, int width, int height, int rough_iris_center_x,
                           int rough_iris_center_y, int rough_iris_diameter) {
  lower_seg_trial_ = 40;
  upper_seg_trial_ = 30;

  int min_iris_center_x_fine = rough_iris_center_x - 4;
  int max_iris_center_x_fine = rough_iris_center_x + 4;
  int min_iris_center_y_fine = rough_iris_center_y - 4;
  int max_iris_center_y_fine = rough_iris_center_y + 4;
  int min_radius_fine = rough_iris_diameter / 2 - 4;
  int max_radius_fine = rough_iris_diameter / 2 + 4;
  if (max_radius_fine > max_radius_) {
    max_radius_fine = max_radius_;
  }

  Point tan_pts[33];
  int best_response = 0;
  for (int r = min_radius_fine; r < max_radius_fine; r++) {
    for (int i_x = min_iris_center_x_fine; i_x < max_iris_center_x_fine; i_x++) {
      for (int i_y = min_iris_center_y_fine; i_y < max_iris_center_y_fine; i_y++) {
        int right_edge_points_used = 0;
        int left_edge_points_used = 0;
        int trial_response = 0;
        int left_edge_response = 0;
        int right_edge_response = 0;

        for (int j_x = 0; j_x < lower_seg_trial_; j_x += 2) {
          int r_x = i_x + trig_vals_[j_x][r + 16][0];
          int u_y = i_y - trig_vals_[j_x][r + 16][1];
          int l_x = i_x - trig_vals_[j_x][r + 16][0];
          int l_y = i_y + trig_vals_[j_x][r + 16][1];

          int w_f = 1;
          if (r_x < width && l_x >= 0) {
            if (u_y >= 0 && j_x < upper_seg_trial_) {
              int point_response = 0;

              for (int r_ctr = 0; r_ctr < 16; r_ctr++, r_ctr++) {
                w_f = 4 - (r_ctr / 4);
                point_response += w_f *
                    ((int) raw_bytes[i_x + trig_vals_[j_x][r + r_ctr][0]][i_y - trig_vals_[j_x][r + r_ctr][1]] -
                        (int) raw_bytes[i_x + trig_vals_[j_x][r - r_ctr][0]][i_y - trig_vals_[j_x][r - r_ctr][1]]);
              }
              if (point_response < 0) {
                point_response = -point_response;
              }
              if (point_response > max_point_response_) {
                point_response = max_point_response_;
              }

              right_edge_response += point_response;
              right_edge_points_used++;
            }
            if (l_y < width && j_x < lower_seg_trial_) {
              int point_response = 0;

              for (int r_ctr = 0; r_ctr < 16; r_ctr++) {
                w_f = 4 - (r_ctr / 4);
                point_response += w_f *
                    ((int) raw_bytes[i_x + trig_vals_[j_x][r + r_ctr][0]][i_y + trig_vals_[j_x][r + r_ctr][1]] -
                        (int) raw_bytes[i_x + trig_vals_[j_x][r - r_ctr][0]][i_y + trig_vals_[j_x][r - r_ctr][1]]);
              }
              if (point_response < 0) {
                point_response = -point_response;
              }
              if (point_response > max_point_response_) {
                point_response = max_point_response_;
              }

              right_edge_response += point_response;
              right_edge_points_used++;
            }
          }
          if (l_x >= 0 && r_x <= width) {
            if (u_y >= 0 && j_x < upper_seg_trial_) {
              int point_response = 0;

              for (int r_ctr = 0; r_ctr < 16; r_ctr++) {
                w_f = 4 - (r_ctr / 4);
                point_response += w_f *
                    ((int) raw_bytes[i_x - trig_vals_[j_x][r + r_ctr][0]][i_y - trig_vals_[j_x][r + r_ctr][1]] -
                        (int) raw_bytes[i_x - trig_vals_[j_x][r - r_ctr][0]][i_y - trig_vals_[j_x][r - r_ctr][1]]);
              }
              if (point_response < 0) {
                point_response = -point_response;
              }
              if (point_response > max_point_response_) {
                point_response = max_point_response_;
              }

              left_edge_response += point_response;
              left_edge_points_used++;
            }
            if (l_y < height - 1 && j_x < lower_seg_trial_) {
              int point_response = 0;

              for (int r_ctr = 0; r_ctr < 16; r_ctr++) {
                w_f = 4 - (r_ctr / 4);
                point_response += w_f *
                    ((int) raw_bytes[i_x - trig_vals_[j_x][r + r_ctr][0]][i_y + trig_vals_[j_x][r + r_ctr][1]] -
                        (int) raw_bytes[i_x - trig_vals_[j_x][r - r_ctr][0]][i_y + trig_vals_[j_x][r - r_ctr][1]]);
              }
              if (point_response < 0) {
                point_response = -point_response;
              }
              if (point_response > max_point_response_) {
                point_response = max_point_response_;
              }

              left_edge_response += point_response;
              left_edge_points_used++;
            }
          }
        }
        trial_response = left_edge_response + right_edge_response;

        if (trial_response > best_response) {
          best_response = trial_response;
          iris_center_x_ = i_x;
          iris_center_y_ = i_y;
          iris_radius_ = r;
        }
      }
    }
  }

  int i_c_x = iris_center_x_;
  int i_c_y = iris_center_y_;
  int i_radius = iris_radius_;
  n_iso_iris_diameter_value_ = 1.0;
  if (iris_radius_ < 25) {
    n_iso_iris_diameter_value_ = 0.0;
  }
  if (iris_radius_ > max_radius_ - 16) {
    n_iso_iris_diameter_value_ = 0.0;
    throw new std::out_of_range("Iris Radius exceeds max allowable");
  }
  if (iris_radius_ < 17) {
    n_iso_iris_diameter_value_ = 0.0;
    throw new std::out_of_range("Invalid Iris Diameter value");
  }

  for (int j = 0; j < max_segs_; j++) {
    for (int r_ctr = 0; r_ctr < 33; r_ctr++) {
      tan_pts[r_ctr].X = i_c_x + trig_vals_[j][i_radius + r_ctr - 16][0];
      tan_pts[r_ctr].Y = i_c_y - trig_vals_[j][i_radius + r_ctr - 16][1];
    }

    int point_response = 0;
    for (int r_ctr = 1; r_ctr < 16; r_ctr++) {
      if (tan_pts[16 + r_ctr].X > 0 && tan_pts[16 + r_ctr].X < width && tan_pts[16 + r_ctr].Y > 0
          && tan_pts[16 + r_ctr].Y < height) {
        int outer_point_val = (int) raw_bytes[tan_pts[16 + r_ctr].X][tan_pts[16 + r_ctr].Y];
        int inner_point_val = (int) raw_bytes[tan_pts[16 - r_ctr].X][tan_pts[16 - r_ctr].Y];

        // Version 2.2.1 (Change to limit effects of single differences)
        int response = outer_point_val - inner_point_val;
        if (response > 75) {
          response = 75;
        }

        point_response = point_response + response;
      } else {
        point_response = 0;
      }
    }
    isgs_vals_[j] = point_response;
  }

  int isgs_total = 0;
  int pts_used = 0;
  for (int j_ctr = 0; j_ctr < 40; j_ctr++) {

    isgs_total = isgs_total + isgs_vals_[255 - j_ctr];
    pts_used++;
    isgs_total = isgs_total + isgs_vals_[128 + j_ctr];
    pts_used++;
  }
  for (int j_ctr = 0; j_ctr < 20; j_ctr++) {
    isgs_total = isgs_total + isgs_vals_[j_ctr];
    pts_used++;
    isgs_total = isgs_total + isgs_vals_[127 - j_ctr];
    pts_used++;
  }

  isgs_diff_mean_avg_ = isgs_total / (pts_used * 15);
  if (isgs_diff_mean_avg_< 0) {
      isgs_diff_mean_avg_ = 0;
  }
}

void MFilter::FindPupilCenter(int **edge_v, int width, int height, int iris_center_x,
                              int iris_center_y, int iris_diameter) {
  // Note that Find Pupil Center operates on a downsized ( 1/4 ) size edge image.
  // The iris_center_x and iris_center_y are already downscaled from full image coordinates
  // Modified 2022_04_14 to set default pupil center location to iris center location
  pupil_cx_ = iris_center_x * 4;
  pupil_cy_ = iris_center_y * 4;
  //default - Full scale default radius of 20 pixels
  pupil_rad_ = 5 * 4;
  //default - Full scale pupil diameter of 40 pixels
  pupil_d_ = pupil_rad_ * 2;
  // Only used if the trials based upon iris location fail to produce any values
  // End modification
  try {
    // Modified 2022_04_14 to use the iris_diameter passed to the function instead of a global value

    // Rough pupil finder uses the downscaled (1/4) edge image
    // iris_diameter is the downscaled iris diameter passed to the rough pupil finding routine
    // Limit the possible pupil radius to 60% of the iris diameter
    max_p_r_trial_ = iris_diameter / 2 * 6 / 10;

    // only use 55 angular segments of 256/4 (64) for each side of the pupil
    p_upper_seg_trial_ = 55;
    p_lower_seg_trial_ = 55;
    int p_max_seg_trial = p_upper_seg_trial_;
    if (p_lower_seg_trial_ > p_max_seg_trial) {
      p_max_seg_trial = p_lower_seg_trial_;
    }

    // define the region of interest in the pupil search to be at least 10% of the pupil diameter from the image margins
    int iris_reduced = iris_diameter / 10;

    int roi_top = iris_center_y - iris_reduced;
    if (roi_top < 0) {
      roi_top = 0;
    }
    int roi_left = iris_center_x - iris_reduced;
    if (roi_left < 0) {
      roi_left = 0;
    }
    int roi_right = iris_center_x + iris_reduced;
    if (roi_right > width - 1) {
      roi_right = width - 1;
    }
    int roi_bottom = iris_center_y + iris_reduced;
    if (roi_bottom > height - 1) {
      roi_bottom = height - 1;
    }

    Point pupil_center_pt[224];
    int p_ctr = 0;
    for (int r = 3; r < 25; r += 3) {
      for (int j = 0; j < 61; j += 10) {
        pupil_center_pt[p_ctr].X = trig_vals_[j][r][0];
        pupil_center_pt[p_ctr].Y = trig_vals_[j][r][1];
        pupil_center_pt[p_ctr + 56].X = trig_vals_[j][r][0];
        pupil_center_pt[p_ctr + 56].Y = -trig_vals_[j][r][1];
        pupil_center_pt[p_ctr + 112].X = -trig_vals_[j][r][0];
        pupil_center_pt[p_ctr + 112].Y = trig_vals_[j][r][1];
        pupil_center_pt[p_ctr + 168].X = -trig_vals_[j][r][0];
        pupil_center_pt[p_ctr + 168].Y = -trig_vals_[j][r][1];
        p_ctr++;
      }
    }

    int best_response = 0;
    int ipgs_response = 0;
    int ipgs_points_used = 0;
    int pt_val = 0;
    for (int r = min_p_r_trial_; r < max_p_r_trial_; r++) {
      for (int i_x = roi_left; i_x < roi_right; i_x++) {
        for (int i_y = roi_top; i_y < roi_bottom; i_y++) {
          int dark_pupil_val = 0;
          int pt_dark_val = 0;
          int trial_x = 0;
          int trial_y = 0;

          for (p_ctr= 0; p_ctr < 224; p_ctr++) {
            trial_x = i_x + pupil_center_pt[p_ctr].X;
            trial_y = i_y - pupil_center_pt[p_ctr].Y;

            if (trial_x > 0 && trial_x < width && trial_y < height && trial_y > 0) {
              pt_val = (int) ds_raw_bin_[trial_x][trial_y];
              // experimental kb
              // Limit pupil values to max_pupil_pt_intensity to reduce effects of specularities with 255 intensity
              int max_pupil_pt_intensity = 100;
              if (pt_val > max_pupil_pt_intensity) {
                pt_val = max_pupil_pt_intensity;
              }
              // pt_dark_val is the difference between the maximum_pupil_pt_intensity and the intensity of the pupil point
              // this value is higher if the point is darker
              pt_dark_val = max_pupil_pt_intensity - pt_dark_val;

              if (pt_dark_val < 0) {
                pt_dark_val = 0;
              }
              // Use a weighting factor of 1/5 to limit the cumulative response held by dark_pupil_val (accumulator)
              pt_dark_val /= 5;
              dark_pupil_val = dark_pupil_val + pt_dark_val;
            }
          }
          ipgs_response = 0;
          ipgs_points_used = 0;

          int right_edge_points_used = 0;
          int left_edge_points_used = 0;
          int trial_response = 0;
          int left_edge_response = 0;
          int right_edge_response = 0;

          for (int j = 0; j < p_max_seg_trial; j++, j++) {
            int r_x = i_x + trig_vals_[j][r][0];
            int u_y = i_y - trig_vals_[j][r][1];
            int l_x = i_x - trig_vals_[j][r][0];
            int l_y = i_y + trig_vals_[j][r][1];

            if (r_x + 4 < width - 1 && r_x - 4 >= 0) {
              if (u_y >= 0 && j < p_upper_seg_trial_) {
                right_edge_response += 2 * edge_v[r_x][u_y];
                right_edge_response += edge_v[r_x - 1][u_y];
                right_edge_response += edge_v[r_x + 1][u_y];
                right_edge_points_used++;
              }
              if (l_y < height - 1 && j < p_lower_seg_trial_) {
                right_edge_response += 2 * edge_v[r_x][l_y];
                right_edge_response += edge_v[r_x - 1][l_y];
                right_edge_response += edge_v[r_x + 1][l_y];
                right_edge_points_used++;
              }
            }
            if (l_x - 4 >= 0 && r_x + 4 <= width - 1) {
              if (u_y >= 0 && j < p_upper_seg_trial_) {
                left_edge_response += 2 * edge_v[l_x][u_y];
                left_edge_response += edge_v[l_x - 1][u_y];
                left_edge_response += edge_v[l_x + 1][u_y];
                left_edge_points_used++;
              }
              if (l_y < height - 1 && j < p_lower_seg_trial_) {
                left_edge_response += 2 * edge_v[l_x][l_y];
                left_edge_response += edge_v[l_x - 1][l_y];
                left_edge_response += edge_v[l_x + 1][l_y];
                left_edge_points_used++;
              }
            }
          }

          int edge_response;
          if (left_edge_points_used + right_edge_points_used) {
            edge_response =
                (left_edge_response + right_edge_response) * 100 / (left_edge_points_used + right_edge_points_used);
          } else {
            edge_response = 0;
          }

          trial_response = edge_response + dark_pupil_val;
          if (trial_response > best_response) {
            best_response = trial_response;
            pupil_cx_ = (i_x - 1) * 4;
            pupil_cy_ = (i_y) * 4;
            pupil_rad_ = (r) * 4;
            pupil_d_ = r * 8;
          }

          int pup_ic_dist = (int) sqrt((double) ((iris_center_y - i_y) * (iris_center_y - i_y)
              + (iris_center_x - i_x) * (iris_center_x - i_x)));
          if ((pup_ic_dist + r) > (iris_diameter / 2) * 65 / 100
              || (i_y + r) > (iris_center_y + (iris_diameter / 2) * 65 / 100)
              || (i_y - r) < (i_y - (iris_center_y / 2) * 65 / 100)) {
            i_y = roi_bottom;
          }
        }

        if ((i_x + r) > (iris_center_x + (iris_diameter / 2) * 65 / 100)
            || (iris_center_x - r) < (i_x - (iris_diameter / 2) * 65 / 100)) {
          i_x = roi_right;
        }
      }
    }

    int rough_pupil_center_x = pupil_cx_;
    int rough_pupil_center_y = pupil_cy_;
    int rough_pupil_radius = pupil_rad_;
    int rough_pupil_diameter = pupil_rad_ * 2;

    FindFinePupil(rough_pupil_center_x, rough_pupil_center_y, width * 4, height * 4);

    ipgs_response = 0;
    ipgs_points_used = 0;
    int i_x = pupil_cx_ / 4;
    int i_y = pupil_cy_ / 4;
    int r = pupil_rad_ / 4;
    for (int j_x = -20; j_x < 20; j_x++) {
      for (int j_y = -5; j_y < 5; j_y++) {
        if ((i_x + j_x - r) > 0 && ((i_x + j_x + r) < (width - 1)) && (i_y + j_y) > 0 && ((i_y + j_y) < (height - 1))) {
          ipgs_response += edge_v[i_x - r + j_x][i_y + j_y];
          ipgs_points_used++;
          ipgs_response += edge_v[i_x + r + j_x][i_y + j_y];
          ipgs_points_used++;
        }
      }
    }

    iris_pupil_gs_diff_ = (double) ipgs_response / (double) ipgs_points_used;
    iris_pupil_diameter_ratio_ = (double) pupil_rad_ / (double) iris_radius_;
  } catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
    pupil_cx_ = iris_center_x_;
    pupil_cy_ = iris_center_y;
    pupil_rad_ = 30;
    pupil_d_ = iris_radius_ * 2;
    iris_pupil_diameter_ratio_ = .2;
    iris_pupil_gs_diff_ = 0;
  }
}

void MFilter::FindFinePupil(int rough_pupil_center_x, int rough_pupil_center_y, int width, int height) {
  try {
    int min_pupil_center_x_fine = pupil_cx_ - 16;
    int max_pupil_center_x_fine = pupil_cx_ + 16;
    int min_pupil_center_y_fine = pupil_cy_ - 16;
    int max_pupil_center_y_fine = pupil_cy_ + 16;
    int min_radius_fine = 16;
    int max_radius_fine = pupil_rad_ + 16;
    if (max_radius_fine > max_p_r_trial_ * 4) {
      max_radius_fine = max_p_r_trial_ * 4;
    }

    int seg_best_r[256];
    for (int counter = 0; counter < 256; counter++) {
      seg_best_r[counter] = 0;
    }

    int best_response = 0;
    for (int r = min_radius_fine; r < max_radius_fine; r++) {
      for (int i_x = min_pupil_center_x_fine; i_x < max_pupil_center_x_fine; i_x++) {
        for (int i_y = min_pupil_center_y_fine; i_y < max_pupil_center_y_fine; i_y++) {
          int trial_response = 0;

          for (int j = 0; j < 256; j++) {
            int point_response = 0;
            seg_best_r[j] = 0;

            int r_x = i_x + (int) trig_vals_[j][r + 16][0];
            int u_y = i_y - (int) trig_vals_[j][r + 16][1];

            if (r_x < width && r_x >= 0 && u_y >= 0 && u_y < height) {
              int outer_point_val = 0;
              int inner_point_val = 0;
              int inner_point_x = 0;
              int inner_point_y = 0;
              int outer_point_x = 0;
              int outer_point_y = 0;

              int max_pupil_point_intensity = 50;
              int max_iris_point_intensity = 220;
              int partial_point_response = 0;

              for (int r_ctr = 1; r_ctr < 8; r_ctr++, r_ctr++) {
                int trig_val_inner_x = trig_vals_[j][r - r_ctr][0];
                int trig_val_outer_x = trig_vals_[j][r + r_ctr][0];
                int trig_val_outer_y = -trig_vals_[j][r + r_ctr][1];
                int trig_val_inner_y = trig_vals_[j][r + r_ctr][1];

                outer_point_x = i_x + trig_vals_[j][r + r_ctr][0];
                inner_point_x = i_x + trig_vals_[j][r - r_ctr][0];
                outer_point_y = i_y - trig_vals_[j][r + r_ctr][1];
                inner_point_y = i_y + trig_vals_[j][r - r_ctr][1];
                if (outer_point_x > 0 && outer_point_x < width && outer_point_y > 0 && outer_point_y < height &&
                    inner_point_x > 0 && inner_point_x < width && inner_point_y > 0 && inner_point_y < height) {
                  outer_point_val = (int) raw_bin_[outer_point_x][outer_point_y];
                  inner_point_val = (int) raw_bin_[inner_point_x][inner_point_y];
                } else {
                  outer_point_val = 0;
                  inner_point_val = 0;
                }
                if (outer_point_val < max_iris_point_intensity && inner_point_val < max_pupil_point_intensity) {
                  partial_point_response = outer_point_val - inner_point_val;
                  if (partial_point_response > 0) {
                    point_response = point_response + partial_point_response;
                  }
                }
              }
            }
            trial_response = trial_response + point_response;
          }
          if (trial_response > best_response) {
            best_response = trial_response;
            pupil_cx_ = i_x;
            pupil_cy_ = i_y;
            pupil_rad_ = r;
          }

          std::string i_x_str = std::to_string(i_x);
          std::string i_y_str = std::to_string(i_y);
          std::string r_str = std::to_string(r);
          std::string radius_str = std::to_string(pupil_rad_);
          std::string pupil_center_x_str = std::to_string(pupil_cx_);
          std::string pupil_center_y_str = std::to_string(pupil_cy_);
          std::string response_str = std::to_string(trial_response);
          std::string best_response_str = std::to_string(best_response);
        }
      }

    }
    int pupil_center_x = pupil_cx_;
    int pupil_center_y = pupil_cy_;
    int pupil_radius = pupil_rad_;
    pupil_center_x_ = pupil_cx_;
    pupil_center_y_ = pupil_cy_;
    pupil_radius_ = pupil_rad_;

    for (int j = 0; j < 256; j++) {
      int best_radial_response = 0;
      int inner_point_x = 0;
      int inner_point_y = 0;
      int outer_point_x = 0;
      int outer_point_y = 0;
      int outer_point_val = 120;
      int inner_point_val = 10;
      int max_pupil_point_intensity = 50;
      int max_iris_point_intensity = 150;
      int err_x_negative_value_exceeded = 0;
      int err_width_exceeded = 0;
      int err_y_negative_value_exceeded = 0;
      int err_height_exceeded = 0;
      for (int r = min_radius_fine; r < max_radius_fine; r++) {
        int point_response = 0;

        for (int rr = 1; rr < 16; rr++) {
          outer_point_x = pupil_cx_ + trig_vals_[j][r + rr][0];
          inner_point_x = pupil_cx_ - trig_vals_[j][r - rr][0];
          outer_point_y = pupil_cy_ - trig_vals_[j][r + rr][1];
          inner_point_y = pupil_cy_ + trig_vals_[j][r + rr][1];
          // modified 2022_03_10 to avoid segfaults
          if (outer_point_x >= 0 && outer_point_x < width && outer_point_y >= 0 && outer_point_y < height
              && inner_point_x >= 0 && inner_point_y >= 0) {
            outer_point_val = (int) raw_bin_[outer_point_x][outer_point_y];
            inner_point_val = (int) raw_bin_[inner_point_x][inner_point_y];
            if (outer_point_val < max_iris_point_intensity && inner_point_val < max_pupil_point_intensity) {
              point_response = point_response + (outer_point_val - inner_point_val);
            }
          }
          // end of modified code

        }
        if (point_response > best_radial_response) {
          seg_best_r[j] = r;
          best_radial_response = point_response;
        }
      }
    }

    int num_angles_used = 0;
    double total_deviation = 0.0;
    for (int j = 0; j < 256; j++) {
      dimless_r_[j] = 0.0;

      if (seg_best_r[j] > 0) {
        dimless_r_[j] = (double) seg_best_r[j] / (double) pupil_radius;
        double deviation = dimless_r_[j] - 1.0f;
        total_deviation = total_deviation + std::abs(deviation);
        num_angles_used++;
      }
    }

    if (num_angles_used > 0) {
      pupil_circularity_avg_deviation_ = total_deviation / (double) num_angles_used;
    }
  } catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
  }
}

void MFilter::FindOcclusions(uint8_t **raw_img, int width, int height, int iris_center_x, int iris_center_y,
                             int iris_radius, int pupil_center_x, int pupil_center_y, int pupil_radius) {
  int edge_lim = 80;
  int spec_lim = 225;
  int spec_def_val = 200;
  double z_val = 5.8f;

  int upper_right_gray[100];
  int upper_left_gray[100];
  int lower_left_gray[100];
  int lower_right_gray[100];
  int u_r_g_pts[100];
  int u_l_g_pts[100];
  int l_r_g_pts[100];
  int l_l_g_pts[100];
  int u_l_occ_pts[100];
  int u_r_occ_pts[100];
  int l_l_occ_pts[100];
  int l_r_occ_pts[100];

  Point u_l_g_start[100];
  Point u_r_g_start[100];
  Point l_l_g_start[100];
  Point l_r_g_start[100];
  Point u_l_g_stop[100];
  Point u_r_g_stop[100];
  Point l_l_g_stop[100];
  Point l_r_g_stop[100];

  memset(upper_right_gray, 0, sizeof(int) * 100);
  memset(upper_left_gray, 0, sizeof(int) * 100);
  memset(lower_left_gray, 0, sizeof(int) * 100);
  memset(lower_right_gray, 0, sizeof(int) * 100);
  memset(u_r_g_pts, 0, sizeof(int) * 100);
  memset(u_l_g_pts, 0, sizeof(int) * 100);
  memset(l_r_g_pts, 0, sizeof(int) * 100);
  memset(l_l_g_pts, 0, sizeof(int) * 100);
  memset(u_r_occ_pts, 0, sizeof(int) * 100);
  memset(u_l_occ_pts, 0, sizeof(int) * 100);
  memset(l_r_occ_pts, 0, sizeof(int) * 100);
  memset(l_l_occ_pts, 0, sizeof(int) * 100);
  memset(u_l_g_start, 0, sizeof(Point) * 100);
  memset(u_r_g_start, 0, sizeof(Point) * 100);
  memset(l_l_g_start, 0, sizeof(Point) * 100);
  memset(l_r_g_start, 0, sizeof(Point) * 100);
  memset(u_l_g_stop, 0, sizeof(Point) * 100);
  memset(u_r_g_stop, 0, sizeof(Point) * 100);
  memset(l_l_g_stop, 0, sizeof(Point) * 100);
  memset(l_r_g_stop, 0, sizeof(Point) * 100);

  for (int c_ctr = 0; c_ctr < width; c_ctr++) {
    for (int d_ctr = 0; d_ctr < height; d_ctr++) {
      masked_raw_bin_[c_ctr][d_ctr] = raw_bin_[c_ctr][d_ctr];
    }
  }

  int pupil_radius_squared = pupil_radius * pupil_radius;
  for (int i_ctr = 0; i_ctr < 100; i_ctr++) {
    int i_rel_x = (int) floor(0.5 + sqrt((double) (10000 - i_ctr * i_ctr)));

    for (int x = 0; x < i_rel_x - 8; x++) {
      int i_y_upper = iris_center_y - iris_radius * i_ctr / 100;
      int i_y_lower = iris_center_y + iris_radius * i_ctr / 100;

      int x_chord_pt = iris_center_x + x * iris_radius / 100;

      if (x_chord_pt < width) {
        int x_test = x_chord_pt - pupil_center_x;

        if (i_y_upper >= 0) {
          int y_test = i_y_upper - pupil_center_y;

          if ((x_test * x_test) + (y_test * y_test) > pupil_radius_squared + 20 * 20) {
            if (x_chord_pt > 0 && x_chord_pt < width && i_y_upper > 0 && i_y_upper < height - 1) {
              int pt_val = (int) raw_img[x_chord_pt][i_y_upper];
              if (pt_val > spec_lim) {
                pt_val = spec_def_val;
              }

              upper_right_gray[i_ctr] = upper_right_gray[i_ctr] + pt_val;
              u_r_g_pts[i_ctr]++;
              if (u_r_g_start[i_ctr].X == 0) {
                u_r_g_start[i_ctr].X = x_chord_pt;
                u_r_g_start[i_ctr].Y = i_y_upper;
              }
              u_r_g_stop[i_ctr].X = x_chord_pt;
              u_r_g_stop[i_ctr].Y = i_y_upper;

              int pt_edge = (int) vert_edge_bin_[x_chord_pt / 4][i_y_upper / 4] + val_f16_edge_[x_chord_pt][i_y_upper];
              if ((int) raw_img[x_chord_pt][i_y_upper] > spec_lim || pt_edge > edge_lim) {
                masked_raw_bin_[x_chord_pt][i_y_upper] = 0xFF;
                u_r_occ_pts[i_ctr]++;
              }
            }
          }
        }
        if (i_y_lower <= height - 1) {
          int y_test = i_y_lower - pupil_center_y;

          if ((x_test * x_test) + (y_test * y_test) > pupil_radius_squared + 20 * 20) {
            if (x_chord_pt > 0 && x_chord_pt < width && i_y_lower > 0 && i_y_lower < height) {
              int pt_val = (int) raw_img[x_chord_pt][i_y_lower];
              if (pt_val > spec_lim) {
                pt_val = spec_def_val;
              }

              lower_right_gray[i_ctr] = lower_right_gray[i_ctr] + pt_val;
              l_r_g_pts[i_ctr]++;
              if (l_r_g_start[i_ctr].X == 0) {
                l_r_g_start[i_ctr].X = x_chord_pt;
                l_r_g_start[i_ctr].Y = i_y_lower;
              }
              l_r_g_stop[i_ctr].X = x_chord_pt;
              l_r_g_stop[i_ctr].Y = i_y_lower;

              int pt_edge = (int) vert_edge_bin_[x_chord_pt / 4][i_y_lower / 4]
                  + (int) edge_raw_bin_[x_chord_pt / 4][i_y_lower / 4];
              if (raw_img[x_chord_pt][i_y_lower] > spec_lim || pt_edge > edge_lim) {
                masked_raw_bin_[x_chord_pt][i_y_lower] = 0xFF;
                l_r_occ_pts[i_ctr]++;
              }
            }
          }
        }
      }
    }

    for (int x = 0; x < i_rel_x - 8; x++) {
      int i_y_upper = iris_center_y - iris_radius * i_ctr / 100;
      int i_y_lower = iris_center_y + iris_radius * i_ctr / 100;

      int x_chord_pt = iris_center_x - x * iris_radius / 100;

      if (x_chord_pt < width - 1) {
        int x_test = x_chord_pt - pupil_center_x;

        if (i_y_upper >= 0) {
          int y_test = i_y_upper - pupil_center_y;

          if ((x_test * x_test) + (y_test * y_test) > pupil_radius_squared + 20 * 20) {
            if (x_chord_pt > 0 && x_chord_pt < width && i_y_upper > 0 && i_y_upper < height) {
              int pt_val = (int) raw_img[x_chord_pt][i_y_upper];
              if (pt_val > spec_lim) {
                pt_val = spec_def_val;
              }

              upper_left_gray[i_ctr] = upper_left_gray[i_ctr] + pt_val;
              u_l_g_pts[i_ctr]++;
              if (u_l_g_stop[i_ctr].X == 0) {
                u_l_g_stop[i_ctr].X = x_chord_pt;
                u_l_g_stop[i_ctr].Y = i_y_upper;
              }
              u_l_g_start[i_ctr].X = x_chord_pt;
              u_l_g_start[i_ctr].Y = i_y_upper;

              int pt_edge = (int) vert_edge_bin_[x_chord_pt / 4][i_y_upper / 4]
                  + (int) edge_raw_bin_[x_chord_pt / 4][i_y_upper / 4];
              if (raw_img[x_chord_pt][i_y_upper] > spec_lim || pt_edge > edge_lim) {
                masked_raw_bin_[x_chord_pt][i_y_upper] = 0xFF;
                u_l_occ_pts[i_ctr]++;
              }
            }
          }
        }
        if (i_y_lower <= height - 1) {
          int y_test = i_y_lower - pupil_center_y;

          if ((x_test * x_test) + (y_test * y_test) > pupil_radius_squared + 20 * 20) {
            if (x_chord_pt > 0 && x_chord_pt < width && i_y_lower > 0 && i_y_lower < height) {
              int pt_val = (int) raw_img[x_chord_pt][i_y_lower];
              if (pt_val > spec_lim) {
                pt_val = spec_def_val;
              }

              lower_left_gray[i_ctr] = lower_left_gray[i_ctr] + pt_val;
              l_l_g_pts[i_ctr]++;
              if (l_l_g_stop[i_ctr].X == 0) {
                l_l_g_stop[i_ctr].X = x_chord_pt;
                l_l_g_stop[i_ctr].Y = i_y_lower;
              }
              l_l_g_start[i_ctr].X = x_chord_pt;
              l_l_g_start[i_ctr].Y = i_y_lower;

              int pt_edge = (int) vert_edge_bin_[x_chord_pt / 4][i_y_lower / 4]
                  + (int) edge_raw_bin_[x_chord_pt / 4][i_y_lower / 4];
              if (raw_img[x_chord_pt][i_y_lower] > spec_lim || pt_edge > edge_lim) {
                masked_raw_bin_[x_chord_pt][i_y_lower] = 0xFF;
                l_l_occ_pts[i_ctr]++;
              }
            }
          }
        }
      }
    }
  }

  for (int k_ctr = 0; k_ctr < 100; k_ctr++) {
    if (u_r_g_pts[k_ctr] > 0) {
      upper_right_gray[k_ctr] /= u_r_g_pts[k_ctr];
    }
    if (l_r_g_pts[k_ctr] > 0) {
      lower_right_gray[k_ctr] /= l_r_g_pts[k_ctr];
    }
    if (u_l_g_pts[k_ctr] > 0) {
      upper_left_gray[k_ctr] /= u_l_g_pts[k_ctr];
    }
    if (l_l_g_pts[k_ctr] > 0) {
      lower_left_gray[k_ctr] /= l_l_g_pts[k_ctr];
    }
  }

  int lower_right_reference = 0;
  int lower_left_reference = 0;
  for (int m_ctr = 0; m_ctr < 50; m_ctr++) {
    lower_right_reference += lower_right_gray[m_ctr];
    lower_left_reference += lower_left_gray[m_ctr];
  }
  int lower_right_ref_avg = lower_right_reference / 50;
  int lower_left_ref_avg = lower_left_reference / 50;

  int lower_right_var = 0;
  int lower_left_var = 0;
  for (int m_ctr = 5; m_ctr < 50; m_ctr++) {
    lower_right_var +=
        (lower_right_gray[m_ctr] - lower_right_ref_avg) * (lower_right_gray[m_ctr] - lower_right_ref_avg);
    lower_left_var += (lower_left_gray[m_ctr] - lower_left_ref_avg) * (lower_left_gray[m_ctr] - lower_left_ref_avg);
  }

  int combined_low_ref = 0;
  for (int m_ctr = 5; m_ctr < 45; m_ctr++) {
    combined_low_ref += lower_right_gray[m_ctr] + lower_left_gray[m_ctr];
  }

  int combined_var = 0;
  int combined_low_avg = combined_low_ref / ((40) * 2);
  for (int m_ctr = 5; m_ctr < 45; m_ctr++) {
    combined_var += (lower_right_gray[m_ctr] - combined_low_avg) * (lower_right_gray[m_ctr] - combined_low_avg);
    combined_var += (lower_left_gray[m_ctr] - combined_low_avg) * (lower_left_gray[m_ctr] - combined_low_avg);
  }
  int std_dev_combined = (int) floor(0.5 + sqrt((double) (combined_var) / 90.0));

  int left_avail_pixels = 0;
  int right_avail_pixels = 0;
  int left_total_pixels = 0;
  int right_total_pixels = 0;
  int l_r_block = 0;
  int u_r_block = 0;
  int l_l_block = 0;
  int u_l_block = 0;
  for (int m_ctr = 0; m_ctr < 100; m_ctr++) {
    if (l_r_block < 3) {
      if (abs(lower_right_gray[m_ctr] - lower_right_ref_avg) > (int) ((double) std_dev_combined * z_val)) {
        l_r_block++;
        right_total_pixels += l_r_g_pts[m_ctr];

        for (int a_ctr = l_r_g_start[m_ctr].X; a_ctr < l_r_g_stop[m_ctr].X; a_ctr++) {
          masked_raw_bin_[a_ctr][l_r_g_start[m_ctr].Y] = 0;
        }
      } else {
        l_r_block = 0;
        right_total_pixels += l_r_g_pts[m_ctr];
        right_avail_pixels += l_r_g_pts[m_ctr] - l_r_occ_pts[m_ctr];
      }
    } else {
      right_total_pixels += l_r_g_pts[m_ctr];

      for (int a_ctr = l_r_g_start[m_ctr].X; a_ctr < l_r_g_stop[m_ctr].X; a_ctr++) {
        masked_raw_bin_[a_ctr][l_r_g_start[m_ctr].Y] = 0;
      }
    }

    if (u_r_block < 3) {
      if (abs(upper_right_gray[m_ctr] - lower_right_ref_avg) > (int) ((double) std_dev_combined * z_val)) {
        u_r_block++;
        right_total_pixels += u_r_g_pts[m_ctr];

        for (int a_ctr = u_r_g_start[m_ctr].X; a_ctr < u_r_g_stop[m_ctr].X; a_ctr++) {
          masked_raw_bin_[a_ctr][u_r_g_start[m_ctr].Y] = 0x00;
        }
      } else {
        right_total_pixels += u_r_g_pts[m_ctr];
        right_avail_pixels += u_r_g_pts[m_ctr] - u_r_occ_pts[m_ctr];
        u_r_block = 0;
      }
    } else {
      right_total_pixels += u_r_g_pts[m_ctr];

      for (int a_ctr = u_r_g_start[m_ctr].X; a_ctr < u_r_g_stop[m_ctr].X; a_ctr++) {
        masked_raw_bin_[a_ctr][u_r_g_start[m_ctr].Y] = 0;
      }
    }

    if (l_l_block < 3) {
      if (abs(lower_left_gray[m_ctr] - lower_left_ref_avg) > (int) ((double) std_dev_combined * z_val)) {

        l_l_block++;
        left_total_pixels += l_l_g_pts[m_ctr];

        for (int a_ctr = l_l_g_start[m_ctr].X; a_ctr < l_l_g_stop[m_ctr].X; a_ctr++) {
          masked_raw_bin_[a_ctr][l_l_g_start[m_ctr].Y] = 0;
        }
      } else {
        left_total_pixels += l_l_g_pts[m_ctr];
        left_avail_pixels += l_l_g_pts[m_ctr] - l_l_occ_pts[m_ctr];
        l_l_block = 0;
      }
    } else {
      left_total_pixels += l_l_g_pts[m_ctr];

      for (int a_ctr = l_l_g_start[m_ctr].X; a_ctr < l_l_g_stop[m_ctr].X; a_ctr++) {
        masked_raw_bin_[a_ctr][l_l_g_start[m_ctr].Y] = 0;
      }
    }

    if (u_l_block < 3) {
      if (abs(upper_left_gray[m_ctr] - lower_left_ref_avg) > (int) ((double) std_dev_combined * z_val)) {
        u_l_block++;
        left_total_pixels += u_l_g_pts[m_ctr];

        for (int a_ctr = u_l_g_start[m_ctr].X; a_ctr < u_l_g_stop[m_ctr].X; a_ctr++) {
          masked_raw_bin_[a_ctr][u_l_g_start[m_ctr].Y] = 0x00;
        }
      } else {
        left_total_pixels += u_l_g_pts[m_ctr];
        left_avail_pixels += u_l_g_pts[m_ctr] - u_l_occ_pts[m_ctr];
        u_l_block = 0;
      }
    } else {
      left_total_pixels += u_l_g_pts[m_ctr];

      for (int a_ctr = u_l_g_start[m_ctr].X; a_ctr < u_l_g_stop[m_ctr].X; a_ctr++) {
        masked_raw_bin_[a_ctr][u_l_g_start[m_ctr].Y] = 0;
      }
    }
  }

  int combined_avail_percent = 0;
  if (left_total_pixels + right_total_pixels > 0) {
    combined_avail_percent =
        100.0 * (left_avail_pixels + right_avail_pixels) / (left_total_pixels + right_total_pixels);
  }
  if (combined_avail_percent > 100) {
    combined_avail_percent = 100.0;
  }
  if (combined_avail_percent < 0) {
    combined_avail_percent = 0.0;
  }

  usable_iris_area_percent_ = combined_avail_percent;
}

void MFilter::RawByteFrameToTwoDimArray(const uint8_t *frame_bytes, uint8_t **raw_img, int width, int height) {
  try {
    int counter = 0;
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        int mod_y = counter / width;
        int mod_x = counter % width;
        raw_img[x][y] = frame_bytes[counter];
        counter++;
      }
    }
  } catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
  }
}

void MFilter::CheckMargins(int width, int height, int iris_center_x, int iris_center_y, int iris_diameter) {
  double b_margin_metric_deduct = 0.0;
  double
      b_margin = ((double) (height - iris_center_y) - ((double) iris_diameter / 2.0)) / ((double) iris_diameter / 2.0);
  if (b_margin < 0.20) {
    b_margin_metric_deduct = (1.0 - (5.0 * b_margin)) / 6.0;
  }

  double l_margin_metric_deduct = 0.0;
  double l_margin = ((double) iris_center_x - ((double) iris_diameter / 2.0)) / ((double) iris_diameter / 2.0);
  if (l_margin < 0.60) {
    l_margin_metric_deduct = (3.0 - (5.0 * l_margin)) / 8.0;
  }

  double r_margin_metric_deduct = 0.0;
  double
      r_margin = ((double) (width - iris_center_x) - ((double) iris_diameter / 2.0)) / ((double) iris_diameter / 2.0);
  if (r_margin < 0.20) {
    r_margin_metric_deduct = (3.0 - (5.0 * r_margin)) / 8.0;
  }

  double t_margin_metric_deduct = 0.0;
  double t_margin = ((double) iris_center_y - ((double) iris_diameter / 2.0)) / ((double) iris_diameter / 2.0);
  if (t_margin < 0.20) {
    t_margin_metric_deduct = (1.0 - (5.0 * t_margin)) / 6.0;
  }

  double t_b_margin_metric_total_deduct;
  if (t_margin_metric_deduct < b_margin_metric_deduct) {
    t_b_margin_metric_total_deduct = b_margin_metric_deduct;
  } else {
    t_b_margin_metric_total_deduct = t_margin_metric_deduct;
  }
  if (t_b_margin_metric_total_deduct >= 1.0) {
    t_b_margin_metric_total_deduct = 1.0;
  }

  double l_r_margin_metric_total_deduct;
  if (l_margin_metric_deduct < r_margin_metric_deduct) {
    l_r_margin_metric_total_deduct = r_margin_metric_deduct;
  } else {
    l_r_margin_metric_total_deduct = l_margin_metric_deduct;
  }
  if (l_r_margin_metric_total_deduct >= 1.0) {
    l_r_margin_metric_total_deduct = 1.0;
  }

  double overall_margin_deduct = t_b_margin_metric_total_deduct + l_r_margin_metric_total_deduct;
  if (overall_margin_deduct > 1.0) {
    overall_margin_deduct = 1.0;
  }

  overall_margin_ = 1.0 - overall_margin_deduct;
}

void MFilter::DownSize(uint8_t **raw_img, int width, int height, int downsize_scale, uint8_t **ds_img) {
  int ds_width = width / downsize_scale;
  int ds_height = height / downsize_scale;

  for (int i_x = 0; i_x < ds_width; i_x++) {
    int orig_x = i_x * downsize_scale;

    for (int i_y = 0; i_y < ds_height; i_y++) {
      int orig_y = i_y * downsize_scale;

      int ds_point = 0;
      for (int j_x = 0; j_x < downsize_scale; j_x++) {
        for (int j_y = 0; j_y < downsize_scale; j_y++) {
          ds_point += (int) raw_img[j_x + orig_x][j_y + orig_y];
        }
      }

      ds_img[i_x][i_y] = (uint8_t) (ds_point / (downsize_scale * downsize_scale));
    }
  }
}

void MFilter::EdgeMap(uint8_t **raw_img, int width, int height) {
  for (int i_y = 4; i_y < height - 5; i_y++) {
    for (int k_x = 0; k_x < 4; k_x++) {
      edge_val_bin_[k_x][i_y] = 0;
      edge_val_bin_[width - 5 + k_x][i_y] = 0;

      pupil_edge_bin_[k_x][i_y] = 0;
      pupil_edge_bin_[width - 5 + k_x][i_y] = 0;

      edge_raw_bin_[k_x][i_y] = 0;
      edge_raw_bin_[width - 5 + k_x][i_y] = 0;
    }

    for (int i_x = 4; i_x < width - 4; i_x++) {
      edge_val_bin_[i_x][i_y] = 0;
      pupil_edge_bin_[i_x][i_y] = 0;

      for (int j_x = 0; j_x < 4; j_x++) {
        int img_pt_1 = (int) raw_img[i_x - j_x - 1][i_y];
        if (img_pt_1 > 200) {
          img_pt_1 = 200;
        }

        int img_pt_2 = (int) raw_img[i_x + j_x][i_y];
        if (img_pt_2 > 200) {
          img_pt_2 = 200;
        }

        edge_val_bin_[i_x][i_y] += img_pt_1 - img_pt_2;
      }
      if (edge_val_bin_[i_x][i_y] < 0) {
        edge_val_bin_[i_x][i_y] *= -1;
      }
      if (edge_val_bin_[i_x][i_y] > 140) {
        edge_val_bin_[i_x][i_y] = 140;
      }
      edge_val_bin_[i_x][i_y] /= 3;
      edge_raw_bin_[i_x][i_y] = (uint8_t) edge_val_bin_[i_x][i_y];
    }
  }
}

void MFilter::GetVertEdges(uint8_t **raw_img, int width, int height) {
  for (int i_x = 0; i_x < width; i_x++) {
    for (int k_y = 0; k_y < 4; k_y++) {
      vert_edge_[i_x][k_y] = 0;
      vert_edge_[i_x][height - k_y - 1] = 0;
    }

    for (int j_y = 4; j_y < height - 8; j_y++) {
      int point_val = (int) raw_img[i_x][j_y - 2] + (int) raw_img[i_x][j_y - 1];
      point_val = point_val - (int) raw_img[i_x][j_y + 1] - (int) raw_img[i_x][j_y];
      if (point_val < 0) {
        point_val = -point_val;
      }
      if (point_val > 255) {
        point_val = 255;
      }

      // Version 2.2.1 change
      vert_edge_[i_x][j_y] = (uint8_t) point_val;
      if (point_val > 20) {
        point_val = 20;
      }
      vert_edge_[i_x][j_y] = point_val;
    }
  }
}
