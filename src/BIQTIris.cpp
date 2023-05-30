// #######################################################################
// NOTICE
//
// This software (or technical data) was produced for the U.S. Government
// under contract, and is subject to the Rights in Data-General Clause
// 52.227-14, Alt. IV (DEC 2007).
//
// Copyright 2022 The MITRE Corporation. All Rights Reserved.
// #######################################################################

#include "BIQTIris.h"

/**
 *  Creates a BIQTIris instance
 */
BIQTIris::BIQTIris() {
  // Initialize metadata
  std::string biqt_home = std::string(getenv("BIQT_HOME")) +
      std::string("/providers/BIQTIris/descriptor.json");
  std::ifstream desc_file(biqt_home.c_str(), std::ifstream::binary);
  desc_file >> DescriptorObject;

  // Initialize module
  mfo.Initialize();
}

/**
 * Evaluates the iris images.
 *
 * @param file the input file.
 *
 * @return The result of the evaluation.
 */
Provider::EvaluationResult BIQTIris::evaluate(const std::string &file) {
  // Initialize some variables
  Provider::EvaluationResult eval_result;
  Provider::QualityResult quality_result;
  try {
    // Read the image
    cv::Mat img = cv::imread(file, cv::IMREAD_GRAYSCALE);
    uint8_t *raw_img = img.data;

    if (img.rows == 0 || img.cols == 0 || img.data == nullptr || !img.isContinuous()) {
      std::cerr << "Error when attempting to read '" << file << "' "
                << std::endl;
      eval_result.errorCode = 1;
      return eval_result;
    }

    if (img.cols < mfo.min_width_ || img.rows < mfo.min_height_) {
      // March 2023 minor revisions
	  // Changed from if (img.rows < mfo.min_width_ || img.cols < mfo.min_height_) 
      // BIQT Iris has been observed to segfault on images that are very small. This block is intended
      // to prevent that.
      std::cerr << "File '" << file << "' is too small to process ("
                << img.cols << "x" << img.rows << ")." << std::endl;
      eval_result.errorCode = 1;
      return eval_result;
    }
    if (img.cols > mfo.max_width_ || img.rows > mfo.max_height_) {
      std::cerr << "File '" << file << "' is too large to process ("
                << img.cols << "x" << img.rows << ")." << std::endl;
      eval_result.errorCode = 2;
      return eval_result;
    }

    // Calculate quality scores
    double overall_quality = mfo.GetQualityFromImageFrame(raw_img, img.cols, img.rows);

    // Quality Attributes Map
    quality_result.metrics["quality"] = overall_quality;
    quality_result.metrics["contrast"] = mfo.GetContrastScore();
    quality_result.metrics["sharpness"] = mfo.GetDefocusScore();

    quality_result.metrics["iris_sclera_gs"] = mfo.GetISGSDiffMeanAvg();
    quality_result.metrics["iris_pupil_gs"] = mfo.GetIrisPupilGSDiff();

    quality_result.metrics["pupil_circularity_avg_deviation"] = mfo.GetPupilCircularityDeviationAvg();

    // (Normalized) Quality Attributes Map
    quality_result.metrics["normalized_contrast"] = mfo.GetNContrast();
    quality_result.metrics["normalized_sharpness"] = mfo.GetNDefocus();

    quality_result.metrics["normalized_iris_diameter"] = mfo.GetNIrisID();
    quality_result.metrics["normalized_iris_sclera_gs"] = mfo.GetNISGSMean();
    quality_result.metrics["normalized_iris_pupil_gs"] = mfo.GetNIPGSDiff();
    quality_result.metrics["normalized_iso_usable_iris_area"] = mfo.GetNIrisVis();

    // ISO Metrics Map
    quality_result.metrics["iso_usable_iris_area"] = mfo.GetUsableIrisAreaPercent();
    quality_result.metrics["iso_iris_sclera_contrast"] = mfo.GetISOIrisScleraContrast();
    quality_result.metrics["iso_iris_pupil_contrast"] = mfo.GetISOIrisPupilContrast();
    quality_result.metrics["iso_pupil_boundary_circularity"] = mfo.GetISOPupilBoundaryCircularity();
    quality_result.metrics["iso_greyscale_utilization"] = mfo.GetISOGreyscaleUtilization();
    quality_result.metrics["iso_iris_pupil_ratio"] = mfo.GetISOPIRatio();
    quality_result.metrics["iso_iris_pupil_concentricity"] = mfo.GetISOIPConcentricity();
    quality_result.metrics["iso_margin_adequacy"] = mfo.GetISOMarginAdequacy();
    quality_result.metrics["iso_sharpness"] = mfo.GetISOSharpness();
    // Normalized ISO Metrics Map
    quality_result.metrics["normalized_iso_sharpness"] = mfo.GetNormalizedISOSharpness();
    quality_result.metrics["normalized_iso_iris_pupil_ratio"] = mfo.GetNormalizedISOPIRatio();
    quality_result.metrics["normalized_iso_iris_pupil_contrast"] = mfo.GetNormalizedISOIrisPupilContrast();
    quality_result.metrics["normalized_iso_iris_sclera_contrast"] = mfo.GetNormalizedISOIrisScleraContrast();
    quality_result.metrics["normalized_iso_margin_adequacy"] = mfo.GetNormalizedISOMarginAdequacy();
    quality_result.metrics["normalized_iso_greyscale_utilization"] = mfo.GetNormalizedISOGreyscaleUtilization();
    quality_result.metrics["normalized_iso_iris_pupil_concentricity"] = mfo.GetNormalizedISOIPConcentricity();
    quality_result.metrics["normalized_iso_iris_diameter"] = mfo.GetNormalizedISOIrisDiameter();
    // Overall ISO quality metric
    quality_result.metrics["iso_overall_quality"] = mfo.GetIsoOverallQuality();

    // Features Map
    quality_result.features["image_width"] = img.cols;
    quality_result.features["image_height"] = img.rows;
    quality_result.features["iris_center_x"] = mfo.GetIrisCenterX();
    quality_result.features["iris_center_y"] = mfo.GetIrisCenterY();
    quality_result.features["iris_diameter"] = mfo.GetIrisRadius() * 2;
    quality_result.features["pupil_center_x"] = mfo.GetPupilCenterX();
    quality_result.features["pupil_center_y"] = mfo.GetPupilCenterY();
    quality_result.features["pupil_diameter"] = mfo.GetPupilRadius() * 2;
    quality_result.features["pupil_radius"] = mfo.GetPupilRadius();

    // Concatenate certain metrics together (as a string) and put it in the eval_result.
    std::string delim = ",";
    std::string summary;
    summary.append(file);
    summary.append(delim);

    std::size_t found = file.find_last_of("/\\");
    summary.append(file.substr(found + 1)).append(delim);

    std::string met_val = std::to_string(mfo.GetIrisCenterX());
    summary.append(met_val).append(delim);

    met_val = std::to_string(mfo.GetIrisCenterY());
    summary.append(met_val).append(delim);

    met_val = std::to_string(mfo.GetIrisRadius());
    summary.append(met_val).append(delim);

    met_val = std::to_string(mfo.GetPupilCenterX());
    summary.append(met_val).append(delim);

    met_val = std::to_string(mfo.GetPupilCenterY());
    summary.append(met_val).append(delim);

    met_val = std::to_string(mfo.GetPupilRadius());
    summary.append(met_val);

    eval_result.message = summary;
    eval_result.errorCode = (overall_quality < 0) ? 1 : 0;
    eval_result.errorCode = (overall_quality < 0) ? 1 : 0;
    eval_result.qualityResult.push_back(std::move(quality_result));
    return eval_result;
  }
  catch (std::exception &ex) {
    std::cerr << "exception caught: " << ex.what() << std::endl;
    std::cerr << "File '" << file << "' could not be read." << std::endl;
    eval_result.errorCode = 3;
    return eval_result;
  }
}

/**
 * Runs BIQTIris with the given parameters.
 * @param cFilePath The path to the input file.
 * @return the result status.
 */
DLL_EXPORT const char *provider_eval(const char *cFilePath) {
  BIQTIris p;

  // Calculate result
  std::string filePath(cFilePath);
  Provider::EvaluationResult result = p.evaluate(filePath);
  return Provider::serializeResult(result);
}
