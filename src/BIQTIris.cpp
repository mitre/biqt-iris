// #######################################################################
// NOTICE
//
// This software (or technical data) was produced for the U.S. Government
// under contract, and is subject to the Rights in Data-General Clause
// 52.227-14, Alt. IV (DEC 2007).
//
// Copyright 2019 The MITRE Corporation. All Rights Reserved.
// #######################################################################

#include <fstream>
#include <iostream>
#include <map>

#include "BIQTIris.h"

/**
 *  Creates a BIQTIris instance
 */
BIQTIris::BIQTIris()
{
    // Initialize metadata
    std::string biqt_home = getenv("BIQT_HOME");
    std::ifstream desc_file(biqt_home + "/providers/BIQTIris/descriptor.json",
                            std::ifstream::binary);
    desc_file >> DescriptorObject;

    // Initialize module
    mfo.initialize();
}

/**
 * Evaluates the iris images.
 *
 * @param file the input file.
 *
 * @return The result of the evaluation.
 */
Provider::EvaluationResult BIQTIris::evaluate(const std::string &file)
{
    // Initialize some variables
    Provider::EvaluationResult eval_result;
    Provider::QualityResult quality_result;

    // Read the image
    cv::Mat img = cv::imread(file, 64);
    unsigned char *rawimage = img.data;

    if ((img.rows == 0) || (img.cols == 0) || img.data == nullptr ||
        !img.isContinuous()) {
        std::cerr << "Error when attempting to read '" << file << "' "
                  << std::endl;
        eval_result.errorCode = 1;
        return eval_result;
    }

    // Calculate quality scores
    if (quality_result.metrics["fast_quality"] =
            mfo.fastQuality(rawimage, (int)img.cols, (int)img.rows) < 0) {
        return eval_result;
    }

    // Metrics Map
    double overallQuality =
        mfo.getQualityFromImageFrame(rawimage, img.cols, img.rows);
    quality_result.metrics["quality"] = overallQuality;
    quality_result.metrics["sharpness"] = mfo.getDefocusScore();
    quality_result.metrics["normalized_sharpness"] = mfo.getNDefocus();
    quality_result.metrics["contrast"] = mfo.getContrastScore();
    quality_result.metrics["normalized_contrast"] = mfo.getNContrast();
    quality_result.metrics["iris_sclera_gs"] = mfo.getISGSDiffMeanAvg();
    quality_result.metrics["normalized_iris_sclera_gs"] = mfo.getNISGSMean();
    quality_result.metrics["iris_pupil_gs"] = mfo.getIrisPupilGSDiff();
    quality_result.metrics["normalized_iris_pupil_gs"] = mfo.getNIPGSDiff();
    quality_result.metrics["iris_diameter"] = mfo.getIrisRadius() * 2;
    quality_result.metrics["normalized_iris_diameter"] = mfo.getNIrisID();
    quality_result.metrics["normalized_percent_visible_iris"] =
        mfo.getNIrisVis();
    quality_result.metrics["percent_visible_iris"] =
        mfo.getUsableIrisAreaPercent();

    // Features Map
    quality_result.features["image_width"] = img.cols;
    quality_result.features["image_height"] = img.rows;
    quality_result.features["iris_center_x"] = mfo.getIrisCenterX();
    quality_result.features["iris_center_y"] = mfo.getIrisCenterY();

    eval_result.qualityResult.push_back(std::move(quality_result));

    if (overallQuality < 0) {
        eval_result.errorCode = 1;
    }
    else {
        eval_result.errorCode = 0;
    }
    return eval_result;
}

/**
 * Runs BIQTIris with the given parameters.
 * @param cFilePath The path to the input file.
 * @return the result status.
 */
DLL_EXPORT const char *provider_eval(const char *cFilePath)
{
    BIQTIris p;

    // Calculate result
    std::string filePath(cFilePath);
    Provider::EvaluationResult result = p.evaluate(filePath);
    return Provider::serializeResult(result);
}
