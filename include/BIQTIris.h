// #######################################################################
// NOTICE
//
// This software (or technical data) was produced for the U.S. Government
// under contract, and is subject to the Rights in Data-General Clause
// 52.227-14, Alt. IV (DEC 2007).
//
// Copyright 2022 The MITRE Corporation. All Rights Reserved.
// #######################################################################

#ifndef BIQTIRIS_H
#define BIQTIRIS_H

#include "ImageOps.h"
#include "ProviderInterface.h"
#include "opencv2/highgui/highgui.hpp"
#include <json/json.h>
#include <json/value.h>
#include <map>
#include <string>
#include <utility>

class BIQTIris : public Provider {

 private:
  // Provider object
  MFilter mfo;

 public:
  BIQTIris();

  Provider::EvaluationResult evaluate(const std::string &file) override;
};

#endif