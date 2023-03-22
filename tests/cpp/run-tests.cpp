// Copyright (C) 2020-2021 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#define CATCH_CONFIG_MAIN

#include "catch.hpp"

const std::string RootPath{"./"};

// Run these tests. You can comment out any test you don't want run
#include "tests/TestBoundingBox.h"
#include "tests/TestBuilding.h"
#include "tests/TestCityModel.h"
// #include "tests/TestDatamodel.h"
#include "tests/TestGEOS.h"
#include "tests/TestGrid.h"
#include "tests/TestPointcloud.h"
#include "tests/TestPolygon.h"
#include "tests/TestProtobufConversion.h"
#include "tests/TestRoofSegmentation.h"
// #include "tests/TestUUID.h"
#include "tests/TestUtils.h"
