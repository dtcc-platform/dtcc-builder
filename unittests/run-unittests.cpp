// Copyright (C) 2020 Anders Logg

#define CATCH_CONFIG_MAIN

#include "catch.hpp"

//using namespace DTCC;

TEST_CASE("Foo")
{
  SECTION("Bar")
  {
    double x = 1.0;
    REQUIRE(x == Approx(1.0));
  }
}
