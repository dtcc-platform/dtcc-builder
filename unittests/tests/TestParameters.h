#include "JSON.h"

using namespace DTCC;

TEST_CASE("Parse Parameters")
{
  Parameters p;
  p["DataDirectory"] = "359a538b-4616-4c33-b9b0-3870776fa28a";
  JSON::Write(p, "Parameters.json");
  Parameters q;
  JSON::Read(q, "Parameters.json");
  const std::string d = q["DataDirectory"];
  REQUIRE(d == "359a538b-4616-4c33-b9b0-3870776fa28a");
}