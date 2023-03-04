#include "Hashing.h"
#include "Utils.h"

using namespace DTCC_BUILDER;

TEST_CASE("Hashing")
{
  SECTION("Hash Point2D")
  {
    Point2D p(1, 2);
    info(Hashing::Hex(Hashing::Hash(p)));
  }

  SECTION("Hash Point3D")
  {
    Point3D p(1, 2, 3);
    info(Hashing::Hex(Hashing::Hash(p)));
  }
}

TEST_CASE("ISO 8559-1 to UTF-8")
{
  std::string testStr("G\345ngv\344g");
  REQUIRE(Utils::Iso88591ToUtf8(testStr) == "Gångväg");
}

TEST_CASE("Get Filename from path")
{
  std::string path = "/path/to/fileName.json";
  std::string pathFileOnly = "fileName.json";
  std::string pathNoFile = "/path/to/dir/";

  REQUIRE(Utils::GetFilename(path) == "fileName.json");
  REQUIRE(Utils::GetFilename(path, true) == "fileName");
  REQUIRE(Utils::GetFilename(pathFileOnly) == "fileName.json");
  REQUIRE(Utils::GetFilename(pathFileOnly, true) == "fileName");
  REQUIRE(Utils::GetFilename(pathNoFile) == "dir");
}
