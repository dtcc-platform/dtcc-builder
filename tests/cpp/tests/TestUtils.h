#include "Hashing.h"
#include "Utils.h"

using namespace DTCC_BUILDER;

TEST_CASE("Hashing")
{
  SECTION("hash Point2D")
  {
    Vector2D p(1, 2);
    info(Hashing::hex(Hashing::hash(p)));
  }

  SECTION("hash Point3D")
  {
    Vector3D p(1, 2, 3);
    info(Hashing::hex(Hashing::hash(p)));
  }
}

TEST_CASE("ISO 8559-1 to UTF-8")
{
  std::string testStr("G\345ngv\344g");
  REQUIRE(Utils::iso88591_to_utf8(testStr) == "Gångväg");
}

TEST_CASE("Get Filename from path")
{
  std::string path = "/path/to/file_name.json";
  std::string pathFileOnly = "file_name.json";
  std::string pathNoFile = "/path/to/dir/";

  REQUIRE(Utils::get_filename(path) == "file_name.json");
  REQUIRE(Utils::get_filename(path, true) == "file_name");
  REQUIRE(Utils::get_filename(pathFileOnly) == "file_name.json");
  REQUIRE(Utils::get_filename(pathFileOnly, true) == "file_name");
  REQUIRE(Utils::get_filename(pathNoFile) == "dir");
}
