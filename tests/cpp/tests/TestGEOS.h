#include "GEOS.h"
#include "model/Polygon.h"

using namespace DTCC_BUILDER;

// Check sum for comparing against sandbox/MergePolygonsGEOS.py
double CheckSum(const Polygon &polygon)
{
  double s = 0.0;
  for (const auto &p : polygon.Vertices)
    s += p.x + p.y;
  std::cout << "Checksum: " << s << std::endl;
  return s;
}

TEST_CASE("GEOS::MergePolygons")
{
  const double TOL = 0.1;
  GEOS::Init();

  // Note: These are the same tests as in the sandbox scripts
  // MergePolygons.py and MergePolygonsGEOS.py.

  SECTION("Test case 0")
  {
    std::cout << "Test case 0" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    for (const auto &p : p0.Vertices)
      p1.Vertices.push_back(p + Vector2D(1.1, 0.5));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(14.3));
  }

  SECTION("Test case 1")
  {
    std::cout << "Test case 1" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    for (const auto &p : p0.Vertices)
      p1.Vertices.push_back(p + Vector2D(0.6, 0.5));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(12.4));
  }

  SECTION("Test case 2")
  {
    std::cout << "Test case 2" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    p1.Vertices.push_back(Point2D(1.1, 0.5));
    p1.Vertices.push_back(Point2D(1.5, 0));
    p1.Vertices.push_back(Point2D(2, 0.5));
    p1.Vertices.push_back(Point2D(1.5, 1));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(14.109756));
  }

  SECTION("Test case 3")
  {
    std::cout << "Test case 3" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    for (const auto &p : p0.Vertices)
      p1.Vertices.push_back(p + Vector2D(1.0, 0.0));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(6.0));
  }

  SECTION("Test case 4")
  {
    std::cout << "Test case 4" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(551.02099997, 57.5619951));
    p0.Vertices.push_back(Point2D(557.94999997, 184.41399511));
    p0.Vertices.push_back(Point2D(545.64399997, 185.08599511));
    p0.Vertices.push_back(Point2D(539.38399997, 70.5609951));
    p0.Vertices.push_back(Point2D(530.07099997, 71.0789951));
    p0.Vertices.push_back(Point2D(529.47099997, 58.7409951));

    Polygon p1;
    p1.Vertices.push_back(Point2D(529.47099997, 58.7409951));
    p1.Vertices.push_back(Point2D(530.07099997, 71.0789951));
    p1.Vertices.push_back(Point2D(460.38899997, 74.87799509));
    p1.Vertices.push_back(Point2D(462.36899997, 111.58199508));
    p1.Vertices.push_back(Point2D(449.93899997, 112.26099508));
    p1.Vertices.push_back(Point2D(447.26099997, 63.22899508));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(4873.52996));
  }

  SECTION("Test case 5")
  {
    std::cout << "Test case 5" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(0, 0));
    p0.Vertices.push_back(Point2D(1, 0));
    p0.Vertices.push_back(Point2D(1, 1));
    p0.Vertices.push_back(Point2D(0, 1));

    Polygon p1;
    for (const auto &p : p0.Vertices)
      p1.Vertices.push_back(p + Vector2D(1.1, 0.5));
    p1.Vertices[0] -= Vector2D(0.1, 0);

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(14.344555));
  }

  SECTION("Test case 6")
  {
    std::cout << "Test case 6" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(338.82099997, 326.20099506));
    p0.Vertices.push_back(Point2D(335.96199997, 268.23299506));
    p0.Vertices.push_back(Point2D(346.62299997, 267.70599506));
    p0.Vertices.push_back(Point2D(346.53899997, 263.96099506));
    p0.Vertices.push_back(Point2D(354.08399997, 263.68399506));
    p0.Vertices.push_back(Point2D(357.65099997, 328.92299506));
    p0.Vertices.push_back(Point2D(350.08599997, 329.20599506));
    p0.Vertices.push_back(Point2D(349.69199997, 325.66799506));

    Polygon p1;
    p1.Vertices.push_back(Point2D(291.75699997, 346.98899504));
    p1.Vertices.push_back(Point2D(295.37799997, 346.74599504));
    p1.Vertices.push_back(Point2D(295.12199997, 342.06399505));
    p1.Vertices.push_back(Point2D(331.63597417, 340.09247421));
    p1.Vertices.push_back(Point2D(331.30799997, 333.47099505));
    p1.Vertices.push_back(Point2D(334.88199997, 329.60499505));
    p1.Vertices.push_back(Point2D(334.81099997, 327.91199505));
    p1.Vertices.push_back(Point2D(336.40799997, 325.82299505));
    p1.Vertices.push_back(Point2D(338.79999997, 325.71299506));
    p1.Vertices.push_back(Point2D(338.81999997, 326.14699506));
    p1.Vertices.push_back(Point2D(349.68999997, 325.61399506));
    p1.Vertices.push_back(Point2D(350.08399997, 329.20399506));
    p1.Vertices.push_back(Point2D(350.85899997, 329.17599506));
    p1.Vertices.push_back(Point2D(350.91399997, 330.10899506));
    p1.Vertices.push_back(Point2D(347.11499997, 334.38699506));
    p1.Vertices.push_back(Point2D(347.46419003, 339.23785189));
    p1.Vertices.push_back(Point2D(347.46599997, 339.26299506));
    p1.Vertices.push_back(Point2D(346.98131866, 339.28840791));
    p1.Vertices.push_back(Point2D(347.79299997, 354.31499506));
    p1.Vertices.push_back(Point2D(330.95999188, 355.22284556));
    p1.Vertices.push_back(Point2D(330.95999997, 355.22299505));
    p1.Vertices.push_back(Point2D(331.11599997, 358.10699505));
    p1.Vertices.push_back(Point2D(325.27699997, 358.42299505));
    p1.Vertices.push_back(Point2D(325.23293055, 357.60827579));
    p1.Vertices.push_back(Point2D(318.91799997, 357.95799505));
    p1.Vertices.push_back(Point2D(318.86933569, 357.14365328));
    p1.Vertices.push_back(Point2D(318.86299997, 357.14399505));
    p1.Vertices.push_back(Point2D(308.43015642, 357.70677084));
    p1.Vertices.push_back(Point2D(308.42599997, 357.70699505));
    p1.Vertices.push_back(Point2D(308.42132643, 357.62019097));
    p1.Vertices.push_back(Point2D(292.37499997, 358.47599504));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(17905.235323));
  }

  SECTION("Test case 7")
  {
    std::cout << "Test case 7" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(689.40299996, 59.61399514));
    p0.Vertices.push_back(Point2D(689.35699996, 64.63399514));
    p0.Vertices.push_back(Point2D(676.57299996, 65.25199513));
    p0.Vertices.push_back(Point2D(676.74399996, 59.23399514));

    Polygon p1;
    p1.Vertices.push_back(Point2D(676.56899996, 65.37299513));
    p1.Vertices.push_back(Point2D(676.57299996, 65.25199513));
    p1.Vertices.push_back(Point2D(689.35699996, 64.63399514));
    p1.Vertices.push_back(Point2D(689.35399996, 64.92899514));
    p1.Vertices.push_back(Point2D(689.96199996, 69.64799514));
    p1.Vertices.push_back(Point2D(677.26999996, 71.23799513));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(4489.04597));
  }

  SECTION("Test case 8")
  {
    std::cout << "Test case 8" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(83.58099997, 249.880995));
    p0.Vertices.push_back(Point2D(88.10899997, 249.639995));
    p0.Vertices.push_back(Point2D(88.28699997, 252.991995));
    p0.Vertices.push_back(Point2D(83.75999997, 253.231995));

    Polygon p1;
    p1.Vertices.push_back(Point2D(71.24799997, 250.721995));
    p1.Vertices.push_back(Point2D(71.48399997, 254.38399499));
    p1.Vertices.push_back(Point2D(75.36099997, 254.137995));
    p1.Vertices.push_back(Point2D(75.12499997, 250.480995));
    p1.Vertices.push_back(Point2D(83.93599997, 256.027995));
    p1.Vertices.push_back(Point2D(64.79599997, 257.07299499));
    p1.Vertices.push_back(Point2D(64.44799997, 250.89299499));
    p1.Vertices.push_back(Point2D(83.59299997, 249.847995));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(2326.688678));
  }

  SECTION("Test case 9")
  {
    std::cout << "Test case 9" << std::endl;
    std::cout << "-----------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(539.39399997, 411.2619951));
    p0.Vertices.push_back(Point2D(534.37199997, 411.5369951));
    p0.Vertices.push_back(Point2D(533.99599997, 404.7629951));
    p0.Vertices.push_back(Point2D(539.02399997, 404.4879951));

    Polygon p1;
    p1.Vertices.push_back(Point2D(526.26499997, 393.8649951));
    p1.Vertices.push_back(Point2D(553.67099997, 392.2179951));
    p1.Vertices.push_back(Point2D(554.30599997, 403.8439951));
    p1.Vertices.push_back(Point2D(533.99599997, 404.7629951));
    p1.Vertices.push_back(Point2D(534.94199997, 421.9599951));
    p1.Vertices.push_back(Point2D(523.44399997, 422.4769951));
    p1.Vertices.push_back(Point2D(522.29099997, 404.43099509));
    p1.Vertices.push_back(Point2D(526.81199997, 404.2289951));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(9424.87874));
  }

  SECTION("Test case 10")
  {
    std::cout << "Test case 10" << std::endl;
    std::cout << "------------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(526.26499997, 393.8649951));
    p0.Vertices.push_back(Point2D(553.67099997, 392.2179951));
    p0.Vertices.push_back(Point2D(554.30599997, 403.8439951));
    p0.Vertices.push_back(Point2D(539.02658749, 404.53536783));
    p0.Vertices.push_back(Point2D(539.39399997, 411.2619951));
    p0.Vertices.push_back(Point2D(534.37199997, 411.5369951));
    p0.Vertices.push_back(Point2D(534.36864511, 411.53717965));
    p0.Vertices.push_back(Point2D(534.94199997, 421.9599951));
    p0.Vertices.push_back(Point2D(523.44399997, 422.4769951));
    p0.Vertices.push_back(Point2D(522.29099997, 404.43099509));
    p0.Vertices.push_back(Point2D(526.81199997, 404.2289951));

    Polygon p1;
    p1.Vertices.push_back(Point2D(523.16699997, 418.0989951));
    p1.Vertices.push_back(Point2D(521.35499997, 418.21899509));
    p1.Vertices.push_back(Point2D(520.64999997, 407.1569951));
    p1.Vertices.push_back(Point2D(522.46099997, 407.0369951));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(13163.017939));
  }

  SECTION("Test case 11")
  {
    std::cout << "Test case 11" << std::endl;
    std::cout << "------------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(568.54499997, 407.83799511));
    p0.Vertices.push_back(Point2D(554.56599997, 408.6009951));
    p0.Vertices.push_back(Point2D(552.50099997, 370.9129951));
    p0.Vertices.push_back(Point2D(566.47899997, 370.14999511));

    Polygon p1;
    p1.Vertices.push_back(Point2D(526.26499997, 393.8649951));
    p1.Vertices.push_back(Point2D(553.67099997, 392.2179951));
    p1.Vertices.push_back(Point2D(554.30599997, 403.8439951));
    p1.Vertices.push_back(Point2D(539.02658749, 404.53536783));
    p1.Vertices.push_back(Point2D(539.39399997, 411.2619951));
    p1.Vertices.push_back(Point2D(534.37199997, 411.5369951));
    p1.Vertices.push_back(Point2D(534.36864511, 411.53717965));
    p1.Vertices.push_back(Point2D(534.94199997, 421.9599951));
    p1.Vertices.push_back(Point2D(523.44399997, 422.4769951));
    p1.Vertices.push_back(Point2D(523.16429105, 418.09917449));
    p1.Vertices.push_back(Point2D(521.35499997, 418.21899509));
    p1.Vertices.push_back(Point2D(520.64999997, 407.1569951));
    p1.Vertices.push_back(Point2D(522.45751798, 407.03722582));
    p1.Vertices.push_back(Point2D(522.29099997, 404.43099509));
    p1.Vertices.push_back(Point2D(526.81199997, 404.2289951));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(16962.607815));
  }

  SECTION("Test case 12")
  {
    std::cout << "Test case 12" << std::endl;
    std::cout << "------------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(552.50099997, 370.9129951));
    p0.Vertices.push_back(Point2D(566.47899997, 370.14999511));
    p0.Vertices.push_back(Point2D(568.54499997, 407.83799511));
    p0.Vertices.push_back(Point2D(554.56599997, 408.6009951));
    p0.Vertices.push_back(Point2D(554.30535647, 403.84403036));

    Polygon p1;
    p1.Vertices.push_back(Point2D(554.42399997, 370.7089951));
    p1.Vertices.push_back(Point2D(554.40499997, 368.7879951));
    p1.Vertices.push_back(Point2D(567.26199997, 368.04799511));
    p1.Vertices.push_back(Point2D(567.85999997, 380.42599511));
    p1.Vertices.push_back(Point2D(567.16399997, 380.36599511));
    p1.Vertices.push_back(Point2D(567.03499997, 380.26199511));
    p1.Vertices.push_back(Point2D(566.47899997, 370.14499511));
    p1.Vertices.push_back(Point2D(554.54399997, 370.8039951));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(7542.28551));
  }

  SECTION("Test case 13")
  {
    std::cout << "Test case 13" << std::endl;
    std::cout << "------------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(292.37499997, 358.47599504));
    p0.Vertices.push_back(Point2D(291.75699997, 346.98899504));
    p0.Vertices.push_back(Point2D(295.37799997, 346.74599504));
    p0.Vertices.push_back(Point2D(295.12199997, 342.06399505));
    p0.Vertices.push_back(Point2D(346.97999997, 339.26399506));
    p0.Vertices.push_back(Point2D(347.79299997, 354.31499506));
    p0.Vertices.push_back(Point2D(330.95999997, 355.22299505));
    p0.Vertices.push_back(Point2D(331.11599997, 358.10699505));
    p0.Vertices.push_back(Point2D(325.27699997, 358.42299505));
    p0.Vertices.push_back(Point2D(325.15126955, 356.09859427));
    p0.Vertices.push_back(Point2D(325.12099997, 355.53899505));
    p0.Vertices.push_back(Point2D(318.79399997, 355.88299505));
    p0.Vertices.push_back(Point2D(308.35499997, 356.44199505));
    p0.Vertices.push_back(Point2D(308.42499997, 357.61999505));

    Polygon p1;
    p1.Vertices.push_back(Point2D(338.79999997, 325.71299506));
    p1.Vertices.push_back(Point2D(338.81999997, 326.14699506));
    p1.Vertices.push_back(Point2D(349.68999997, 325.61399506));
    p1.Vertices.push_back(Point2D(350.08399997, 329.20399506));
    p1.Vertices.push_back(Point2D(350.85899997, 329.17599506));
    p1.Vertices.push_back(Point2D(350.91399997, 330.10899506));
    p1.Vertices.push_back(Point2D(347.11499997, 334.38699506));
    p1.Vertices.push_back(Point2D(347.46599997, 339.26299506));
    p1.Vertices.push_back(Point2D(331.63599997, 340.09299505));
    p1.Vertices.push_back(Point2D(331.30799997, 333.47099505));
    p1.Vertices.push_back(Point2D(334.88199997, 329.60499505));
    p1.Vertices.push_back(Point2D(334.81099997, 327.91199505));
    p1.Vertices.push_back(Point2D(336.40799997, 325.82299505));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
    REQUIRE(CheckSum(pm) == Approx(16767.20206));
  }

  SECTION("Test case 14")
  {
    std::cout << "Test case 14" << std::endl;
    std::cout << "------------" << std::endl;

    Polygon p0;
    p0.Vertices.push_back(Point2D(231.14399994979613, 150.44799979962409));
    p0.Vertices.push_back(Point2D(230.8629999498953, 156.68399979919195));
    p0.Vertices.push_back(Point2D(227.78699994989438, 156.25299980025738));
    p0.Vertices.push_back(Point2D(228.67399994982406, 150.07099980022758));
    p0.Vertices.push_back(Point2D(231.14399994979613, 150.44799979962409));

    Polygon p1;
    p1.Vertices.push_back(Point2D(207.7159999498399, 139.2019997993484));
    p1.Vertices.push_back(Point2D(217.76999994984362, 140.09099979978055));
    p1.Vertices.push_back(Point2D(216.67799994983943, 151.3609997993335));
    p1.Vertices.push_back(Point2D(214.94699994981056, 154.25899980030954));
    p1.Vertices.push_back(Point2D(215.61399995000102, 155.51899980101734));
    p1.Vertices.push_back(Point2D(227.78699994989438, 156.25299980025738));
    p1.Vertices.push_back(Point2D(226.24299994989997, 166.89199980068952));
    p1.Vertices.push_back(Point2D(207.03799994988367, 164.323999799788));
    p1.Vertices.push_back(Point2D(205.54999994981335, 162.44599980022758));
    p1.Vertices.push_back(Point2D(207.7159999498399, 139.2019997993484));

    Polygon pm = GEOS::MergePolygons(p0, p1, TOL);
  }
}
