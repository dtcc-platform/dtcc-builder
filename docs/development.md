# Development

This section contains documentation relevant for developers of DTCC
Builder.

## Coding style

DTCC Builder uses Microsoft C# coding style (for both C++ and Python
code):

```
ClassName
MemberFunction
PublicMemberVariable
privateMemberVariable
argumentName
variableName
```

Code formatting is enforced using
[ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) as defined
in the top-level `.clang-format` file. The style is based on the
default LLVM style with minimal modifications.

Algorithms should be implemented as static functions in separate
classes (for example in the class `MeshGenerator` rather than in the
class `Mesh`). This means that pure data classes (like `Mesh`) can be
kept clean with only class data and functions for data access and
initialization.

## Processing data from Lantmäteriet

DTCC Builder uses the following data sources from Lantmäteriet:

* Point clouds (Lantmäteriet:Laserdata NH 2019 (laz); EPSG:3006)
* Property maps (Lantmäteriet:Fastighetskartan Bebyggelse; EPSG:3006)
* Road network maps (Lantmäteriet:Vägkartan; EPSG:3006)

Chalmers has a license for downloading data from `http://zeus.slu.se`.

Point cloud data come in the form of a number square grids large
enough to cover the requested domain. Each point cloud is compressed
as a RAR file. Uncompress it to get the LAS point cloud file, for
example:

    unrar e 09B008_64050_3225_25.rar

This will create the file 09B008_64050_3225_25.las. The lower left
corner will in this example be at EPSG:3006 coordinates (6405000,
322500).

Property map data and the road network comes in the form of SHP files
(with corresponding SHX, DBF and PRJ files). The files of interest are
the ones named `by_get.*` and `vl_*`.

The unit of length is metres relative to the SWEREF99 TM (EPSG:3006)
coordinate system.

## Tips & tricks

If you are using Windows, you might want to make sure that Git does
not convert Unix-style file endings on checkout. This can be
accomplished by:

    git config --global dtcc-builder.autocrlf false

## Remote development in VS Code

On the left-side menu, go to Remote Explorer, on the SSH line press the + sign and add `username@develop.dtcc.chalmers.se` for the case of user `username` and the host being `develop.dtcc.chalmers.se`

Testing new CI

![image](https://user-images.githubusercontent.com/125367195/231126612-d6031bce-ca2d-4340-b0e5-9e728da57238.png)

Then after `develop.dtcc.chalmers.se` appears in the list, click on the connect to current (right arrow) or new window (plus with folder) signs respectively.

![image](https://user-images.githubusercontent.com/125367195/231126959-d1dc6498-576b-42aa-95ce-aff8df80c110.png)

Then you can use `Open...` for opening a folder/file (eg `/home/username/dtcc-builder`) and `Terminal-> New terminal` to have a new terminal connected.


