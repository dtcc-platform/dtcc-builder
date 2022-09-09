# Development

This section contains documentation relevant for developers of DTCC
Buiulder.

## Design

### Code organization

The DTCC Platform is organized as a collection of independent but
interoperable components. Each component may be implemented using
different libraries, and languages (C++, Python, ...) but follows a
common naming scheme and provides a standardized command-line
interface.

Common C++ code that is used across components is *header only* and is
placed in the common directory `include`. The common code should have
no (or minimal) external dependencies.

### Coding style

The DTCC Platform uses Microsoft C# coding style (for both C++ and
Python code):

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
