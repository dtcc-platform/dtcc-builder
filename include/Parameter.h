// Copyright (C) 2021 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_PARAMETER_H
#define DTCC_PARAMETER_H

#include "Logging.h"

namespace DTCC
{

/// Parameter type
enum class ParameterType
{
  Bool,
  Int,
  Float,
  String
};

/// Parameter is used to store the value of a parameter with
/// dynamic type and value.

class Parameter : public Printable
{
public:
  /// Parameter key
  std::string Key{};

  /// Parameter type
  ParameterType Type{ParameterType::Bool};

  // Number of times the parameter has been accessed
  mutable size_t AccessCount{0};

  /// Create empty parameter (defaults to bool)
  explicit Parameter() {}

  /// Create parameter
  explicit Parameter(ParameterType type, const std::string &key)
      : Type(type), Key(key)
  {
  }

  /// Set bool value
  const Parameter &operator=(bool value)
  {
    if (Type != ParameterType::Bool)
      Error("Unable to set parameter \"" + Key + "\"; not a bool parameter");
    valueBool = value;
  }

  /// Set int value
  const Parameter &operator=(int value)
  {
    if (Type != ParameterType::Int)
      Error("Unable to set parameter \"" + Key + "\"; not an int parameter");
    valueInt = value;
  }

  /// Set float value
  const Parameter &operator=(double value)
  {
    if (Type != ParameterType::Float)
      Error("Unable to set parameter \"" + Key + "\"; not a float parameter");
    valueFloat = value;
  }

  /// Set string value
  const Parameter &operator=(const std::string &value)
  {
    if (Type != ParameterType::String)
      Error("Unable to set parameter \"" + Key + "\"; not a string parameter");
    valueString = value;
  }

  /// Set string value (handle string literals, will otherwise go to bool)
  const Parameter &operator=(const char *value)
  {
    if (Type != ParameterType::String)
      Error("Unable to set parameter \"" + Key + "\"; not a string parameter");
    valueString = std::string(value);
  }

  /// Return bool value
  operator bool() const
  {
    if (Type != ParameterType::Bool)
      Error("Unable to access parameter \"" + Key + "\"; not a bool parameter");
    AccessCount++;
    return valueBool;
  }

  /// Return int value
  operator int() const
  {
    if (Type != ParameterType::Int)
      Error("Unable to access parameter \"" + Key + "\"; not an int parameter");
    AccessCount++;
    return valueInt;
  }

  /// Return unsigned int value
  operator size_t() const
  {
    if (Type != ParameterType::Int)
      Error("Unable to access parameter \"" + Key + "\"; not an int parameter");
    if (valueInt < 0)
      Error("Unable to access parameter \"" + Key +
            "\" as unsigned int; value is negative");
    AccessCount++;
    return static_cast<size_t>(valueInt);
  }

  /// Return float value
  operator double() const
  {
    if (Type != ParameterType::Float)
      Error("Unable to access parameter \"" + Key +
            "\"; not a float parameter");
    AccessCount++;
    return valueFloat;
  }

  /// Return string value
  operator std::string() const
  {
    if (Type != ParameterType::String)
      Error("Unable to access parameter \"" + Key +
            "\"; not a string parameter");
    return valueString;
  }

  /// Return type string
  std::string TypeString() const
  {
    switch (Type)
    {
    case ParameterType::Bool:
      return "Bool";
    case ParameterType::Int:
      return "Int";
    case ParameterType::Float:
      return "Float";
    case ParameterType::String:
      return "String";
    default:
      return "Unknown type";
    }
  }

  /// Return value string
  std::string ValueString() const
  {
    switch (Type)
    {
    case ParameterType::Bool:
      return str(valueBool);
    case ParameterType::Int:
      return str(valueInt);
    case ParameterType::Float:
      return str(valueFloat);
    case ParameterType::String:
      return valueString;
    default:
      return "Unknown value";
    }
  }

  /// Pretty-print
  std::string __str__() const override
  {
    return "Parameter " + Key + " = " + ValueString() + " (" + TypeString() +
           ")";
  }

private:
  // Bool value (optional)
  bool valueBool{};

  // Int value (optional)
  int valueInt{};

  // Float value (optional)
  double valueFloat{};

  // String value (optional)
  std::string valueString{};
};

} // namespace DTCC

#endif
