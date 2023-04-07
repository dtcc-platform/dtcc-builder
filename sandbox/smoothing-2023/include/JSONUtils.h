// Copyright (C) 2020-2021 ReSpace AB, Anton J Olsson
// Licensed under the MIT License

#include <Logging.h>
#include <nlohmann/json.hpp>

/// Read JSON data from file
static void Read(nlohmann::json &json, std::string fileName)
{
  info("JSON: Reading from file " + fileName + "...");
  std::ifstream f(fileName);
  if (!f)
    error("Unable to read from file " + fileName);
  f >> json;
}

/// Write JSON data to file. indent: number of spaces in indent (-1 = no indent)
static void
Write(const nlohmann::json &json, std::string fileName, int indent = -1)
{
  info("JSON: Writing to file " + fileName + "...");
  std::ofstream f(fileName);
  if (!f)
    error("Unable to write to file " + fileName);
  f << json.dump(indent);
}

/// Read object from file
template <class T> static void Read(T &t, std::string fileName)
{
  nlohmann::json json{};
  Read(json, fileName);
  Deserialize(t, json);
}

/// Write object to file
template <class T>
static void Write(const T &t, std::string fileName, int indent = -1)
{
  nlohmann::json json{};
  Serialize(t, json);
  Write(json, fileName, indent);
}

/// Write object and origin (offset) to file
template <class T>
static void Write(const T &t, std::string fileName, const Point2D origin)
{
  nlohmann::json json{};
  Serialize(t, json, origin);
  Write(json, fileName);
}

/// Deserialize object from string
template <class T> static void ReadString(T &t, const std::string &jsonString)
{
  nlohmann::json json{};
  json = nlohmann::json::parse(jsonString);
  Deserialize(t, json);
}

/// Serialize object to string
template <class T> static void WriteString(const T &t, std::string &jsonString)
{
  nlohmann::json json{};
  Serialize(t, json);
  jsonString = json.dump();
}

/// Get type of JSON file
static void CheckType(std::string typeName, const nlohmann::json &json)
{
  if (json.find("Type") == json.end() || json["Type"] != typeName)
    error("Unable to read JSON data; expecting type " + typeName + ".");
}

/// Check json document for key
static bool HasKey(std::string key, const nlohmann::json &json)
{
  return !(json.find(key) == json.end());
}

/// Try reading boolean value
static bool ToBool(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Try reading integer value
static int ToInt(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Try reading integer value
static size_t ToUnsignedInt(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    error("Missing field '" + key + "' in JSON file.");
  long int value = json[key];
  if (value < 0)
    error("Expecting non-negative integer for field '" + key + "'");
  return static_cast<size_t>(value);
}

/// Try reading double value
static double ToDouble(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Try reading string value
static std::string ToString(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Get JSON object in JSON array by key/value pair in object.
/// \tparam T value's type
/// \param key key in key/value/pair
/// \param value value in key/value/pair
/// \param jsonArray JSON array containing object
/// \return JSON object, if found
template <typename T>
static nlohmann::json GetObjectByAttribute(const std::string key,
                                           T value,
                                           const nlohmann::json &jsonArray)
{
  for (const auto &jsonObj : jsonArray)
    if (jsonObj[key] == value)
      return jsonObj;
  error("No object with key " + key + " and value " + value + " in JSON array");
  return nullptr;
}
