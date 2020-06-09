// Copyright (C) 2020 ReSpace AB
// Licensed under the MIT License

/// Read JSON data from file
static void Read(nlohmann::json& json, std::string fileName)
{
  Info("JSON: Reading from file " + fileName + "...");
  std::ifstream f(fileName);
  if (!f)
    Error("Unable to read from file " + fileName);
  f >> json;
}

/// Write JSON data to file
static void Write(const nlohmann::json& json, std::string fileName)
{
  Info("JSON: Writing to file " + fileName + "...");
  std::ofstream f(fileName);
  if (!f)
    Error("Unable to write to file " + fileName);
  f << json;
}

/// Read object from file
template <class T>
static void Read(T& t, std::string fileName)
{
  nlohmann::json json{};
  Read(json, fileName);
  Deserialize(t, json);
}

/// Write object to file
template <class T>
static void Write(const T& t, std::string fileName)
{
  nlohmann::json json{};
  Serialize(t, json);
  Write(json, fileName);
}

/// Deserialize object from string
template <class T>
static void ReadString(T& t, const std::string& jsonString)
{
  nlohmann::json json{};
  json = nlohmann::json::parse(jsonString);
  Deserialize(t, json);
}

/// Serialize object to string
template <class T> static
void WriteString(const T& t, std::string& jsonString)
{
  nlohmann::json json{};
  Serialize(t, json);
  jsonString = json.dump();
}

/// Get type of JSON file
static void CheckType(std::string typeName, const nlohmann::json& json)
{
  if (json.find("Type") == json.end() || json["Type"] != typeName)
    Error("Unable to read JSON data; expecting type " + typeName + ".");
}

/// Try reading boolean value
static bool ToBool(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    Error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Try reading integer value
static int ToInt(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    Error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Try reading integer value
static size_t ToUnsignedInt(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    Error("Missing field '" + key + "' in JSON file.");
  long int value = json[key];
  if (value < 0)
    Error("Expecting non-negative integer for field '" + key + "'");
  return static_cast<size_t>(value);
}

/// Try reading double value
static double ToDouble(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    Error("Missing field '" + key + "' in JSON file.");
  return json[key];
}

/// Try reading string value
static std::string ToString(std::string key, const nlohmann::json &json)
{
  if (json.find(key) == json.end())
    Error("Missing field '" + key + "' in JSON file.");
  return json[key];
}