//
// Created by Anton J Olsson on 2020-09-25.
//

#ifndef CORE_XMLPARSER_H
#define CORE_XMLPARSER_H

#include "Logging.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include <pugixml.hpp>
#include <string>

namespace DTCC
{

const char *node_types[] = {"null",  "document", "element", "pcdata",
                            "cdata", "comment",  "pi",      "declaration"};

class XMLParser

{

public:
  static nlohmann::json GetJsonFromXML(const char *filePath)
  {

    nlohmann::json json;

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filePath);
    if (!result)
      PrintError(result, doc, filePath);

    ParseNode(doc, json);

    std::cout << json.dump(4) << std::endl;

    return json;
  }

private:

  static void ParseNode(pugi::xml_node &node, nlohmann::json &json)
  {
    for (pugi::xml_attribute attr : node.attributes())
    {
      InsertJsonValue(attr.name(), attr.value(), json);
    }

    for (pugi::xml_node child : node.children())
    {
      nlohmann::json jsonChild = {};
      if (json.count(child.name()) > 0)
      {
        ParseNode(child, jsonChild);
        if (json[child.name()].is_array())
        {
          json[child.name()].push_back(jsonChild);
        }
        else
        {
          std::string jsonString = json[child.name()].dump();
          nlohmann::json sameNameChild = nlohmann::json::parse(jsonString);
          json.erase(child.name());
          json[child.name()] = {jsonChild, sameNameChild};
        }
      }
      else
      {
        json[child.name()] = jsonChild;
        ParseNode(child, json[child.name()]);
      }
    }
  }

  static void MakeNumericAnalysis(const std::string &value,
                                  bool &isNumeric,
                                  bool &isNegative,
                                  bool &hasFraction)
  {
    int pos = value.at(0) == '-' ? 1 : 0;
    isNegative = pos;
    int numPoints = 0;
    size_t i = pos - 1;
    while (++i < value.size())
    {
      if (value[i] == '.')
        numPoints++;
    }
    isNumeric =
        value.find_first_not_of("0123456789.", pos) == std::string::npos &&
        numPoints < 2;
    hasFraction = numPoints > 0;
  }

  static void
  InsertJsonValue(const char *key, const char *value, nlohmann::json &json)
  {
    bool isNumeric, isNegative, hasFraction;
    MakeNumericAnalysis(value, isNumeric, isNegative, hasFraction);
    if (isNumeric)
    {
      if (hasFraction)
        json[key] = std::stod(value);
      else if (isNegative)
        json[key] = std::stoi(value);
      else
        json[key] = (size_t)std::stoul(value);
    }
    else if (std::string(value) == "true")
      json[key] = true;
    else if (std::string(value) == "false")
      json[key] = false;
    else
      json[key] = value;
  }

  static void PrintError(pugi::xml_parse_result result,
                         const pugi::xml_document &doc,
                         const char *path)
  {
    std::stringstream ss;
    ss << "XML [" << path << "] parsed with errors, attr value: ["
       << doc.child("node").attribute("attr").value() << "]\n";
    ss << "Error description: " << result.description() << "\n";
    ss << "Error offset: " << result.offset << " (error at [..."
       << (path + result.offset) << "]\n\n";
    Error(ss.str());
  }
};

} // namespace DTCC

#endif // CORE_XMLPARSER_H
