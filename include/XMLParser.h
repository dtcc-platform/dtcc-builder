// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_XMLPARSER_H
#define CORE_XMLPARSER_H

#include "Logging.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include <pugixml.hpp>
#include <string>

namespace DTCC
{

class XMLParser
{

public:
  static nlohmann::json GetJsonFromXML(const char *filePath,
                                       bool includeRoot = false)
  {
    nlohmann::json json;

    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file(filePath);
    if (!result)
      PrintError(result, doc, filePath);

    if (!includeRoot && !doc.empty())
    {
      pugi::xml_node root;
      root = doc.first_child();
      ParseNode(root, json);
    }
    else
      ParseNode(doc, json);

    std::cout << json.dump(4) << std::endl;

    return json;
  }

private:
  static int GetNumChildren(pugi::xml_node &node)
  {
    size_t n = 0;
    for (pugi::xml_node_iterator it = node.begin(); it != node.end(); ++it)
      n++;
    return n;
  }

  static void ParseNode(pugi::xml_node &node, nlohmann::json &json)
  {
    for (pugi::xml_attribute attr : node.attributes())
      InsertJsonValue(attr.name(), attr.value(), json);
    std::map<std::string, int> childNames;
    for (pugi::xml_node child : node.children())
    {
      if (childNames.count(child.name()) == 0)
        childNames.insert({child.name(), 1});
      else
        childNames[child.name()]++;
    }
    for (pugi::xml_node child : node.children())
    {
      nlohmann::json jsonChild = {};
      if (childNames[child.name()] > 1)
      {
        bool textChild = (GetNumChildren(child) == 1);
        if (!textChild)
          ParseNode(child, jsonChild);
        if (json[child.name()].is_array())
        {
          json[child.name()].push_back(textChild ? child.first_child().value()
                                                 : jsonChild);
        }
        else
          json[child.name()] = {textChild ? child.first_child().value()
                                          : jsonChild};
        continue;
      }
      if (child.type() == pugi::node_pcdata)
      {
        InsertJsonValue("#content", child.value(), json);
        continue;
      }
      pugi::xml_node grandTextChild = GetTextChild(json, child);
      if (grandTextChild != child && GetNumChildren(child) == 1 &&
          child.first_attribute() == nullptr)
      {
        std::string grandTextChildName = grandTextChild.value();
        InsertChildAsJson(json, child, grandTextChild);
        continue;
      }
      else
      {
        json[child.name()] = jsonChild;
        ParseNode(child, json[child.name()]);
      }
    }
  }

  static void InsertChildAsJson(nlohmann::json &json,
                                pugi::xml_node &node,
                                pugi::xml_node &child)
  {
    node.remove_child(child);
    InsertJsonValue(node.name(), child.value(), json);
  }

  static pugi::xml_node GetTextChild(nlohmann::json &json, pugi::xml_node &node)
  {
    for (pugi::xml_node child : node.children())
    {
      if (child.type() == pugi::node_pcdata)
      {
        return child;
      }
    }
    return node;
  }

  static void CreateJsonArray(nlohmann::json &jsonChild,
                              nlohmann::json &json,
                              pugi::xml_node &node)
  {
    std::string jsonString = json[node.name()].dump();
    nlohmann::json sameNameChild = nlohmann::json::parse(jsonString);
    json.erase(node.name());
    json[node.name()] = {jsonChild, sameNameChild};
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
    if (value.size() > 1 && value[0] == '0' && value[1] != '.')
      isNumeric = false;
    else
      isNumeric =
          value.find_first_not_of("0123456789.", pos) == std::string::npos &&
          numPoints < 2;

    hasFraction = numPoints > 0;
  }

  static void StripWhitespace(std::string &value)
  {
    value.erase(0, value.find_first_not_of(" \n\r\t"));
    value.erase(value.find_last_not_of(" \n\r\t") + 1);
  }

  static void
  InsertJsonValue(const char *key, std::string value, nlohmann::json &json)
  {
    StripWhitespace(value);
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
