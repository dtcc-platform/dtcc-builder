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

const char *node_types[] = {"null",  "document", "element", "pcdata",
                            "cdata", "comment",  "pi",      "declaration"};

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

    simple_walker walker;
    // doc.traverse(walker);

    if (!includeRoot && !doc.empty())
    {
      pugi::xml_node root;
      root = doc.first_child();
      ParseNode(root, json, walker);
    }
    else
      ParseNode(doc, json, walker);

    std::cout << json.dump(4) << std::endl;

    // doc.traverse(walker);

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

  struct simple_walker : pugi::xml_tree_walker
  {
    bool for_each(pugi::xml_node &node) override
    {
      for (int i = 0; i < depth(); ++i)
        std::cout << "  "; // indentation

      std::cout << node_types[node.type()] << ": name='" << node.name()
                << "', value='" << node.value() << "', attributes={";
      for (pugi::xml_attribute attr : node.attributes())
      {
        std::cout << attr.name() << "='" << attr.value() << "'";
      }
      std::cout << "}" << std::endl;

      return true; // continue traversal
    }
  };

  static void
  ParseNode(pugi::xml_node &node, nlohmann::json &json, simple_walker &walker)
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
          ParseNode(child, jsonChild, walker);
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
        InsertChildAsJson(json, child, grandTextChild, walker);
        continue;
      }
      else
      {
        json[child.name()] = jsonChild;
        ParseNode(child, json[child.name()], walker);
      }
    }
  }

  static void InsertChildAsJson(nlohmann::json &json,
                                pugi::xml_node &node,
                                pugi::xml_node &child,
                                simple_walker &walker)
  {
    /*node.traverse(walker);*/
    node.remove_child(child);
    InsertJsonValue(node.name(), child.value(), json);
    /*node.traverse(walker);*/
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

  static void
  ParseNode2(pugi::xml_node &node, nlohmann::json &json, simple_walker &walker)
  {
    for (pugi::xml_attribute attr : node.attributes())
      InsertJsonValue(attr.name(), attr.value(), json);
    for (pugi::xml_node child : node.children())
    {

      nlohmann::json jsonChild = {};
      /*if (child.first_attribute() == nullptr && !child.empty() &&
          child.first_child().type() == pugi::node_pcdata)
       {
        InsertJsonValue(child.name(), child.first_child().value(), json);
        continue;
      }*/
      std::string childName = child.name();
      if (json.count(child.name()) > 0)
      {
        ParseNode(child, jsonChild, walker);
        if (json[child.name()].is_array())
        {
          json[child.name()].push_back(jsonChild);
        }
        else
          CreateJsonArray(jsonChild, json, child);
        continue;
      }
      if (child.type() == pugi::node_pcdata)
      {
        InsertJsonValue("#content", child.value(), json);
        continue;
      }
      pugi::xml_node grandTextChild = GetTextChild(json, child);
      if (grandTextChild != child && GetNumChildren(child) == 1)
      {
        std::string grandTextChildName = grandTextChild.value();
        InsertChildAsJson(json, child, grandTextChild, walker);
        continue;
      }
      else
      {
        json[child.name()] = jsonChild;
        ParseNode(child, json[child.name()], walker);
      }
    }
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
