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

    return json;
  }

private:
  static int GetNumOfChildren(pugi::xml_node &node)
  {
    size_t n = 0;
    for (pugi::xml_node_iterator it = node.begin(); it != node.end(); ++it)
      n++;
    return n;
  }

  static std::vector<pugi::xml_node> GetTextChildren(pugi::xml_node node)
  {
    std::vector<pugi::xml_node> children;
    for (pugi::xml_node child : node.children())
    {
      if (child.type() == pugi::node_pcdata)
      {
        children.push_back(child);
      }
    }
    return children;
  }

  static void ParseNode(pugi::xml_node &xmlNode, nlohmann::json &json)
  {
    for (pugi::xml_attribute attr : xmlNode.attributes())
      FormatAndInsertJson(attr.name(), attr.value(), json);
    std::map<std::string, int> childNames = GetChildrensNameCount(xmlNode);
    for (pugi::xml_node xmlChild : xmlNode.children())
    {
      nlohmann::json jsonChild = {};
      // Handle multiple child elements with same name
      if (childNames[xmlChild.name()] > 1)
        OnSameNameChildren(jsonChild, json, xmlChild);
      // Handle element with attributes and text child
      else if (xmlChild.type() == pugi::node_pcdata)
        FormatAndInsertJson("#content", xmlChild.value(), json);
      else
      {
        std::vector<pugi::xml_node> grandTextChildren =
            GetTextChildren(xmlChild);
        // Base case of element of unique type and text node as only child
        if (grandTextChildren.size() == 1 && GetNumOfChildren(xmlChild) == 1 &&
            xmlChild.first_attribute() == nullptr)
        {
          InsertAsJsonAndRemove(json, xmlChild, grandTextChildren[0], false);
        }
        else
        {
          // Handle multiple text nodes in element
          if (grandTextChildren.size() > 1)
            PutTextNodesInArray(jsonChild, grandTextChildren, json, xmlChild);
          json[xmlChild.name()] = jsonChild;
          ParseNode(xmlChild, json[xmlChild.name()]);
        }
      }
    }
  }

  static void OnSameNameChildren(nlohmann::json &jsonChild,
                                 nlohmann::json &json,
                                 pugi::xml_node &xmlNode)
  {
    bool hasOnlyTextChildren =
        (GetNumOfChildren(xmlNode) == GetNumTextChildren(xmlNode) &&
         xmlNode.first_attribute() == nullptr);
    if (!hasOnlyTextChildren)
      ParseNode(xmlNode, jsonChild);
    PutNodeInArray(jsonChild, json, xmlNode, hasOnlyTextChildren);
  }

  static void
  PutTextNodesInArray(nlohmann::json &jsonChild,
                      const std::vector<pugi::xml_node> &nodeChildren,
                      nlohmann::json &json,
                      pugi::xml_node &node)
  {
    nlohmann::json arr;
    for (pugi::xml_node grandChild : nodeChildren)
    {
      std::string valueStr(grandChild.value());
      InsertAsJsonAndRemove(arr, node, grandChild, true);
    }
    jsonChild["#content"] = arr;
  }

  static int GetNumTextChildren(pugi::xml_node &node)
  {
    size_t n = 0;
    for (auto &child : node)
      if (child.type() == pugi::node_pcdata)
        n++;
    return n;
  }

  static void PutNodeInArray(nlohmann::json &jsonChild,
                             nlohmann::json &json,
                             pugi::xml_node &xmlNode,
                             bool hasOnlyTextChildren)
  {
    if (json[xmlNode.name()].is_array())
      json[xmlNode.name()].push_back(
          hasOnlyTextChildren ? xmlNode.first_child().value() : jsonChild);
    else
    {
      std::string childName = xmlNode.name();
      json[xmlNode.name()] = {
          hasOnlyTextChildren ? xmlNode.first_child().value() : jsonChild};
    }
  }

  static std::map<std::string, int>
  GetChildrensNameCount(const pugi::xml_node &node)
  {
    std::map<std::string, int> childNames;
    for (pugi::xml_node child : node.children())
    {
      if (childNames.count(child.name()) == 0)
        childNames.insert({child.name(), 1});
      else
        childNames[child.name()]++;
    }
    return childNames;
  }

  static void InsertAsJsonAndRemove(nlohmann::json &json,
                                    pugi::xml_node &node,
                                    pugi::xml_node &child,
                                    bool insertInArray = false)
  {
    node.remove_child(child);
    FormatAndInsertJson(node.name(), child.value(), json, insertInArray);
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

  template <typename T>
  static void InsertJsonValue(const char *key,
                              T value,
                              nlohmann::json &json,
                              bool jsonIsArray)
  {
    if (jsonIsArray)
      json.push_back(value);
    else
      json[key] = value;
  }

  static void FormatAndInsertJson(const char *key,
                                  std::string value,
                                  nlohmann::json &json,
                                  bool jsonIsArray = false)
  {
    StripWhitespace(value);
    bool isNumeric, isNegative, hasFraction;
    MakeNumericAnalysis(value, isNumeric, isNegative, hasFraction);
    if (isNumeric)
    {
      if (hasFraction)
        InsertJsonValue(key, std::stod(value), json, jsonIsArray);
      else if (isNegative)
        InsertJsonValue(key, std::stoi(value), json, jsonIsArray);
      else
        InsertJsonValue(key, std::stoul(value), json, jsonIsArray);
    }
    else if (std::string(value) == "true")
      InsertJsonValue(key, true, json, jsonIsArray);
    else if (std::string(value) == "false")
      InsertJsonValue(key, false, json, jsonIsArray);
    else
      InsertJsonValue(key, value, json, jsonIsArray);
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
