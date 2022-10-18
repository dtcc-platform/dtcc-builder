// Copyright (C) 2020 Anton J Olsson
// Licensed under the MIT License

#ifndef CORE_XMLPARSER_H
#define CORE_XMLPARSER_H

#include "Logging.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include <pugixml.hpp>
#include <string>

namespace DTCC_BUILDER
{

/// Class for reading an XML file and returning the data as a JSON object, with
/// inferred and formatted JSON data types.
class XMLParser
{

public:
  /// Return JSON object given an XML file. Root XML node can be excluded, as
  /// this creates yet another level of nesting and is often obvious from
  /// object/file name.
  ///
  /// \param filePath Path to the XML file
  /// \param includeRoot Whether the root node should be included in the JSON
  /// object or not
  /// \return File contents as a JSON object
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
  /// Get all children of XML node that have values (i.e. text).
  ///
  /// \param node Node to get children from
  /// \return A vector with children of type node_pcdata
  static std::vector<pugi::xml_node> GetTextChildren(pugi::xml_node node)
  {
    std::vector<pugi::xml_node> children;
    for (pugi::xml_node child : node.children())
      if (child.type() == pugi::node_pcdata)
        children.push_back(child);
    return children;
  }

  /// Recursively insert XML node's content into JSON object.
  ///
  /// \param xmlNode The node whose content will be parsed and inserted
  /// \param json The JSON object to insert the content into
  static void ParseNode(pugi::xml_node &xmlNode, nlohmann::json &json)
  {
    // Add all xmlNode's attributes to json
    for (pugi::xml_attribute attr : xmlNode.attributes())
      FormatAndInsertJson(attr.name(), attr.value(), json);
    std::map<std::string, int> childNames = GetChildrensNameCount(xmlNode);
    for (pugi::xml_node xmlChild : xmlNode.children())
    {
      nlohmann::json jsonChild = {};
      // Handle multiple child nodes with same name
      if (childNames[xmlChild.name()] > 1)
        PutSameNameNodeInArray(jsonChild, json, xmlChild);
      // Handle element with attributes and text child
      else if (xmlChild.type() == pugi::node_pcdata)
        FormatAndInsertJson("#content", xmlChild.value(), json);
      else
      {
        std::vector<pugi::xml_node> grandTextChildren =
            GetTextChildren(xmlChild);
        // Base case of node with unique name and text node as only child
        if (grandTextChildren.size() == 1 && GetNumChildren(xmlChild) == 1 &&
            xmlChild.first_attribute() == nullptr)
        {
          InsertAsJsonAndRemove(json, xmlChild, grandTextChildren[0], false);
        }
        else
        {
          // Handle multiple text nodes in node
          if (grandTextChildren.size() > 1)
            PutTextNodesInArray(jsonChild, grandTextChildren, xmlChild);
          json[xmlChild.name()] = jsonChild;
          ParseNode(xmlChild, json[xmlChild.name()]);
        }
      }
    }
  }

  /// If node has just one (text) child and no attributes, put text into array
  /// (create array first if it doesn't exist). Otherwise, parse node and put
  /// jsonChild with its content in the array.
  ///
  /// \param jsonChild JSON object on same level as xmlNode
  /// \param json jsonChild's parent
  /// \param xmlNode Child with same name as at least one sibling
  static void PutSameNameNodeInArray(nlohmann::json &jsonChild,
                                     nlohmann::json &json,
                                     pugi::xml_node &xmlNode)
  {
    bool hasJustOneValue =
        GetNumChildren(xmlNode) == 1 && xmlNode.first_attribute() == nullptr;
    if (!hasJustOneValue)
      ParseNode(xmlNode, jsonChild);
    if (json[xmlNode.name()].is_array())
      json[xmlNode.name()].push_back(
          hasJustOneValue ? xmlNode.first_child().value() : jsonChild);
    else
      json[xmlNode.name()] = {hasJustOneValue ? xmlNode.first_child().value()
                                              : jsonChild};
  }

  /// Put all children (text nodes) of XML node into JSON array, remove the
  /// children from node and insert array into JSON object.
  ///
  /// \param json The JSON object to insert the array into
  /// \param nodeChildren The nodes to insert in the array
  /// \param node The node whose children are to be inserted in the array
  static void
  PutTextNodesInArray(nlohmann::json &json,
                      const std::vector<pugi::xml_node> &nodeChildren,
                      pugi::xml_node &node)
  {
    nlohmann::json array;
    for (pugi::xml_node child : nodeChildren)
    {
      std::string valueStr(child.value());
      InsertAsJsonAndRemove(array, node, child, true);
    }
    json["#content"] = array;
  }

  /// Get number of children of XML node.
  ///
  /// \param node The XML node
  /// \return The number of children
  static int GetNumChildren(pugi::xml_node &node)
  {
    size_t n = 0;
    for (auto it = node.begin(); it != node.end(); it++)
      n++;
    return n;
  }

  /// Get the name frequency of a node's children.
  ///
  /// \param node The node whose children's names will be counted
  /// \return A map containing all child names and their frequency
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

  /// Insert an XML child node into a JSON object/array and then remove the
  /// child node from its parent.
  ///
  /// \param json The JSON to insert the child node into
  /// \param node The parent node
  /// \param child The node to insert and remove
  /// \param insertInArray Whether the node should be inserted into a JSON array
  /// or object.
  static void InsertAsJsonAndRemove(nlohmann::json &json,
                                    pugi::xml_node &node,
                                    pugi::xml_node &child,
                                    bool insertInArray = false)
  {
    node.remove_child(child);
    FormatAndInsertJson(node.name(), child.value(), json, insertInArray);
  }

  /// Determine whether a string is a number and if so, whether it's negative
  /// and/or has a fraction.
  ///
  /// \param value The string to analyze
  /// \param isNumeric Will be set true if the string is a number
  /// \param isNegative Will be set true if the string starts with '-'
  /// \param hasFraction Will be set true if the string has at least one decimal
  /// point
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

  /// Remove leading and trailing whitespace, newline and tab characters.
  /// Could be moved to Utils.h.
  ///
  /// \param string The string to trim
  static void TrimString(std::string &string)
  {
    string.erase(0, string.find_first_not_of(" \n\r\t"));
    string.erase(string.find_last_not_of(" \n\r\t") + 1);
  }

  /// Insert a value into a JSON object or array.
  ///
  /// \tparam T The type of JSON value to insert
  /// \param key The value's key (null if jsonIsArray)
  /// \param value The value to insert
  /// \param json The JSON object or array to insert value into
  /// \param jsonIsArray Whether json is an array or object
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

  /// Trim value string, determine its JSON type, format and insert it
  /// into a JSON object/array.
  ///
  /// \param key The key of the value
  /// \param value The non-formatted JSON value
  /// \param json The JSON object/array
  /// \param jsonIsArray Whether json is an array or object
  static void FormatAndInsertJson(const char *key,
                                  std::string value,
                                  nlohmann::json &json,
                                  bool jsonIsArray = false)
  {
    TrimString(value);
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

  /// Print an error message if there was an error when parsing the XML.
  ///
  /// \param result Details of the parsing result
  /// \param doc The xml_document that tried to parse the file
  /// \param path The XML file path
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
    error(ss.str());
  }
};

} // namespace DTCC_BUILDER

#endif // CORE_XMLPARSER_H
