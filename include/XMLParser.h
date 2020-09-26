//
// Created by Anton J Olsson on 2020-09-25.
//

#ifndef CORE_XMLPARSER_H
#define CORE_XMLPARSER_H

#include "Logging.h"
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>

namespace DTCC
{

class XMLParser
{
public:
  static nlohmann::json GetJsonFromXML(const std::string &filePath)
  {
    nlohmann::json json;
    std::ifstream xmlFile;
    xmlFile.open(filePath);
    if (!xmlFile.is_open())
    {
      Error("Couldn't open " + filePath);
    }

    Node root;
    CreateNode(root, xmlFile);

    // ParseElement(xmlFile, json);

    return json;
  }

private:
  static const size_t CHAR_STREAM_SIZE = 1000;

  struct Node
  {
    std::string name;
    std::string value;
    std::vector<std::string> children;
  };

  static std::string GetStartTag(std::istream &istream)
  {
    char charTag[CHAR_STREAM_SIZE];
    istream.getline(charTag, CHAR_STREAM_SIZE, '>');
    return std::string(charTag).substr(1);
  }

  static std::array<int, 2> GetNodeContentPositions(std::istream &istream,
                                                    const std::string &nodeName)
  {
    std::array<int, 2> startEndPos{};
    startEndPos[0] = istream.tellg();
    char endTag[CHAR_STREAM_SIZE];
    while (endTag != nodeName)
    {
      if (istream.eof())
        Error("Missing end-tag");

      istream.ignore(CHAR_STREAM_SIZE, '<');
      if (istream.get() != '/')
        continue;
      istream.getline(endTag, CHAR_STREAM_SIZE, '>');
    }

    startEndPos[1] =
        (int)istream.tellg() - nodeName.size() - std::string("</>").size();
    return startEndPos;
  }

  static void SetNodeChildren(Node &node,
                              std::istream &istream,
                              std::array<int, 2> startEndPos)
  {
  }

  static void CreateNode(Node &node, std::istream &istream)
  {
    std::string name = GetStartTag(istream);
    node.name = name;
    std::array<int, 2> startEndPos = GetNodeContentPositions(istream, name);
    istream.seekg(startEndPos[0]);
    if (istream.peek() == '<')
      SetNodeChildren(node, istream, startEndPos);
    else
    {
      char charValue[CHAR_STREAM_SIZE];
      istream.getline(charValue, CHAR_STREAM_SIZE, '<');
      node.value = std::string(charValue);
    }
  }

  /*static std::string GetString(std::ifstream &ifstream, char delimiter,
                               int startIndex)
  {
    char *element = new char;
    ifstream.getline(element, 10000, delimiter);
    std::string elemStr(element);
    return elemStr.substr(startIndex, elemStr.size());
  }

  static std::string ToLowerCase(const std::string& str)
  {
    std::string lowStr;
    for (char i : str)
    {
      lowStr += std::tolower(i);
    }
    return lowStr;
  }*/

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

  /*static void ParseElement(std::ifstream& ifstream, nlohmann::json& json)
  {
    std::string key = GetString(ifstream, '>', 1);
    std::string value = GetString(ifstream, '<', 0);

    bool isNumeric, isNegative, hasFraction;
    MakeNumericAnalysis(value, isNumeric, isNegative, hasFraction);
    if (isNumeric)
    {
      if (hasFraction) json[key] = std::stod(value);
      else if (isNegative) json[key] = std::stoi(value);
      else json[key] = (size_t ) std::stoul(value);
    }
    else if (ToLowerCase(value) == "true")
      json[key] = true;
    else if (ToLowerCase(value) == "false")
      json[key] = false;
    else json[key] = value;
  }*/
};

} // namespace DTCC

#endif // CORE_XMLPARSER_H
