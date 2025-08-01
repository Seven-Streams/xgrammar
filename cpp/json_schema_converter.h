/*!
 *  Copyright (c) 2024 by Contributors
 * \file xgrammar/json_schema_converter.h
 * \brief Convert a JSON schema string to EBNF grammar string.
 */

#ifndef XGRAMMAR_JSON_SCHEMA_CONVERTER_H_
#define XGRAMMAR_JSON_SCHEMA_CONVERTER_H_

#include <picojson.h>

#include <optional>
#include <string>
#include <utility>

namespace xgrammar {

/*!
 * \brief Convert JSON schema string to EBNF grammar string.
 * \param schema The JSON schema string.
 * \param indent The number of spaces for indentation. If set to std::nullopt, the output will be
 * in one line. Default: 2.
 * \param separators Two separators used in the schema: comma and colon. Examples: {",", ":"},
 * {", ", ": "}. If std::nullopt, the default separators will be used: {",", ": "} when the
 * indent is not -1, and {", ", ": "} otherwise. This follows the convention in python
 * json.dumps(). Default: std::nullopt.
 * \param strict_mode Whether to use strict mode. In strict
 * mode, the generated grammar will not allow properties and items that is not specified in the
 * schema. This is equivalent to setting unevaluatedProperties and unevaluatedItems to false.
 *
 * This helps LLM to generate accurate output in the grammar-guided generation with JSON
 * schema. Default: true.
 * \returns The EBNF grammar string.
 */
std::string JSONSchemaToEBNF(
    const std::string& schema,
    bool any_whitespace = true,
    std::optional<int> indent = std::nullopt,
    std::optional<std::pair<std::string, std::string>> separators = std::nullopt,
    bool strict_mode = true
);

/*!
 * \brief Convert JSON schema string to EBNF grammar string.
 * \param schema The JSON schema object.
 * \param indent The number of spaces for indentation. If set to std::nullopt, the output will be
 * in one line. Default: 2.
 * \param separators Two separators used in the schema: comma and colon. Examples: {",", ":"},
 * {", ", ": "}. If std::nullopt, the default separators will be used: {",", ": "} when the
 * indent is not -1, and {", ", ": "} otherwise. This follows the convention in python
 * json.dumps(). Default: std::nullopt. \param strict_mode Whether to use strict mode. In strict
 * mode, the generated grammar will not allow properties and items that is not specified in the
 * schema. This is equivalent to setting unevaluatedProperties and unevaluatedItems to false.
 *
 * This helps LLM to generate accurate output in the grammar-guided generation with JSON
 * schema. Default: true.
 * \returns The EBNF grammar string.
 */
std::string JSONSchemaToEBNF(
    const picojson::value& schema,
    bool any_whitespace = true,
    std::optional<int> indent = std::nullopt,
    std::optional<std::pair<std::string, std::string>> separators = std::nullopt,
    bool strict_mode = true
);

/*!
 * \brief Generate regex pattern for integer/float range.
 * \param start The start of the range (inclusive). If null assume negative infinity.
 * \param end The end of the range (inclusive). If null assume infinity.
 * \returns The regex pattern that matches integers/floats in the given range.
 */
std::string GenerateRangeRegex(std::optional<int64_t> start, std::optional<int64_t> end);

std::string GenerateFloatRangeRegex(std::optional<double> start, std::optional<double> end);

}  // namespace xgrammar

#endif  // XGRAMMAR_JSON_SCHEMA_CONVERTER_H_
