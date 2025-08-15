/*!
 *  Copyright (c) 2024 by Contributors
 * \file xgrammar/structural_tag.cc
 */
#include <xgrammar/xgrammar.h>

#include <sstream>
#include <string>
#include <vector>

#include "grammar_impl.h"
#include "grammar_parser.h"
#include "support/encoding.h"

namespace xgrammar {

std::string PrintTokenByIds(
    const std::vector<int32_t>& token_ids,
    const std::vector<std::pair<int32_t, std::string>>& sorted_decoded_vocab,
    int max_print_num
) {
  std::stringstream ss;
  ss << "[";
  int print_num = std::min(static_cast<int>(token_ids.size()), max_print_num);
  for (int i = 0; i < print_num; ++i) {
    ss << "#" << token_ids[i] << " <";
    for (const auto& [id, str] : sorted_decoded_vocab) {
      if (token_ids[i] == id) {
        ss << EscapeString(str) << ">";
      }
    }
    if (i < print_num - 1) {
      ss << ", ";
    }
  }
  if (static_cast<int>(token_ids.size()) > max_print_num) {
    ss << ", ...";
  }
  ss << "]";
  return ss.str();
}

Grammar _EBNFToGrammarNoNormalization(
    const std::string& ebnf_string, const std::string& root_rule_name
) {
  return ParseEBNF(ebnf_string, root_rule_name);
}

std::string _PrintGrammarFSMs(const Grammar& grammar) {
  std::string result;
  for (int i = 0; i < grammar->NumRules(); i++) {
    result += "Rule " + std::to_string(i) + ": " + grammar->GetRule(i).name + ", FSM: ";
    if (grammar->per_rule_fsms[i].has_value()) {
      result += grammar->per_rule_fsms[i]->ToString();
    } else {
      result += "None";
    }
    result += "\n";
  }
  return result;
}

}  // namespace xgrammar
