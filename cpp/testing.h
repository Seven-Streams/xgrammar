/*!
 *  Copyright (c) 2024 by Contributors
 * \file xgrammar/testing.h
 * \brief The header testing utilities.
 */
#ifndef XGRAMMAR_TESTING_H_
#define XGRAMMAR_TESTING_H_

#include <xgrammar/xgrammar.h>

#include <string>
#include <vector>

namespace xgrammar {

std::string PrintTokenByIds(
    const std::vector<int32_t>& token_ids,
    const std::vector<std::pair<int32_t, std::string>>& sorted_decoded_vocab,
    int max_print_num
);

Grammar _EBNFToGrammarNoNormalization(
    const std::string& ebnf_string, const std::string& root_rule_name
);

std::string _PrintGrammarFSMs(const Grammar& grammar);

}  // namespace xgrammar

#endif  // XGRAMMAR_TESTING_H_
