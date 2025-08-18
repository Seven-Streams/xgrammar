/*!
 *  Copyright (c) 2024 by Contributors
 * \file xgrammar/compiler.cc
 */

#include <xgrammar/compiler.h>

#include <algorithm>
#include <bitset>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

#include "compiled_grammar_impl.h"
#include "earley_parser.h"
#include "fsm.h"
#include "grammar_functor.h"
#include "grammar_impl.h"
#include "support/dynamic_bitset.h"
#include "support/int_set.h"
#include "support/logging.h"
#include "support/thread_pool.h"
#include "support/thread_safe_cache.h"
#include "support/utils.h"
#include "xgrammar/grammar.h"
#include "xgrammar/tokenizer_info.h"

namespace std {

/*! \brief Define the hash function for StructuralTagItem. */
template <>
struct hash<xgrammar::StructuralTagItem> {
  size_t operator()(const xgrammar::StructuralTagItem& tag) const {
    return xgrammar::HashCombine(
        std::hash<std::string>{}(tag.begin),
        std::hash<std::string>{}(tag.schema),
        std::hash<std::string>{}(tag.end)
    );
  }
};

}  // namespace std

namespace xgrammar {

/************** Use GrammarMatcher to generate the AdaptiveTokenMaskCache **************/

/*! \brief The concrete implementation of GrammarMatcherNode. */
class GrammarMatcherForTokenMaskCache : public EarleyParser {
 public:
  GrammarMatcherForTokenMaskCache(
      const Grammar& grammar,
      const ParserState& init_state,
      const TokenizerInfo& tokenizer_info,
      const bool& need_expand = true
  )
      : EarleyParser(grammar, init_state),
        init_rule_id_(init_state.rule_id),
        initial_state_(init_state),
        tokenizer_info_(tokenizer_info) {}
  /*!
   * \brief Get the adaptive token mask for the given ParserState.
   * \param is_root_rule Whether to consider the parent rule. If false, there will be
   * no uncertain tokens. Useful for the root rule.
   */
  AdaptiveTokenMask GetAdaptiveTokenMask(bool is_root_rule);

  static void ClearCache() { crossing_cache_manager_.ClearCache(); }

 private:
  /*! \brief Check if a token can pass the lookahead assertion. */
  bool IsTokenPassLookaheadAssertion(
      const std::string& token, const std::vector<bool>& can_reach_end_stack
  );

  /*! \brief Update the cache with lookahead information. */
  void AdaptCacheWithLookahead(AdaptiveTokenMask& cache, bool is_root_rule);

  /*!
   * \brief Check if speculative calculation will be applied.
   * \return first: whether speculative calculation is applicable.
   * \return second: part of the first character mask,
   * which can be used in speculative calculation.
   */
  std::pair<bool, std::bitset<256>> GetSpeculativeCalculation();

  /*!
   * \brief Get the first character mask.
   */
  void GetFirstCharacterMask(std::bitset<256>& first_character_mask);

  /*!
   * \brief Get the token mask for the given ParserState.
   * \param sorted_decoded_vocab The sorted decoded vocabulary.
   * \param first_char_mask The first character mask.
   * \param is_root_rule Whether to consider the parent rule. If false, there will be
   * no uncertain tokens. Useful for the root rule.
   */
  void GetTokenMaskWithFirstCharacterCheck(
      const std::bitset<256>& first_char_mask,
      bool is_root_rule,
      bool crossing_cache_is_available,
      const std::optional<std::optional<AdaptiveTokenMask>>& crossing_cache
  );

  // The id of the initial rule.
  int32_t init_rule_id_;

  // The initial state of the parser.
  ParserState initial_state_;

  // Temporary data for GetAdaptiveTokenMask.
  std::vector<int32_t> tmp_accepted_indices_;
  std::vector<int32_t> tmp_rejected_indices_;
  std::vector<int32_t> tmp_uncertain_indices_;
  std::vector<int32_t> tmp_rejected_by_lookahead_indices_;
  std::vector<bool> tmp_can_reach_end_stack_;
  std::vector<bool> tmp_can_reach_end_prefix_or_stack_;

  const TokenizerInfo& tokenizer_info_;

  /*! \brief A static crossing cache manager, used for better efficiency. */
  static CrossingCacheManager crossing_cache_manager_;
};

CrossingCacheManager GrammarMatcherForTokenMaskCache::crossing_cache_manager_;

bool GrammarMatcherForTokenMaskCache::IsTokenPassLookaheadAssertion(
    const std::string& token, const std::vector<bool>& can_reach_end_stack
) {
  auto lookahead_assertion_id = grammar_->GetRule(init_rule_id_).lookahead_assertion_id;
  if (lookahead_assertion_id == -1) {
    return true;
  }
  auto lookahead_state =
      ParserState(/*rule_id*/ -1, lookahead_assertion_id, 0, ParserState::kNoPrevInputPos, 0);
  PushStateAndExpand(lookahead_state);
  int token_len = token.size();

  // Find all positions that can come to and end. Then check if the suffix from that position
  // can be accepted by the lookahead assertion.
  for (int i = static_cast<int>(can_reach_end_stack.size()) - 1; i >= 0; --i) {
    if (!can_reach_end_stack[i]) {
      continue;
    }
    int last_accept_pos = i - 1;
    for (int pos = i; pos < token_len; ++pos) {
      if (!Advance(token[pos])) {
        break;
      }
      last_accept_pos = pos;
      // Case 1. The whole rule is finished.
      if (IsCompleted()) {
        // accepted chars: pos - i + 1
        // we need to rollback the pushed initial state as well
        PopLastStates(pos - i + 2);
        return true;
      }
    }
    // Case 2. The whole token is accepted
    if (last_accept_pos == token_len - 1) {
      PopLastStates(last_accept_pos - i + 2);
      return true;
    }
    // Case 3. The token is not accepted. Check the next position.
    PopLastStates(last_accept_pos - i + 1);
  }

  PopLastStates(1);
  return false;
}

// Comparator for std::pair<int32_t, std::string> based on the string value.
class IntStringPairComparator {
 public:
  bool operator()(
      const std::pair<int32_t, std::string>& lhs, const std::pair<int32_t, std::string>& rhs
  ) const {
    return lhs.second < rhs.second;
  }
};

int GetPossibleTokenIntervals(
    const std::vector<std::pair<int32_t, std::string>>& sorted_decoded_vocab,
    const std::bitset<256>& first_char_mask,
    std::vector<std::pair<int32_t, int32_t>>& possible_intervals
) {
  int possible_token_num = 0;
  int matched_size = 0;
  int last_interval_end = -1;
  for (int32_t i = 0; i < 256; i++) {
    if (first_char_mask[i]) {
      if (last_interval_end == -1) {
        last_interval_end = i;
      }
    } else {
      if (last_interval_end != -1) {
        int32_t interval_left_end =
            std::lower_bound(
                sorted_decoded_vocab.begin() + matched_size,
                sorted_decoded_vocab.end(),
                std::make_pair(0, std::string(1, static_cast<uint8_t>(last_interval_end))),
                IntStringPairComparator()
            ) -
            sorted_decoded_vocab.begin();
        int32_t interval_right_end = std::lower_bound(
                                         sorted_decoded_vocab.begin() + interval_left_end,
                                         sorted_decoded_vocab.end(),
                                         std::make_pair(0, std::string(1, static_cast<uint8_t>(i))),
                                         IntStringPairComparator()
                                     ) -
                                     sorted_decoded_vocab.begin();
        possible_intervals.emplace_back(interval_left_end, interval_right_end);
        possible_token_num += interval_right_end - interval_left_end;
        last_interval_end = -1;
        matched_size = interval_right_end;
      }
    }
  }

  if (last_interval_end != -1) {
    // If the last interval is not closed, we need to close it.
    int32_t interval_left_end =
        std::lower_bound(
            sorted_decoded_vocab.begin() + matched_size,
            sorted_decoded_vocab.end(),
            std::make_pair(0, std::string(1, static_cast<uint8_t>(last_interval_end))),
            IntStringPairComparator()
        ) -
        sorted_decoded_vocab.begin();
    possible_intervals.emplace_back(interval_left_end, sorted_decoded_vocab.size());
    possible_token_num += sorted_decoded_vocab.size() - interval_left_end;
  }
  return possible_token_num;
}

std::pair<bool, std::bitset<256>> GrammarMatcherForTokenMaskCache::GetSpeculativeCalculation() {
  using GrammarExprType = Grammar::Impl::GrammarExprType;
  // Check if the initial state is self-recursive-like. If the state is self-recursive-like,
  // and it covers a large part of the vocabulary, we will do speculative calculation in compiling.
  if (!grammar_->per_rule_fsms[init_rule_id_].has_value()) {
    if (initial_state_.sub_element_id == 0) {
      const auto& sequence_expr = grammar_->GetGrammarExpr(initial_state_.sequence_id);
      // A self-recursive-like rule must be a sequence.
      if (sequence_expr.type == GrammarExprType::kSequence) {
        const auto& current_element_expr =
            grammar_->GetGrammarExpr(sequence_expr[initial_state_.element_id]);
        // If the current element is a character class star, then it's self-recursive without doubt.
        if (current_element_expr.type == GrammarExprType::kCharacterClassStar) {
          return {true, {}};
          // If the current element is a character class, and the next element is a rule ref to
          // itself, and the rule only has 2 elements, then it's self-recursive-like.
        } else if (current_element_expr.type == GrammarExprType::kCharacterClass &&
                   sequence_expr.size() == 2 && initial_state_.element_id == 0) {
          const auto& end_element_expr = grammar_->GetGrammarExpr(sequence_expr[1]);
          if (end_element_expr.type == GrammarExprType::kRuleRef &&
              end_element_expr[0] == initial_state_.rule_id) {
            return {true, {}};
          }
        }
      }
    }
    return {false, {}};
  }
  // If the initial state is a FSM, we will check if the FSM is self-recursive-like.
  bool can_be_applied = false;
  std::bitset<256> speculative_mask;
  const auto& fsm = grammar_->per_rule_fsms[init_rule_id_].value();
  XGRAMMAR_DCHECK(initial_state_.element_id < fsm.NumStates());
  for (const auto& edge : fsm.GetFsm().GetEdges(initial_state_.element_id)) {
    if (edge.IsCharRange()) {
      // Case A: The edge is towards itself.
      if (edge.target == initial_state_.element_id) {
        can_be_applied = true;
        for (int ch = edge.min; ch <= edge.max; ++ch) {
          speculative_mask.set(ch);
        }
        continue;
      }

      // Case B: The state is the start state, and there's an edge to another state,
      // which calls the fsm itself.
      if (fsm.GetStart() == initial_state_.element_id) {
        for (const auto& next_edge : fsm.GetFsm().GetEdges(edge.target)) {
          if (next_edge.IsRuleRef() && next_edge.GetRefRuleId() == init_rule_id_) {
            can_be_applied = true;
            for (int ch = edge.min; ch <= edge.max; ++ch) {
              speculative_mask.set(ch);
            }
            break;
          }
        }
      }
    }
  }
  return {can_be_applied, speculative_mask};
}

void GrammarMatcherForTokenMaskCache::GetTokenMaskWithFirstCharacterCheck(
    const std::bitset<256>& first_char_mask,
    bool is_root_rule,
    bool crossing_cache_is_available,
    const std::optional<std::optional<AdaptiveTokenMask>>& crossing_cache
) {
  const auto& sorted_decoded_vocab = tokenizer_info_.GetSortedDecodedVocab();
  const auto& subtree_nodes_range = tokenizer_info_.GetTrieSubtreeNodesRange();
  // the pair (a, b) means [a, b). Intialize the possible intervals.
  std::vector<std::pair<int32_t, int32_t>> possible_intervals;
  int possible_token_num =
      GetPossibleTokenIntervals(sorted_decoded_vocab, first_char_mask, possible_intervals);

  XGRAMMAR_DCHECK(possible_intervals.size() > 0)
      << "There should be at least one possible interval for the first character mask.";

  if (possible_intervals[0].first != 0) {
    for (int i = 0; i < possible_intervals[0].first; ++i) {
      tmp_rejected_indices_.push_back(i);
    }
  }

  bool speculative_calculation = false;
  std::bitset<256> speculative_mask;
  if (init_rule_id_ == -1 || !grammar_->per_rule_fsms[init_rule_id_].has_value()) {
    speculative_calculation =
        GetSpeculativeCalculation().first &&
        (possible_token_num >= static_cast<int>(sorted_decoded_vocab.size() / 4));
    speculative_mask = first_char_mask;
  } else {
    std::tie(speculative_calculation, speculative_mask) = GetSpeculativeCalculation();
  }

  int prev_matched_size = 0;
  int last_rejected_range = 0;
  const std::string* prev_token = nullptr;
  for (size_t interval_idx = 0; interval_idx < possible_intervals.size(); ++interval_idx) {
    const auto& interval = possible_intervals[interval_idx];
    for (int i = interval.first; i < interval.second; ++i) {
      // Check if the current token is in the rejected range. i.e. check if the current token
      // is on the subtree of the rejected token.
      if (i < last_rejected_range) {
        tmp_rejected_indices_.push_back(i);
        continue;
      }

      const auto& token = sorted_decoded_vocab[i].second;
      // This optimization is useful for simple self-recursive rules, like string content.
      if (speculative_calculation) {
        bool all_accepted = true;
        for (char ch : token) {
          // If the first character is not the ascii character or can't be accepted by the
          // first character mask, we need to check them in the parser.
          if (isascii(ch) == 0 || !speculative_mask[static_cast<uint8_t>(ch)]) {
            all_accepted = false;
            break;
          }
        }
        if (all_accepted) {
          tmp_accepted_indices_.push_back(i);
          continue;
        }
      }
      // Many tokens may contain the same prefix, so we will avoid unnecessary matching
      // by finding the longest common prefix with the previous token.
      bool accepted = true;
      if (prev_token != nullptr) {
        int lcp_len =
            std::mismatch(token.begin(), token.end(), prev_token->begin(), prev_token->end())
                .first -
            token.begin();
        if (lcp_len > prev_matched_size) {
          // Case 1. The common prefix is rejected by the matcher in the last token. Reject
          // directly.
          accepted = false;
        } else if (lcp_len < prev_matched_size) {
          // Case 2. The common prefix is shorter than the previous matched size. Rollback
          // the non-common part.
          PopLastStates(prev_matched_size - lcp_len);
          tmp_can_reach_end_stack_.erase(
              tmp_can_reach_end_stack_.end() - (prev_matched_size - lcp_len),
              tmp_can_reach_end_stack_.end()
          );
          tmp_can_reach_end_prefix_or_stack_.erase(
              tmp_can_reach_end_prefix_or_stack_.end() - (prev_matched_size - lcp_len),
              tmp_can_reach_end_prefix_or_stack_.end()
          );
        }
        prev_matched_size = std::min(prev_matched_size, lcp_len);
      }

      prev_token = &token;

      if (accepted) {
        // Accept the rest chars one by one.
        for (int j = prev_matched_size; j < static_cast<int>(token.size()); ++j) {
          if (!Advance(token[j])) {
            accepted = false;
            break;
          }
          tmp_can_reach_end_stack_.push_back(IsCompleted());
          tmp_can_reach_end_prefix_or_stack_.push_back(
              tmp_can_reach_end_stack_.back() || tmp_can_reach_end_prefix_or_stack_.back()
          );
          prev_matched_size = j + 1;
        }
      }

      bool can_reach_end = tmp_can_reach_end_prefix_or_stack_.back();

      if (accepted) {
        tmp_accepted_indices_.push_back(i);
      } else if (can_reach_end && prev_matched_size > 0) {
        if ((!is_root_rule) && IsTokenPassLookaheadAssertion(token, tmp_can_reach_end_stack_)) {
          tmp_uncertain_indices_.push_back(i);
        } else {
          for (int j = i; j < subtree_nodes_range[i]; j++) {
            tmp_rejected_indices_.push_back(j);
            tmp_rejected_by_lookahead_indices_.push_back(j);
          }
          i = subtree_nodes_range[i] - 1;  // Skip the subtree nodes.
        }
      } else {
        tmp_rejected_indices_.push_back(i);
        last_rejected_range = subtree_nodes_range[i];
      }
    }
    if (interval_idx != possible_intervals.size() - 1) {
      const auto& next_interval = possible_intervals[interval_idx + 1];
      for (int i = interval.second; i < next_interval.first; ++i) {
        tmp_rejected_indices_.push_back(i);
      }
    }
  }

  // Rollback the last matched part.
  PopLastStates(prev_matched_size);

  if (possible_intervals.back().second != static_cast<int>(sorted_decoded_vocab.size())) {
    // If the last interval is not closed, we need to reject the rest tokens.
    for (int i = possible_intervals.back().second;
         i < static_cast<int>(sorted_decoded_vocab.size());
         ++i) {
      tmp_rejected_indices_.push_back(i);
    }
  }
}

void GrammarMatcherForTokenMaskCache::GetFirstCharacterMask(std::bitset<256>& first_character_mask
) {
  first_character_mask.reset();
  const auto& sequence = grammar_->GetGrammarExpr(initial_state_.sequence_id);
  if (!grammar_->per_rule_fsms[init_rule_id_].has_value()) {
    const auto& sub_sequence = grammar_->GetGrammarExpr(sequence[initial_state_.element_id]);
    switch (sub_sequence.type) {
      case Grammar::Impl::GrammarExprType::kByteString: {
        first_character_mask[sub_sequence[initial_state_.sub_element_id]] = true;
        break;
      }
      case xgrammar::Grammar::Impl::GrammarExprType::kCharacterClass:
      case xgrammar::Grammar::Impl::GrammarExprType::kCharacterClassStar: {
        if (initial_state_.sub_element_id == 0) {
          bool is_negative = sub_sequence[0];
          for (int i = 1; i < sub_sequence.size(); i += 2) {
            int left_char = static_cast<uint8_t>(sub_sequence[i]);
            int right_char = static_cast<uint8_t>(sub_sequence[i + 1]);
            for (int c = left_char; c <= right_char; ++c) {
              first_character_mask[c] = true;
            }
          }
          if (is_negative) {
            first_character_mask = ~first_character_mask;
          }
          break;
        }
        // Otherwise, it's matching a UTF-8 character. We can optimize the matching process
        // here.
        for (size_t i = 0x80; i < 0xC0; ++i) {
          first_character_mask[i] = true;
        }
        break;
      }
      default: {
        XGRAMMAR_LOG(FATAL) << "Unsupported grammar expr type: " << static_cast<int>(sequence.type);
      }
    }
  } else {
    const auto& fsm = grammar_->per_rule_fsms[init_rule_id_].value();
    const auto& edges = fsm.GetFsm().GetEdges(initial_state_.element_id);
    for (const auto& edge : edges) {
      if (edge.IsCharRange()) {
        for (int c = edge.min; c <= edge.max; ++c) {
          first_character_mask[c] = true;
        }
      }
    }
  }
}

AdaptiveTokenMask GrammarMatcherForTokenMaskCache::GetAdaptiveTokenMask(bool is_root_rule) {
  tmp_accepted_indices_.clear();
  tmp_rejected_indices_.clear();
  tmp_uncertain_indices_.clear();
  tmp_rejected_by_lookahead_indices_.clear();
  tmp_can_reach_end_prefix_or_stack_.clear();
  tmp_can_reach_end_stack_.clear();
  // For every character in the current token, stores whether it is possible to reach the end of
  // the rule when matching until this character. Store it in a stack for later rollback.
  tmp_can_reach_end_stack_.push_back(false);
  tmp_can_reach_end_prefix_or_stack_.push_back(false);

  // Try to get the crossing cache.
  bool crossing_cache_is_available = grammar_->per_rule_fsm_hashes[init_rule_id_].has_value();
  std::optional<AdaptiveTokenMask> crossing_cache = std::nullopt;
  const auto& original_to_new_id = grammar_->per_rule_fsm_new_state_ids[init_rule_id_];
  uint64_t fsm_hash = -1;
  uint32_t new_state_id = -1;
  if (crossing_cache_is_available) {
    fsm_hash = grammar_->per_rule_fsm_hashes[init_rule_id_].value();
    auto get_new_state_id = std::find_if(
        original_to_new_id->begin(),
        original_to_new_id->end(),
        [&](const auto& original_new_pair) {
          return original_new_pair.first == initial_state_.element_id;
        }
    );
    XGRAMMAR_DCHECK(get_new_state_id != original_to_new_id->end());
    new_state_id = get_new_state_id->second;
    crossing_cache = crossing_cache_manager_.GetCache(
        fsm_hash, new_state_id, tokenizer_info_.GetTokenizerHash()
    );
    // If the rule doesn't have a lookahead, then it is exactly the same fsm.
    if (crossing_cache.has_value()) {
      AdaptCacheWithLookahead(crossing_cache.value(), is_root_rule);
      return std::move(crossing_cache.value());
    }
  }
  std::bitset<256> first_character_mask;
  GetFirstCharacterMask(first_character_mask);

  GetTokenMaskWithFirstCharacterCheck(
      first_character_mask, is_root_rule, crossing_cache_is_available, crossing_cache
  );

  auto return_value = AdaptiveTokenMask(
      tokenizer_info_.GetVocabSize(),
      tokenizer_info_.GetSortedDecodedVocab(),
      tmp_accepted_indices_,
      tmp_rejected_indices_,
      tmp_uncertain_indices_
  );
  if (crossing_cache_is_available && crossing_cache == std::nullopt) {
    // We can add a cache.
    // All the tokens rejected by lookahead should be uncertain.
    IntsetUnion(&tmp_uncertain_indices_, tmp_rejected_by_lookahead_indices_);
    std::vector<int32_t> rejected_indices_without_lookahead;
    rejected_indices_without_lookahead.reserve(
        tmp_rejected_indices_.size() - tmp_rejected_by_lookahead_indices_.size()
    );
    std::set_difference(
        tmp_rejected_indices_.begin(),
        tmp_rejected_indices_.end(),
        tmp_rejected_by_lookahead_indices_.begin(),
        tmp_rejected_by_lookahead_indices_.end(),
        std::back_inserter(rejected_indices_without_lookahead)
    );
    crossing_cache_manager_.AddCache(
        fsm_hash,
        new_state_id,
        tokenizer_info_.GetTokenizerHash(),
        AdaptiveTokenMask(
            tokenizer_info_.GetVocabSize(),
            tokenizer_info_.GetSortedDecodedVocab(),
            tmp_accepted_indices_,
            rejected_indices_without_lookahead,
            tmp_uncertain_indices_
        )
    );
  }
  return return_value;
}

void GrammarMatcherForTokenMaskCache::AdaptCacheWithLookahead(
    AdaptiveTokenMask& cache, bool is_root_rule
) {
  const auto& sorted_decoded_vocab = tokenizer_info_.GetSortedDecodedVocab();
  const auto& subtree_nodes_range = tokenizer_info_.GetTrieSubtreeNodesRange();
  const std::string* prev_token = nullptr;
  int prev_matched_size = 0;
  int last_rejected_range = 0;
  for (const auto& uncertain_index : cache.uncertain_indices) {
    const auto& token = sorted_decoded_vocab[uncertain_index].second;
    // Many tokens may contain the same prefix, so we will avoid unnecessary matching
    // by finding the longest common prefix with the previous token.
    bool accepted = true;
    if (uncertain_index < last_rejected_range) {
      tmp_rejected_indices_.push_back(uncertain_index);
      continue;
    }
    if (prev_token != nullptr) {
      int lcp_len =
          std::mismatch(token.begin(), token.end(), prev_token->begin(), prev_token->end()).first -
          token.begin();
      if (lcp_len > prev_matched_size) {
        // Case 1. The common prefix is rejected by the matcher in the last token. Reject
        // directly.
        accepted = false;
      } else if (lcp_len < prev_matched_size) {
        // Case 2. The common prefix is shorter than the previous matched size. Rollback
        // the non-common part.
        PopLastStates(prev_matched_size - lcp_len);
        tmp_can_reach_end_stack_.erase(
            tmp_can_reach_end_stack_.end() - (prev_matched_size - lcp_len),
            tmp_can_reach_end_stack_.end()
        );
        tmp_can_reach_end_prefix_or_stack_.erase(
            tmp_can_reach_end_prefix_or_stack_.end() - (prev_matched_size - lcp_len),
            tmp_can_reach_end_prefix_or_stack_.end()
        );
      }
      prev_matched_size = std::min(prev_matched_size, lcp_len);
    }

    prev_token = &token;

    if (accepted) {
      // Accept the rest chars one by one.
      for (int j = prev_matched_size; j < static_cast<int>(token.size()); ++j) {
        if (!Advance(token[j])) {
          accepted = false;
          break;
        }
        tmp_can_reach_end_stack_.push_back(IsCompleted());
        tmp_can_reach_end_prefix_or_stack_.push_back(
            tmp_can_reach_end_stack_.back() || tmp_can_reach_end_prefix_or_stack_.back()
        );
        prev_matched_size = j + 1;
      }
    }

    bool can_reach_end = tmp_can_reach_end_prefix_or_stack_.back();

    // All the tokens are at least uncertain!
    XGRAMMAR_CHECK(!accepted);

    if (can_reach_end && !is_root_rule &&
        IsTokenPassLookaheadAssertion(token, tmp_can_reach_end_stack_) && prev_matched_size > 0) {
      continue;
    } else {
      tmp_rejected_indices_.push_back(uncertain_index);
      last_rejected_range = subtree_nodes_range[uncertain_index];
    }
  }

  // Update the uncertain tokens.
  tmp_uncertain_indices_.reserve(cache.uncertain_indices.size() - tmp_rejected_indices_.size());
  std::set_difference(
      cache.uncertain_indices.begin(),
      cache.uncertain_indices.end(),
      tmp_rejected_indices_.begin(),
      tmp_rejected_indices_.end(),
      std::back_inserter(tmp_uncertain_indices_)
  );

  cache.uncertain_indices = tmp_uncertain_indices_;
  switch (cache.store_type) {
    case AdaptiveTokenMask::StoreType::kAccepted:
    case AdaptiveTokenMask::StoreType::kAcceptedBitset: {
      // Rejected tokens will be more. So it's still kAcceptedBitset or kAccept.
      break;
    }
    case AdaptiveTokenMask::StoreType::kRejected: {
      // Thus, rejected tokens will be more. We need to re-check its type.
      IntsetUnion(&cache.rejected_indices, tmp_rejected_indices_);
      int accepted_token_size = sorted_decoded_vocab.size() - cache.rejected_indices.size() -
                                cache.uncertain_indices.size();
      if (cache.rejected_indices.size() >= AdaptiveTokenMask::USE_BITSET_THRESHOLD &&
          accepted_token_size >= AdaptiveTokenMask::USE_BITSET_THRESHOLD) {
        // In this case, we need to switch to bitset representation.
        cache.store_type = AdaptiveTokenMask::StoreType::kAcceptedBitset;
        cache.accepted_bitset = DynamicBitset(tokenizer_info_.GetVocabSize());
        cache.accepted_bitset.Set();
        for (const auto& reject_index : cache.rejected_indices) {
          cache.accepted_bitset.Reset(sorted_decoded_vocab[reject_index].first);
        }
        for (const auto& uncertain_index : cache.uncertain_indices) {
          cache.accepted_bitset.Reset(sorted_decoded_vocab[uncertain_index].first);
        }
        for (const auto& stop_id : tokenizer_info_.GetStopTokenIds()) {
          cache.accepted_bitset.Reset(stop_id);
        }
        for (const auto& special_id : tokenizer_info_.GetSpecialTokenIds()) {
          cache.accepted_bitset.Reset(special_id);
        }
        cache.rejected_indices.clear();
      } else {
        if (accepted_token_size < static_cast<int>(cache.rejected_indices.size())) {
          // Switch to kAccepted.
          cache.store_type = AdaptiveTokenMask::StoreType::kAccepted;

          // Use a dynamicBitset to reset the state.
          DynamicBitset bit_set(sorted_decoded_vocab.size());
          for (const auto& reject_index : cache.rejected_indices) {
            bit_set.Set(reject_index);
          }
          for (const auto& uncertain_index : cache.uncertain_indices) {
            bit_set.Set(uncertain_index);
          }
          for (int i = 0; i < static_cast<int>(sorted_decoded_vocab.size()); ++i) {
            if (!bit_set[i]) {
              cache.accepted_indices.push_back(i);
            }
          }
          cache.rejected_indices.clear();
        } else {
          // Just keep kRejected.
        }
      }
      break;
    }
    default: {
      XGRAMMAR_LOG(FATAL) << "Unsupported store type: " << static_cast<int>(cache.store_type);
    }
  }
}

/******************* GrammarCompiler::Impl *******************/

using SchemaKey =
    std::tuple<std::string, bool, std::optional<int>, std::pair<std::string, std::string>, bool>;
using StructuralTagKey = std::tuple<std::vector<StructuralTagItem>, std::vector<std::string>>;
using GrammarKey = std::pair<std::string, std::string>;

class GrammarCompiler::Impl {
 public:
  Impl(
      const TokenizerInfo& tokenizer_info,
      int max_threads,
      bool cache_enabled,
      long long max_memory_bytes
  )
      : tokenizer_info_(tokenizer_info),
        max_threads_(max_threads),
        cache_enabled_(cache_enabled),
        compile_builtin_json_grammar_cache_([&] { return CompileJson(); }),
        compile_cache_(static_cast<std::size_t>(max_memory_bytes), *this) {}

  /*!
   * \brief Build the fsm for each rule and store in the grammar.
   */
  void BuildFSM(Grammar grammar);

  /*! \brief Multi-thread compile the grammar. */
  CompiledGrammar MultiThreadCompileGrammar(Grammar grammar);

  /*! \brief Compile the built-in JSON grammar. */
  CompiledGrammar CompileJson();

  /*!
   * \brief Compile different types of grammars.
   * \attention This template function is marked as deleted.
   * User must explicitly specialize the template to support new key types.
   */
  template <typename KeyType>
  CompiledGrammar Compute(const KeyType& key) = delete;

  /*! \brief Forwards the key to the corresponding compile function. */
  template <typename KeyType>
  CompiledGrammar operator()(const KeyType& key) {
    return Compute<KeyType>(key);
  }

  CompiledGrammar CompileBuiltinJSONGrammar();

  CompiledGrammar CompileJSONSchema(
      const std::string& schema,
      bool any_whitespace,
      std::optional<int> indent,
      std::optional<std::pair<std::string, std::string>> separators,
      bool strict_mode = true
  );

  CompiledGrammar CompileStructuralTag(
      const std::vector<StructuralTagItem>& tags, const std::vector<std::string>& triggers
  );

  CompiledGrammar CompileRegex(const std::string& regex);

  CompiledGrammar CompileGrammar(const Grammar& grammar);

  void ClearCache();

  long long GetCacheSizeBytes() const;
  long long CacheLimitBytes() const;

 private:
  using GrammarExprType = Grammar::Impl::GrammarExprType;
  using MultipleKey = std::variant<SchemaKey, StructuralTagKey, std::string, GrammarKey>;

  struct Computer {
    Computer(Impl& compiler) : compiler(compiler) {}
    // dispatch the key to the corresponding compile function
    CompiledGrammar operator()(const MultipleKey& key) const { return std::visit(compiler, key); }
    GrammarCompiler::Impl& compiler;
  };

  struct SizeEstimator {
    std::size_t operator()(const CompiledGrammar& value) const { return value.MemorySizeBytes(); }
  };

  /*! \brief The vocabulary associated with this storage class. */
  const TokenizerInfo tokenizer_info_;
  /*! \brief The maximum number of threads to use. */
  const int max_threads_;
  /*! \brief Whether the cache is enabled. */
  const bool cache_enabled_;

  ThreadSafeCache<CompiledGrammar> compile_builtin_json_grammar_cache_;
  ThreadSafeLRUCache<MultipleKey, CompiledGrammar, Computer, SizeEstimator> compile_cache_;
};

CompiledGrammar GrammarCompiler::Impl::MultiThreadCompileGrammar(Grammar grammar) {
  auto compiled_grammar_impl = std::make_shared<CompiledGrammar::Impl>();

  compiled_grammar_impl->grammar = grammar;
  compiled_grammar_impl->tokenizer_info = tokenizer_info_;
  grammar->allow_empty_rule_ids = AllowEmptyRuleAnalyzer::Apply(compiled_grammar_impl->grammar);
  RepetitionNormalizer::Apply(&compiled_grammar_impl->grammar);
  GrammarFSMBuilder::Apply(&compiled_grammar_impl->grammar);
  GrammarFSMHasher::Apply(&compiled_grammar_impl->grammar);
  if (tokenizer_info_.GetVocabSize() == 0) {
    return CompiledGrammar(compiled_grammar_impl);
  }
  // Step 3. Compute the adaptive token mask cache
  // The token mask cache is computed for these positions in the grammar:
  // 1. All character class or character class star (with last_utf8_bytes=0, 1, 2, 3)
  // 2. All byte strings (with element_in_string=0, 1, 2, ...)
  // since other positions will be expanded to the above positions

  // TODO(Charlie): Figure out how to support ThreadPool and std::mutex in WebAssembly.
  // Only declare ThreadPool and mutex if max_threads > 1, so when max_threads = 1, we do
  // not need ThreadPool or std::mutex, which throws error in runtime in WebAssembly.
  std::optional<ThreadPool> thread_pool;
  std::optional<std::mutex> adaptive_token_mask_cache_mutex;

  if (max_threads_ > 1) {
    thread_pool.emplace(max_threads_);
    adaptive_token_mask_cache_mutex.emplace();
  }

  auto add_adaptive_token_mask = [&](const ParserState& state, bool is_root_rule) {
    auto grammar_matcher = GrammarMatcherForTokenMaskCache(grammar, state, tokenizer_info_, false);
    auto cur_adaptive_token_mask_cache = grammar_matcher.GetAdaptiveTokenMask(is_root_rule);
    if (max_threads_ > 1) {
      std::lock_guard<std::mutex> lock(adaptive_token_mask_cache_mutex.value());
      compiled_grammar_impl->adaptive_token_mask_cache[state] = cur_adaptive_token_mask_cache;
    } else {
      compiled_grammar_impl->adaptive_token_mask_cache[state] = cur_adaptive_token_mask_cache;
    }
  };

  auto add_task_adaptive_token_mask = [&](const ParserState& state, bool is_root_rule) {
    // Execute depending on whether we use thread_pool
    if (max_threads_ > 1) {
      thread_pool->Execute([add_adaptive_token_mask, state, is_root_rule]() {
        add_adaptive_token_mask(state, is_root_rule);
      });
    } else {
      add_adaptive_token_mask(state, is_root_rule);
    }
  };

  auto root_rule_id = grammar->GetRootRuleId();

  for (int32_t rule_id = 0; rule_id < static_cast<int>(grammar->NumRules()); ++rule_id) {
    auto rule = grammar->GetRule(rule_id);
    auto rule_body = grammar->GetGrammarExpr(rule.body_expr_id);
    const auto& rule_fsm = grammar->per_rule_fsms[rule_id];
    if (rule_fsm.has_value()) {
      auto cur_stack_element =
          ParserState(rule_id, rule.body_expr_id, 0, ParserState::kNoPrevInputPos, 0);
      std::unordered_set<int> reachable_states;
      rule_fsm->GetReachableStates(&reachable_states);
      for (int i : reachable_states) {
        cur_stack_element.element_id = i;
        if (!rule_fsm->IsScanableState(i)) {
          continue;
        }
        add_task_adaptive_token_mask(cur_stack_element, rule_id == root_rule_id);
      }
      continue;
    }
    XGRAMMAR_DCHECK(rule_body.type == GrammarExprType::kChoices);
    for (auto sequence_id : rule_body) {
      const auto& sequence = grammar->GetGrammarExpr(sequence_id);
      if (sequence.type == GrammarExprType::kEmptyStr) {
        continue;
      }
      XGRAMMAR_DCHECK(sequence.type == GrammarExprType::kSequence);
      auto state = ParserState(rule_id, sequence_id, 0, ParserState::kNoPrevInputPos, 0);
      for (int element_id = 0; element_id < sequence.size(); ++element_id) {
        state.element_id = element_id;
        auto element = grammar->GetGrammarExpr(sequence[element_id]);
        if (element.type == GrammarExprType::kRuleRef || element.type == GrammarExprType::kRepeat) {
          continue;
        }
        if (element.type == GrammarExprType::kByteString) {
          for (int idx = 0; idx < element.size(); ++idx) {
            state.sub_element_id = idx;
            add_task_adaptive_token_mask(state, rule_id == root_rule_id);
          }
        } else {
          XGRAMMAR_DCHECK(
              element.type == GrammarExprType::kCharacterClassStar ||
              element.type == GrammarExprType::kCharacterClass
          );
          for (int left_utf8_bytes = 0; left_utf8_bytes <= 3; ++left_utf8_bytes) {
            state.sub_element_id = left_utf8_bytes;
            add_task_adaptive_token_mask(state, rule_id == root_rule_id);
          }
        }
      }
    }
  }

  if (max_threads_ > 1) {
    thread_pool->Join();
  }

  return CompiledGrammar(compiled_grammar_impl);
}

CompiledGrammar GrammarCompiler::Impl::CompileJson() {
  return MultiThreadCompileGrammar(Grammar::BuiltinJSONGrammar());
}

template <>
CompiledGrammar GrammarCompiler::Impl::Compute<SchemaKey>(const SchemaKey& key) {
  const auto& [schema, any_whitespace, indent, separators, strict_mode] = key;
  auto grammar = Grammar::FromJSONSchema(schema, any_whitespace, indent, separators, strict_mode);
  return MultiThreadCompileGrammar(grammar);
}

template <>
CompiledGrammar GrammarCompiler::Impl::Compute<StructuralTagKey>(const StructuralTagKey& key) {
  const auto& [tags, triggers] = key;
  return MultiThreadCompileGrammar(Grammar::FromStructuralTag(tags, triggers));
}

template <>
CompiledGrammar GrammarCompiler::Impl::Compute<std::string>(const std::string& key) {
  return MultiThreadCompileGrammar(Grammar::FromRegex(key));
}

template <>
CompiledGrammar GrammarCompiler::Impl::Compute<GrammarKey>(const GrammarKey& key) {
  const auto& [grammar_str, root_rule_name] = key;
  return MultiThreadCompileGrammar(Grammar::FromEBNF(grammar_str, root_rule_name));
}

CompiledGrammar GrammarCompiler::Impl::CompileBuiltinJSONGrammar() {
  if (!cache_enabled_) {
    return MultiThreadCompileGrammar(Grammar::BuiltinJSONGrammar());
  }
  return compile_builtin_json_grammar_cache_.Get();
}

CompiledGrammar GrammarCompiler::Impl::CompileJSONSchema(
    const std::string& schema,
    bool any_whitespace,
    std::optional<int> indent,
    std::optional<std::pair<std::string, std::string>> separators,
    bool strict_mode
) {
  if (!cache_enabled_) {
    return MultiThreadCompileGrammar(
        Grammar::FromJSONSchema(schema, any_whitespace, indent, separators, strict_mode)
    );
  }
  auto separators_value = separators.value_or(
      (indent == std::nullopt) ? std::make_pair(", ", ": ") : std::make_pair(",", ": ")
  );
  auto key = std::make_tuple(schema, any_whitespace, indent, separators_value, strict_mode);
  return compile_cache_.Get(key);
}

CompiledGrammar GrammarCompiler::Impl::CompileStructuralTag(
    const std::vector<StructuralTagItem>& tags, const std::vector<std::string>& triggers
) {
  if (!cache_enabled_) {
    return MultiThreadCompileGrammar(Grammar::FromStructuralTag(tags, triggers));
  }
  auto key = std::make_tuple(tags, triggers);
  return compile_cache_.Get(key);
}

CompiledGrammar GrammarCompiler::Impl::CompileRegex(const std::string& regex) {
  if (!cache_enabled_) {
    return MultiThreadCompileGrammar(Grammar::FromRegex(regex));
  }
  return compile_cache_.Get(regex);
}

CompiledGrammar GrammarCompiler::Impl::CompileGrammar(const Grammar& grammar) {
  if (!cache_enabled_) {
    return MultiThreadCompileGrammar(grammar);
  }
  auto key = std::make_pair(grammar.ToString(), grammar->GetRootRule().name);
  return compile_cache_.Get(key);
}

void GrammarCompiler::Impl::ClearCache() {
  compile_builtin_json_grammar_cache_.Clear();
  compile_cache_.Clear();
  GrammarMatcherForTokenMaskCache::ClearCache();
}

long long GrammarCompiler::Impl::GetCacheSizeBytes() const {
  return static_cast<long long>(compile_cache_.MemorySize());
}

long long GrammarCompiler::Impl::CacheLimitBytes() const {
  const auto size = compile_cache_.MaxMemorySize();
  if (size == compile_cache_.UNLIMITED_SIZE) return -1;
  return static_cast<long long>(size);
}

/******************* GrammarCompiler *******************/

GrammarCompiler::GrammarCompiler(
    const TokenizerInfo& tokenizer_info,
    int max_threads,
    bool cache_enabled,
    long long max_memory_bytes
)
    : pimpl_(std::make_shared<Impl>(tokenizer_info, max_threads, cache_enabled, max_memory_bytes)) {
  if (max_memory_bytes < -1) {
    XGRAMMAR_LOG(FATAL) << "Invalid max_memory_bytes: " << max_memory_bytes << ". "
                        << "It should be -1 (unlimited) or a non-negative integer.";
  }
}

CompiledGrammar GrammarCompiler::CompileJSONSchema(
    const std::string& schema,
    bool any_whitespace,
    std::optional<int> indent,
    std::optional<std::pair<std::string, std::string>> separators,
    bool strict_mode
) {
  return pimpl_->CompileJSONSchema(schema, any_whitespace, indent, separators, strict_mode);
}

CompiledGrammar GrammarCompiler::CompileBuiltinJSONGrammar() {
  return pimpl_->CompileBuiltinJSONGrammar();
}

CompiledGrammar GrammarCompiler::CompileStructuralTag(
    const std::vector<StructuralTagItem>& tags, const std::vector<std::string>& triggers
) {
  return pimpl_->CompileStructuralTag(tags, triggers);
}

CompiledGrammar GrammarCompiler::CompileRegex(const std::string& regex) {
  return pimpl_->CompileRegex(regex);
}

CompiledGrammar GrammarCompiler::CompileGrammar(const Grammar& grammar) {
  return pimpl_->CompileGrammar(grammar);
}

void GrammarCompiler::ClearCache() { pimpl_->ClearCache(); }

long long GrammarCompiler::GetCacheSizeBytes() const { return pimpl_->GetCacheSizeBytes(); }

long long GrammarCompiler::CacheLimitBytes() const { return pimpl_->CacheLimitBytes(); }

}  // namespace xgrammar
