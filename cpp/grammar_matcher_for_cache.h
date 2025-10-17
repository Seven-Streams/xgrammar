/*!
 *  Copyright (c) 2025 by Contributors
 * \file xgrammar/grammar_matcher_for_cache.h
 * \brief The header for the grammar matcher for the cache.
 */

#ifndef XGRAMMAR_GRAMMAR_MATCHER_FOR_CACHE_H_
#define XGRAMMAR_GRAMMAR_MATCHER_FOR_CACHE_H_

#include <bitset>
#include <cstdint>
#include <unordered_set>

#include "compiled_grammar_impl.h"
#include "earley_parser.h"
#include "grammar_functor.h"
#include "xgrammar/tokenizer_info.h"

namespace xgrammar {
/*! \brief The concrete implementation of GrammarMatcherNode. */
class GrammarMatcherForTokenMaskCache : public EarleyParser {
 public:
  GrammarMatcherForTokenMaskCache(
      const Grammar& grammar,
      const ParserState& init_state,
      const std::unordered_map<int32_t, DynamicBitset>&
          tag_dispatch_rule_id_to_second_slicing_bitset,
      const TokenizerInfo& tokenizer_info,
      const bool& need_expand = false
  )
      : EarleyParser(grammar, init_state, need_expand),
        init_rule_id_(init_state.rule_id),
        initial_state_(init_state),
        tokenizer_info_(tokenizer_info),
        tag_dispatch_rule_id_to_second_slicing_bitset_(tag_dispatch_rule_id_to_second_slicing_bitset
        ) {}
  /*!
   * \brief Get the adaptive token mask for the given ParserState.
   * \param is_root_rule Whether to consider the parent rule. If false, there will be
   * no uncertain tokens. Useful for the root rule.
   */
  AdaptiveTokenMask GetAdaptiveTokenMask(bool is_root_rule);

  static void ClearCache() { crossing_cache_manager_.ClearCache(); }

 private:
  /*! \brief Check if a token can pass the lookahead assertion. */
  std::pair</*acceptable*/ bool, /*can reach end*/ bool> IsTokenPassLookaheadAssertion(
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
  bool GetTokenMaskWithFirstCharacterCheck(
      const std::bitset<256>& first_char_mask, bool is_root_rule, bool crossing_cache_is_available
  );

  int GetLengthOfString(
      int current_state, std::unordered_set<int32_t>& accepted_str_size, int accepted_character
  );

  std::bitset<256> GetCurrentAcceptedCharacters(int current_state);

  int GetForceLength(
      int current_state,
      const std::bitset<256>& character_ranges,
      std::unordered_set<int32_t>& visited_states
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
  std::vector<int32_t> tmp_accepted_by_lookahead_indices_;
  std::vector<bool> tmp_can_reach_end_stack_;
  std::vector<bool> tmp_can_reach_end_prefix_or_stack_;

  const TokenizerInfo& tokenizer_info_;

  /*! \brief A static crossing cache manager, used for better efficiency. */
  static CrossingCacheManager crossing_cache_manager_;

  /*!
   \brief This is a mapping from TagDispatch rule id to the bitset used for second slicing.
   \note If a rule is a TagDispatch rule, then there will be an AC automaton for its triggers.
    Which means that it can accept a lot of tokens. However, it will be slow to check a lot of
    tokens. The DynamicBitset here is used to do a second slicing: if a token's substr(1, n - 1)
    can be accepted by the start state of the AC automaton, then it will be True in the bitset.
    When we check a token, we first check if its first character can transit to the start state.
    If yes, then we check if it is in the bitset. If yes, then we accept it directly.
  */
  const std::unordered_map<int32_t, DynamicBitset>& tag_dispatch_rule_id_to_second_slicing_bitset_;
};
}  // namespace xgrammar

#endif  // XGRAMMAR_GRAMMAR_MATCHER_FOR_CACHE_H_
