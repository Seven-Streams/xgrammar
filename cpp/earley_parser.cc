/*!
 *  Copyright (c) 2025 by Contributors
 * \file xgrammar/earley_parser.cc
 */

#include "earley_parser.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <ctime>
#include <utility>
#include <vector>

#include "fsm.h"
#include "grammar_impl.h"
#include "support/encoding.h"
#include "support/logging.h"
#include "xgrammar/grammar.h"

namespace xgrammar {

using GrammarExprType = Grammar::Impl::GrammarExprType;

using GrammarExpr = Grammar::Impl::GrammarExpr;

bool EarleyParser::IsCompleted() const { return is_completed_.back(); }

void EarleyParser::PopLastStates(int32_t cnt) {
  if (stop_token_is_accepted_) {
    stop_token_is_accepted_ = false;
  }
  if (cnt >= static_cast<int32_t>(rule_id_to_completeable_states_.size())) {
    XGRAMMAR_LOG(FATAL) << "The number of states to be popped is larger than the size of states.";
  }
  rule_id_to_completeable_states_.erase(
      rule_id_to_completeable_states_.end() - cnt, rule_id_to_completeable_states_.end()
  );
  is_completed_.erase(is_completed_.end() - cnt, is_completed_.end());
  scanable_state_history_.PopBack(cnt);
}

void EarleyParser::Complete(const ParserState& state, const GrammarExpr& grammar_expr) {
  // Check if a rule is completed.
  if (state.rule_start_pos == ParserState::kNoPrevInputPos) {
    // assert: if a root rule can achieve here, then it must be completed.
    tmp_accept_stop_token_ = true;
    return;
  }
  // Check all the possible parent states.
  const auto& parent_states_map = rule_id_to_completeable_states_[state.rule_start_pos];
  for (auto parent_state_iter = parent_states_map.lower_bound(state.rule_id);
       parent_state_iter != parent_states_map.end() && parent_state_iter->first == state.rule_id;
       parent_state_iter++) {
    const auto& parent_state = parent_state_iter->second;
    const auto& parent_expr = grammar_->GetGrammarExpr(parent_state.sequence_id);
    if (parent_state.rule_id == -1 || !grammar_->per_rule_fsms[parent_state.rule_id].has_value()) {
      const auto& element_expr = grammar_->GetGrammarExpr(parent_expr[parent_state.element_id]);
      // The new rule is not referenced by a fsm.
      XGRAMMAR_DCHECK(
          element_expr.type == GrammarExprType::kRuleRef ||
          element_expr.type == GrammarExprType::kRepeat
      );
      if (element_expr.type == GrammarExprType::kRuleRef) {
        Enqueue(ParserState{
            parent_state.rule_id,
            parent_state.sequence_id,
            parent_state.element_id + 1,
            parent_state.rule_start_pos,
            0
        });
        continue;
      }
      XGRAMMAR_DCHECK(element_expr.type == GrammarExprType::kRepeat);
      if (state.rule_start_pos ==
              static_cast<int32_t>(rule_id_to_completeable_states_.size() - 1) &&
          std::binary_search(
              grammar_->allow_empty_rule_ids.begin(),
              grammar_->allow_empty_rule_ids.end(),
              element_expr[0]
          )) {
        // It means that the subrule of the repeat is empty, and we have already detected it.
        // We shouldn't add it into the queue.
        continue;
      }
      // The parent state is a repeat, we need to increase the repeat count.
      auto new_state = parent_state;
      const int32_t& min_repeat_count = element_expr[1];
      const int32_t& max_repeat_count = element_expr[2];
      new_state.repeat_count++;
      // The repeat rule can be completed, and we advance the state. Don't forget to
      // reset the repeat count.
      if (new_state.repeat_count >= min_repeat_count) {
        Enqueue(ParserState{
            parent_state.rule_id,
            parent_state.sequence_id,
            parent_state.element_id + 1,
            parent_state.rule_start_pos,
            0
        });
      }
      // If the repeat count is less than the max repeat count, we can continue to
      // visit the repeat state for another round.
      if (new_state.repeat_count < max_repeat_count) {
        Enqueue(new_state);
      }
      continue;
    }
    // If the rule is referenced by a fsm, we need to advance the fsm.
    XGRAMMAR_DCHECK(grammar_->per_rule_fsms[parent_state.rule_id].has_value());
    const auto& current_fsm = grammar_->per_rule_fsms[parent_state.rule_id].value();
    for (const auto edge : current_fsm->GetEdges(parent_state.element_id)) {
      if (edge.IsRuleRef() && edge.GetRefRuleId() == state.rule_id) {
        Enqueue(
            {parent_state.rule_id,
             parent_state.sequence_id,
             edge.target,
             parent_state.rule_start_pos,
             0}
        );
      }
    }
  }
}

std::pair</* scanable */ bool, /* completable */ bool> EarleyParser::Predict(
    const ParserState& state, const GrammarExpr& grammar_expr
) {
  // Check if it's the tag dispatch.
  if (state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value()) {
    // Try to expand the fsm.
    ExpandNextRuleRefElement(state, grammar_expr, nullptr);
    const auto& fsm = grammar_->per_rule_fsms[state.rule_id].value();
    return std::make_pair(fsm.IsScanableState(state.element_id), fsm.IsEndState(state.element_id));
  }
  XGRAMMAR_DCHECK(
      grammar_expr.type == GrammarExprType::kSequence ||
      grammar_expr.type == GrammarExprType::kEmptyStr
  );
  if (state.element_id == grammar_expr.size()) {
    // The rule is completed.
    return std::make_pair(false, true);
  }
  const auto& element_expr = grammar_->GetGrammarExpr(grammar_expr[state.element_id]);
  switch (element_expr.type) {
    case GrammarExprType::kRuleRef: {
      ExpandNextRuleRefElement(state, grammar_expr, &element_expr);
      return std::make_pair(false, false);
    }
    case GrammarExprType::kCharacterClassStar: {
      if (state.sub_element_id == 0) {
        Enqueue(ParserState{
            state.rule_id, state.sequence_id, state.element_id + 1, state.rule_start_pos, 0
        });
      }
      return std::make_pair(true, false);
    }
    case GrammarExprType::kRepeat: {
      const int32_t& min_repeat_count = element_expr[1];
      const int32_t& max_repeat_count = element_expr[2];
      // If the current repeat count is less than the max repeat count,
      // we can expand the next rule reference element.
      XGRAMMAR_DCHECK(state.repeat_count <= max_repeat_count);
      ExpandNextRuleRefElement(state, grammar_expr, &element_expr);
      if (state.repeat_count >= min_repeat_count) {
        Enqueue(ParserState{
            state.rule_id, state.sequence_id, state.element_id + 1, state.rule_start_pos, 0
        });
      }
      return std::make_pair(false, false);
    }
    default: {
      return std::make_pair(true, false);
    }
  }
}

void EarleyParser::Scan(const ParserState& state, const uint8_t ch) {
  const auto& cur_rule = grammar_->GetGrammarExpr(state.sequence_id);

  if (state.rule_id == -1 || (!grammar_->per_rule_fsms[state.rule_id].has_value())) {
    const auto& element_expr = grammar_->GetGrammarExpr(cur_rule[state.element_id]);
    // The element is a rule reference, we do not need to scan it.
    switch (element_expr.type) {
      case (GrammarExprType::kByteString): {
        AdvanceByteString(state, ch, element_expr);
        break;
      }
      case (GrammarExprType::kCharacterClass): {
        AdvanceCharacterClass(state, ch, element_expr);
        break;
      }
      case (GrammarExprType::kCharacterClassStar): {
        AdvanceCharacterClassStar(state, ch, element_expr);
        break;
      }
      default: {
        XGRAMMAR_LOG(FATAL) << "The element type is not supported! The type is: "
                            << int(element_expr.type);
      }
    }
  } else {
    AdvanceFsm(state, ch, cur_rule);
  }
}

/*!
  \note The workflow of Advance is as follows:
  1. Scan all the states in the latest states. Add all the possible states
  to the next states.
  2. If the next states are empty, then the character is not accepted.
  3. If the next states are not empty, then the character is accepted. Moreover,
  we need to complete and predict the next states.

  \note Thus, when initializing the Earley parser, we need to add the initial state
  to the history_states[0], and perform prediction and completion on the initial state.
*/
bool EarleyParser::Advance(const uint8_t ch) {
  // Initialize the containers.
  XGRAMMAR_DCHECK(tmp_process_state_queue_.empty())
      << "The tmp_process_state_queue_ should be empty before the scan.";
  tmp_states_visited_in_queue_.Clear();
  tmp_states_to_be_added_.clear();
  tmp_accept_stop_token_ = false;
  const auto& latest_states = scanable_state_history_[scanable_state_history_.size() - 1];
  // Scan all the scanable states.
  for (const auto& state : latest_states) {
    Scan(state, ch);
  }

  // Check if the character is accepted.
  if (tmp_process_state_queue_.empty() && tmp_states_to_be_added_.empty()) {
    return false;
  }

  // execute Predict and Complete for all states in the queue until empty.
  rule_id_to_completeable_states_.emplace_back();
  while (!tmp_process_state_queue_.empty()) {
    const auto state = tmp_process_state_queue_.front();
    tmp_process_state_queue_.pop();
    GrammarExpr grammar_expr = grammar_->GetGrammarExpr(state.sequence_id);
    auto [scanable, completable] = Predict(state, grammar_expr);
    if (completable) {
      Complete(state, grammar_expr);
    }
    if (scanable) {
      tmp_states_to_be_added_.push_back(state);
    }
  }

  // Check if the grammar is completed, and add the scannable states to the history.
  is_completed_.push_back(tmp_accept_stop_token_);
  scanable_state_history_.PushBack(tmp_states_to_be_added_);
  return true;
}

EarleyParser::EarleyParser(
    const Grammar& grammar, const ParserState& init_state, const bool need_expand
)
    : grammar_(grammar) {
  // Check if the initial state is valid. If invalid, then we choose the root state as default.
  ParserState init = init_state;
  if (init_state.IsInvalid()) {
    init = ParserState(
        grammar_->GetRootRuleId(),
        ParserState::kUnexpandedRuleStartSequenceId,
        0,
        ParserState::kNoPrevInputPos,
        0
    );
  } else {
    init = init_state;
  }

  // If there is no need to expand the initial state, we only need to add it to the
  // scanable states history.
  if (!need_expand) {
    rule_id_to_completeable_states_.emplace_back();
    is_completed_.push_back(false);
    scanable_state_history_.PushBack({init});
  }

  // Otherwise, we expand the initial state, and process the queue.
  PushStateAndExpand(init);
}

void EarleyParser::PushStateAndExpand(const ParserState& state) {
  tmp_states_visited_in_queue_.Clear();
  tmp_accept_stop_token_ = false;
  tmp_states_to_be_added_.clear();
  rule_id_to_completeable_states_.emplace_back();
  if (state.IsInvalid()) {
    ExpandAndEnqueueUnexpandedState(ParserState{
        grammar_->GetRootRuleId(),
        ParserState::kUnexpandedRuleStartSequenceId,
        0,
        ParserState::kNoPrevInputPos,
        0
    });
  } else {
    // If the rule can't be expanded, we need to add it to the queue.
    if (!ExpandAndEnqueueUnexpandedState(state)) {
      Enqueue(state);
    }
  }
  while (!tmp_process_state_queue_.empty()) {
    const auto state = tmp_process_state_queue_.front();
    tmp_process_state_queue_.pop();
    GrammarExpr grammar_expr = grammar_->GetGrammarExpr(state.sequence_id);
    auto [scanable, completable] = Predict(state, grammar_expr);
    if (completable) {
      Complete(state, grammar_expr);
    }
    if (scanable) {
      tmp_states_to_be_added_.push_back(state);
    }
  }
  is_completed_.push_back(tmp_accept_stop_token_);
  scanable_state_history_.PushBack(tmp_states_to_be_added_);
}

void EarleyParser::Reset() {
  rule_id_to_completeable_states_.clear();
  scanable_state_history_.PopBack(scanable_state_history_.size());
  is_completed_.clear();
  stop_token_is_accepted_ = false;
  XGRAMMAR_DCHECK(tmp_process_state_queue_.empty());
  PushStateAndExpand(ParserState(
      grammar_->GetRootRuleId(),
      ParserState::kUnexpandedRuleStartSequenceId,
      0,
      ParserState::kNoPrevInputPos,
      0
  ));
}

bool EarleyParser::ExpandAndEnqueueUnexpandedState(const ParserState& state) {
  if (state.sequence_id != ParserState::kUnexpandedRuleStartSequenceId) {
    return false;
  }
  // The rule is already expanded, and finished.
  auto cur_rule_id = state.rule_id;
  auto cur_rule_body_id = grammar_->GetRule(cur_rule_id).body_expr_id;
  auto cur_rule_body = grammar_->GetGrammarExpr(cur_rule_body_id);
  // There are two types of an unexpanded rule:
  // 1. The rule is a tag dispatch rule.
  // 2. The rule is a choice, consisting of multiple sequences.
  if (state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value()) {
    Enqueue(ParserState{
        cur_rule_id,
        cur_rule_body_id,
        grammar_->per_rule_fsms[state.rule_id].value().GetStart(),
        ParserState::kNoPrevInputPos,
        0
    });
    return true;
  }
  XGRAMMAR_DCHECK(cur_rule_body.type == GrammarExprType::kChoices);
  for (const auto& sequence_id : cur_rule_body) {
    Enqueue(ParserState{cur_rule_id, sequence_id, 0, ParserState::kNoPrevInputPos, 0});
  }
  return true;
}

void EarleyParser::ExpandNextRuleRefElement(
    const ParserState& state, const GrammarExpr& grammar_expr, const GrammarExpr* sub_grammar_expr
) {
  std::vector<int32_t> ref_rule_ids;
  // Path A. The rule has a corresponding FSM.
  if (state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value()) {
    const auto& current_fsm = grammar_->per_rule_fsms[state.rule_id].value();
    for (const auto& edge : current_fsm->GetEdges(state.element_id)) {
      if (edge.IsRuleRef()) {
        ref_rule_ids.push_back(edge.GetRefRuleId());
      } else if (edge.IsEpsilon()) {
        Enqueue(ParserState{state.rule_id, state.sequence_id, edge.target, state.rule_start_pos, 0}
        );
      }
    }
  } else {
    XGRAMMAR_DCHECK(grammar_expr.type == GrammarExprType::kSequence);
    XGRAMMAR_DCHECK(
        sub_grammar_expr->type == GrammarExprType::kRuleRef ||
        sub_grammar_expr->type == GrammarExprType::kRepeat
    );
    ref_rule_ids.push_back((*sub_grammar_expr)[0]);
  }
  for (const auto& ref_rule_id : ref_rule_ids) {
    {  // Add the reference rule to map.
      if ((state.element_id != grammar_expr.size() - 1) ||
          state.rule_start_pos == ParserState::kNoPrevInputPos ||
          (state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value()) ||
          sub_grammar_expr->type == GrammarExprType::kRepeat) {
        // It's not the right recursion, or it's the root rule.
        auto& states_map = rule_id_to_completeable_states_.back();
        states_map.insert({ref_rule_id, state});
      } else {
        // If it's the right recursion, we need to add the ancestors of the parent state.
        auto& states_map = rule_id_to_completeable_states_.back();
        auto& parent_states_map = rule_id_to_completeable_states_[state.rule_start_pos];
        const auto& range = states_map.equal_range(ref_rule_id);
        const auto in_vec = [&](const ParserState& state_) {
          return std::find_if(range.first, range.second, [&](const auto& s) {
                   return StateEqualForParsing()(s.second, state_);
                 }) != range.second;
        };

        bool no_fsm_parent = true;
        for (auto parent_state_iter = parent_states_map.lower_bound(state.rule_id);
             parent_state_iter != parent_states_map.end() &&
             parent_state_iter->first == state.rule_id;
             parent_state_iter++) {
          const auto& parent_state = parent_state_iter->second;
          XGRAMMAR_DCHECK(ref_rule_id != -1);
          if (grammar_->per_rule_fsms[parent_state.rule_id].has_value()) {
            no_fsm_parent = false;
            break;
          }
        }
        if (no_fsm_parent) {
          for (auto parent_state_iter = parent_states_map.lower_bound(state.rule_id);
               parent_state_iter != parent_states_map.end() &&
               parent_state_iter->first == state.rule_id;
               parent_state_iter++) {
            const auto& parent_state = parent_state_iter->second;
            if (!in_vec(parent_state)) {
              states_map.insert({ref_rule_id, parent_state});
            }
          }
        } else {
          auto& states_map = rule_id_to_completeable_states_.back();
          states_map.insert({ref_rule_id, state});
        }
      }

      // Check if the reference rule is already visited.
      if (IsStateVisitedInQueue({ref_rule_id, -1, -1, -1, -1})) {
        if (std::binary_search(
                grammar_->allow_empty_rule_ids.begin(),
                grammar_->allow_empty_rule_ids.end(),
                ref_rule_id
            )) {
          if (state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value()) {
            const auto& current_fsm = grammar_->per_rule_fsms[state.rule_id].value();
            for (const auto& edge : current_fsm->GetEdges(state.element_id)) {
              if (edge.IsRuleRef() && edge.GetRefRuleId() == ref_rule_id) {
                Enqueue(ParserState{
                    state.rule_id, state.sequence_id, edge.target, state.rule_start_pos, 0
                });
              }
            }
            continue;
          }
          XGRAMMAR_DCHECK(grammar_expr.type == GrammarExprType::kSequence);
          if (sub_grammar_expr->type == GrammarExprType::kRuleRef) {
            Enqueue(ParserState{
                state.rule_id, state.sequence_id, state.element_id + 1, state.rule_start_pos, 0
            });
          }
        }
        continue;
      }

      // If the reference rule is not visited, we need to add it to the queue.
      tmp_states_visited_in_queue_.Insert({ref_rule_id, -1, -1, -1, -1});
      const auto& ref_rule = grammar_->GetRule(ref_rule_id);
      const auto& ref_grammar_expr_id = ref_rule.body_expr_id;
      const auto& ref_grammar_expr = grammar_->GetGrammarExpr(ref_grammar_expr_id);
      XGRAMMAR_DCHECK(
          ref_grammar_expr.type == GrammarExprType::kChoices ||
          grammar_->per_rule_fsms[ref_rule_id].has_value()
      );
      if (grammar_->per_rule_fsms[ref_rule_id].has_value()) {
        tmp_process_state_queue_.push(ParserState{
            ref_rule_id,
            ref_grammar_expr_id,
            grammar_->per_rule_fsms[ref_rule_id]->GetStart(),
            int32_t(rule_id_to_completeable_states_.size()) - 1,
            0
        });
      } else {
        for (const auto& sequence_id : ref_grammar_expr) {
          const auto& sequence = grammar_->GetGrammarExpr(sequence_id);
          if (sequence.type == GrammarExprType::kEmptyStr) {
            if (state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value()) {
              const auto& current_fsm = grammar_->per_rule_fsms[state.rule_id].value();
              XGRAMMAR_DCHECK(state.element_id < current_fsm->NumStates());
              const auto& current_edges = current_fsm->GetEdges(state.element_id);
              for (const auto& edge : current_edges) {
                if (edge.IsRuleRef() && edge.GetRefRuleId() == ref_rule_id) {
                  Enqueue(ParserState{
                      state.rule_id, state.sequence_id, edge.target, state.rule_start_pos, 0
                  });
                }
              }
            } else {
              Enqueue(ParserState{
                  state.rule_id, state.sequence_id, state.element_id + 1, state.rule_start_pos, 0
              });
            }
            continue;
          }
          // Assert: the state can't be repeated. Since the rule_start_pos is the current
          // position, and the rule can only be predicted once.
          tmp_process_state_queue_.push(ParserState{
              ref_rule_id, sequence_id, 0, int32_t(rule_id_to_completeable_states_.size()) - 1, 0
          });
        }
      }
    }
  }
}

void EarleyParser::AdvanceByteString(
    const ParserState& state, const uint8_t ch, const GrammarExpr& sub_rule
) {
  XGRAMMAR_DCHECK(sub_rule.type == GrammarExprType::kByteString);
  XGRAMMAR_DCHECK(sub_rule.size() > state.sub_element_id);
  if (static_cast<uint8_t>(sub_rule[state.sub_element_id]) == ch) {
    auto new_state = state;
    new_state.sub_element_id++;
    if (new_state.sub_element_id == sub_rule.size()) {
      new_state.element_id++;
      new_state.sub_element_id = 0;
      Enqueue(new_state);
      // Assert: In a sequence, the bytestring can't be skipped. So the state can't be repeated.
    } else {
      tmp_states_to_be_added_.push_back(new_state);
    }
  }
  return;
}

void EarleyParser::AdvanceCharacterClass(
    const ParserState& state, const uint8_t ch, const GrammarExpr& sub_sequence
) {
  XGRAMMAR_DCHECK(sub_sequence.type == GrammarExprType::kCharacterClass)
      << "The element type is not supported!";

  // The state is matching a UTF8 character.
  if (state.sub_element_id > 0) {
    if ((ch & 0xC0) == 0x80) {
      auto new_state = state;
      new_state.sub_element_id--;
      // Check if the UTF8 character is completed.
      if (new_state.sub_element_id == 0) {
        new_state.element_id++;
        Enqueue(new_state);
        // Assert: In a sequence, the CharacterClass can't be skipped. So the state can't be
        // repeated. the fllowing tmp_process_state_queue_.push(new_state) is for the same reason.
      } else {
        tmp_states_to_be_added_.push_back(new_state);
      }
    }
    return;
  }
  bool is_negative = static_cast<bool>(sub_sequence[0]);

  // This trick is based on the current structure that character class
  // can't accept a UTF8 character, unless it has a negation.
  if (!isascii(ch)) {
    if (!is_negative) {
      return;
    }
    auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
    if (!accepted) {
      return;
    }

    // A new UTF8 character is accepted.
    XGRAMMAR_DCHECK(num_bytes > 1);
    auto new_state = state;
    new_state.sub_element_id = num_bytes - 1;
    tmp_states_to_be_added_.push_back(new_state);
    return;
  }

  for (int i = 1; i < sub_sequence.size(); i += 2) {
    if (static_cast<uint8_t>(sub_sequence[i]) <= ch &&
        ch <= static_cast<uint8_t>(sub_sequence[i + 1])) {
      if (!is_negative) {
        auto new_state = state;
        new_state.element_id++;
        new_state.sub_element_id = 0;
        Enqueue(new_state);
      }
      return;
    }
  }
  if (is_negative) {
    auto new_state = state;
    new_state.element_id++;
    new_state.sub_element_id = 0;
    Enqueue(new_state);
  }
}

void EarleyParser::AdvanceCharacterClassStar(
    const ParserState& state, const uint8_t ch, const GrammarExpr& sub_sequence
) {
  XGRAMMAR_DCHECK(sub_sequence.type == GrammarExprType::kCharacterClassStar)
      << "The element type is not supported!";

  // The state is matching a UTF8 character.
  if (state.sub_element_id > 0) {
    if ((ch & 0xC0) == 0x80) {
      auto new_state = state;
      new_state.sub_element_id--;
      // Check if the UTF8 character is completed.
      if (new_state.sub_element_id == 0) {
        Enqueue(new_state);
      } else {
        tmp_states_to_be_added_.push_back(new_state);
      }
    }
    return;
  }
  bool is_negative = static_cast<bool>(sub_sequence[0]);

  // This trick is based on the current structure that character class
  // can't accept a UTF8 character, unless it has a negation.
  if (!isascii(ch)) {
    if (!is_negative) {
      return;
    }
    auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
    if (!accepted) {
      return;
    }
    // A new UTF8 character is accepted.
    XGRAMMAR_DCHECK(num_bytes > 1);
    auto new_state = state;
    new_state.sub_element_id = num_bytes - 1;
    tmp_states_to_be_added_.push_back(new_state);
    return;
  }

  for (int i = 1; i < sub_sequence.size(); i += 2) {
    if (static_cast<uint8_t>(sub_sequence[i]) <= ch &&
        ch <= static_cast<uint8_t>(sub_sequence[i + 1])) {
      if (!is_negative) {
        Enqueue(state);
      }
      return;
    }
  }
  if (is_negative) {
    Enqueue(state);
  }
}

void EarleyParser::AdvanceFsm(
    const ParserState& state, const uint8_t ch, const GrammarExpr& cur_sequence
) {
  XGRAMMAR_DCHECK(state.rule_id != -1 && grammar_->per_rule_fsms[state.rule_id].has_value());
  auto current_fsm = grammar_->per_rule_fsms[state.rule_id].value();
  current_fsm->GetNextStates(
      state.element_id, ch, FSMEdge::EdgeType::kCharRange, &tmp_fsm_targets_
  );
  for (const auto& next_node : tmp_fsm_targets_) {
    auto new_state = state;
    new_state.element_id = next_node;
    tmp_process_state_queue_.push(new_state);
  }
}

bool RepeatDetector::IsVisited(const ParserState& state) const {
  // If the size is larger than the threshold, then we use the set to check.
  if (size_ > transition_threshold_) {
    return visited_set_.find(state) != visited_set_.end();
  }
  return std::find_if(
             visited_vector_.begin(),
             visited_vector_.begin() + size_,
             [&state](const ParserState& s) { return StateEqualForParsing()(state, s); }
         ) != visited_vector_.begin() + size_;
}

void RepeatDetector::Insert(const ParserState& state) {
  if (size_ == transition_threshold_) {
    for (const auto& s : visited_vector_) {
      visited_set_.insert(s);
    }
  }
  size_++;
  if (size_ > transition_threshold_) {
    visited_set_.insert(state);
  } else {
    visited_vector_[size_ - 1] = state;
  }
}

void RepeatDetector::Clear() {
  if (size_ > transition_threshold_) {
    visited_set_.clear();
  }
  size_ = 0;
}

}  // namespace xgrammar
