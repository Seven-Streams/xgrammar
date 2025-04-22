#include "earley_parser.h"

#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "grammar_data_structure.h"
#include "support/encoding.h"
#include "support/logging.h"
#include "xgrammar/grammar.h"
namespace xgrammar {
constexpr int32_t kUnexpandedRuleStartSequenceId = 128000;

constexpr int32_t kDispatchedTagDispatchElementId = -1;
using RuleExprType = Grammar::Impl::RuleExprType;
using RuleExpr = Grammar::Impl::RuleExpr;
inline bool EarleyParser::IsEndOfGrammar(const State& stack_element) const {
  if (stack_element.parent_pos != State::kNoParent) {
    return false;
  }
  auto seq_expr = grammar_->GetRuleExpr(stack_element.sequence_id);
  if (seq_expr.type == RuleExprType::kTagDispatch) {
    return stack_element.element_id != -1;
  } else {
    return seq_expr.size() == stack_element.element_id;
  }
}

bool EarleyParser::CanReachEnd() const {
  const auto& current_states = history_states.back();
  return std::any_of(current_states.begin(), current_states.end(), [&](const State& state) {
    return IsEndOfGrammar(state);
  });
}

void EarleyParser::PopBackStates(int32_t cnt) {
  if (cnt >= static_cast<int32_t>(states.size())) {
    XGRAMMAR_LOG(FATAL) << "The number of states to be popped is larger than the size of states.";
  }
  states.erase(states.end() - cnt, states.end());
  history_states.erase(history_states.end() - cnt, history_states.end());
  return;
}

inline void EarleyParser::Complete(const State& state) {
  auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
  /*
    Case 1: The current state is the end of the rule.
    Case 2: The current rule can be empty.
    If neither of the above two cases is true, we return.
  */
  if (((cur_rule.size() != state.element_id) || (state.parent_pos == State::kNoParent)) &&
      (std::find(
           grammar_->allow_empty_rule_ids.begin(),
           grammar_->allow_empty_rule_ids.end(),
           state.rule_id
       ) == grammar_->allow_empty_rule_ids.end())) {
    return;
  }
  const auto& parent_states_list = states[state.parent_pos];
  for (const auto& parent_state : parent_states_list) {
    if (!parent_state.predictions.has_value()) {
      continue;
    }
    bool in_vec =
        std::find(
            parent_state.predictions->begin(), parent_state.predictions->end(), state.rule_id
        ) != parent_state.predictions->end();
    if (!in_vec) {
      continue;
    }
    history_states.back().emplace_back(
        parent_state.rule_id,
        parent_state.sequence_id,
        parent_state.element_id + 1,
        parent_state.left_utf8_bytes,
        parent_state.element_in_string,
        parent_state.parent_pos
    );
  }
  return;
}
inline void EarleyParser::Predict(const State& state) {
  // Step 1. Handle unexpanded rules.
  if (state.sequence_id == kUnexpandedRuleStartSequenceId) {
    auto cur_rule_id = state.rule_id;
    auto cur_rule_body_id = grammar_->GetRule(cur_rule_id).body_expr_id;
    auto cur_rule_body = grammar_->GetRuleExpr(cur_rule_body_id);

    if (cur_rule_body.type == RuleExprType::kTagDispatch) {
      history_states.back().emplace_back(
          cur_rule_id,
          cur_rule_body_id,
          grammar_->root_tag_dispatch_fsm->StartNode(),
          state.parent_pos
      );
      return;
    } else {
      XGRAMMAR_DCHECK(cur_rule_body.type == RuleExprType::kChoices);
      for (auto sequence_id : cur_rule_body) {
        auto ref_rule_sequence = grammar_->GetRuleExpr(sequence_id);
        if (ref_rule_sequence.type == RuleExprType::kEmptyStr &&
            state.parent_pos != State::kNoParent) {
          // If the empty string is in a root rule, it indicates the end of the grammar and we
          // just add it as a stack top to indicate the matching ends.
          continue;
        }
        history_states.back().emplace_back(cur_rule_id, sequence_id, 0, state.parent_pos);
      }
      return;
    }
  }

  auto cur_sequence = grammar_->GetRuleExpr(state.sequence_id);

  auto current_element = grammar_->GetRuleExpr(cur_sequence[state.element_id]);

  // Step 3. Iterate into sub rules
  if (current_element.type == RuleExprType::kRuleRef) {
    State sub_rule =
        State(current_element[0], kUnexpandedRuleStartSequenceId, 0, states.size() - 1);
    history_states.back().push_back(sub_rule);
    auto parent = std::find(states.back().begin(), states.back().end(), state);
    if (parent == states.back().end()) {
      auto state_copy = state;
      state_copy.predictions = std::vector<int32_t>({current_element[0]});
    } else {
      parent->predictions->push_back(current_element[0]);
    }
  }
  return;
}
inline bool EarleyParser::IsAccepted(const State& state, const uint8_t& ch) const {
  auto current_sequence = grammar_->GetRuleExpr(state.sequence_id);
  if (current_sequence.type == RuleExprType::kTagDispatch) {
    XGRAMMAR_DCHECK(state.element_id != -1);
    return true;
  }

  auto current_element = grammar_->GetRuleExpr(current_sequence[state.element_id]);
  if (current_element.type == RuleExprType::kCharacterClass ||
      current_element.type == RuleExprType::kCharacterClassStar) {
    if (state.left_utf8_bytes > 0) {
      return (ch & 0xC0) == 0x80;
    }
    auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
    if (!accepted) {
      return false;
    }
    bool is_negative = static_cast<bool>(current_element[0]);
    if (num_bytes > 1) {
      return is_negative;
    }
    for (int i = 1; i < current_element.size(); i += 2) {
      if (current_element[i] <= ch && ch <= current_element[i + 1]) {
        return !is_negative;
      }
    }
    return is_negative;
  } else if (current_element.type == RuleExprType::kByteString) {
    return current_element[state.element_in_string] == ch;
  } else {
    XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in CheckIfAccepted: "
                        << static_cast<int>(current_element.type);
  }
}
inline void MoveToNextPosition(State& state) {
  state.element_id += 1;
  state.element_in_string = 0;
  state.left_utf8_bytes = 0;
  return;
}
inline void EarleyParser::Scan(const State& state, const uint8_t& ch) {
  if (IsAccepted(state, ch)) {
    auto current_sequence = grammar_->GetRuleExpr(state.sequence_id);
    if (current_sequence.type == Grammar::Impl::RuleExprType::kTagDispatch) {
      auto root_tag_dispatch_fsm = grammar_->root_tag_dispatch_fsm;
      if (!root_tag_dispatch_fsm) {
        XGRAMMAR_LOG(FATAL
        ) << "The grammar does not have a root tag dispatch rule; it is not built.";
        XGRAMMAR_UNREACHABLE();
      }
      auto start_node = root_tag_dispatch_fsm->StartNode();
      auto next_node = root_tag_dispatch_fsm->Transition(state.element_id, ch);
      auto new_state = state;
      if (next_node == CompactFSM::NO_TRANSITION) {
        // Case 1. The new char cannot continue to be accepted by the tag dispatch fsm.
        // We try to accept the new char from the start node. If accepted, we go to the target node.
        // If it still cannot be accepted, we stay at the start node.
        auto new_next_node = root_tag_dispatch_fsm->Transition(start_node, ch);
        new_state.element_id =
            new_next_node == CompactFSM::NO_TRANSITION ? start_node : new_next_node;
      } else if (!root_tag_dispatch_fsm->IsEndNode(next_node)) {
        // Case 2. The new char can continue to be accepted by the tag dispatch fsm.
        // We need to update the element id to the next node.
        new_state.element_id = next_node;
      } else {
        // Case 3. The new char can continue to be accepted by the tag dispatch fsm.
        // We need to dispatch the tag dispatch fsm to the end node.
        // We need to create a new stack element to represent the dispatched tag dispatch.
        new_state.element_id = kDispatchedTagDispatchElementId;
        XGRAMMAR_DCHECK(grammar_->tag_dispatch_end_node_to_rule_id.count(next_node))
            << "The end node of the tag dispatch fsm does not correspond to any rule id";
        auto refered_rule_id = grammar_->tag_dispatch_end_node_to_rule_id.at(next_node);
        history_states.back().emplace_back(
            refered_rule_id, kUnexpandedRuleStartSequenceId, 0, state.parent_pos
        );
        return;
      }
    }

    auto current_element = grammar_->GetRuleExpr(current_sequence[state.element_id]);
    State new_state = State(
        state.rule_id,
        state.sequence_id,
        state.element_id,
        state.left_utf8_bytes,
        state.element_in_string,
        state.parent_pos
    );
    switch (current_element.type) {
      case RuleExprType::kCharacterClass: {
        if (state.left_utf8_bytes > 1) {
          new_state.left_utf8_bytes -= 1;
          history_states.back().emplace_back(new_state);
          return;
        } else if (state.left_utf8_bytes == 1) {
          MoveToNextPosition(new_state);
          history_states.back().emplace_back(new_state);
          return;
        }
        // If no left utf8 bytes, check the first byte to find the left bytes needed.
        XGRAMMAR_DCHECK(state.left_utf8_bytes == 0);
        auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
        XGRAMMAR_DCHECK(accepted);
        if (num_bytes > 1) {
          new_state.left_utf8_bytes = num_bytes - 1;
          history_states.back().emplace_back(new_state);
          return;
        }
        MoveToNextPosition(new_state);
        history_states.back().emplace_back(new_state);
        return;
      }
      case RuleExprType::kCharacterClassStar: {
        if (state.left_utf8_bytes >= 1) {
          new_state.left_utf8_bytes -= 1;
        } else {
          XGRAMMAR_DCHECK(state.left_utf8_bytes == 0);
          auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
          XGRAMMAR_DCHECK(accepted);
          new_state.left_utf8_bytes = num_bytes - 1;
        }
        history_states.back().emplace_back(new_state);
        return;
      }
      case RuleExprType::kByteString: {
        if (state.element_in_string + 1 < current_element.size()) {
          new_state.element_in_string += 1;
          history_states.back().emplace_back(new_state);
          return;
        }
        MoveToNextPosition(new_state);
        history_states.back().emplace_back(new_state);
        return;
      }
      default:
        XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in AdvanceStackElementWithChar: "
                            << static_cast<int>(current_element.type);
    }
  }
  return;
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
bool EarleyParser::Advance(const uint8_t& ch) {
  const auto& latest_states = history_states.back();
  history_states.push_back(std::vector<State>());
  for (const auto& state : latest_states) {
    Scan(state, ch);
  }
  if (history_states.back().empty()) {
    history_states.pop_back();
    return false;
  }
  states.push_back(std::vector<State>());
  std::unordered_set<State, StateHash> visited;
  // We need a copy of the states, since we need to rollback the states.
  auto check_queue = history_states.back();
  while (!check_queue.empty()) {
    const auto& state = check_queue.back();
    if (visited.find(check_queue.back()) != visited.end()) {
      check_queue.pop_back();
      continue;
    }
    visited.insert(state);
    Complete(state);
    Predict(state);
    check_queue.pop_back();
  }
  return true;
}
}  // namespace xgrammar
