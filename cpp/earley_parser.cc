#include "earley_parser.h"

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "grammar_data_structure.h"
#include "support/logging.h"
namespace xgrammar {

inline bool EarleyParser::IsEndOfGrammar(const State& stack_element) const {
  if (stack_element.parent_pos != State::kNoParent) {
    return false;
  }
  auto seq_expr = grammar_->GetRuleExpr(stack_element.sequence_id);
  if (seq_expr.type == Grammar::Impl::RuleExprType::kTagDispatch) {
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
  if ((cur_rule.size() != state.element_id) || (state.parent_pos == State::kNoParent)) {
    return;
  }
  const auto& parent_states_list = states[state.parent_pos];
  for (const auto& parent_state : parent_states_list) {
    if (!parent_state.predictions.has_value()) {
      continue;
    }
    bool in_vec =
        std::find(parent_state.predictions->begin(), parent_state.predictions->end(), state) !=
        parent_state.predictions->end();
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
  // TODO:
}
inline bool EarleyParser::Scan(const State& state, const uint8_t& ch) {
  // TODO:
  XGRAMMAR_LOG(FATAL) << "Scan is not implemented yet.";
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
