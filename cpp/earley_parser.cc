#include "earley_parser.h"

#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "grammar_data_structure.h"
#include "support/logging.h"
#include "xgrammar/grammar.h"
namespace xgrammar {
constexpr int32_t kUnexpandedRuleStartSequenceId = 128000;

// constexpr int32_t kDispatchedTagDispatchElementId = -1;
using RuleExprType = Grammar::Impl::RuleExprType;
using RuleExpr = Grammar::Impl::RuleExpr;
inline bool EarleyParser::IsEndOfGrammar(const State& state) const {
  if (state.parent_pos != State::kNoParent) {
    return false;
  }
  auto seq_expr = grammar_->GetRuleExpr(state.sequence_id);
  if (seq_expr.type == RuleExprType::kTagDispatch) {
    return state.element_id != -1;
  } else {
    return seq_expr.size() == state.element_id;
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
    State new_state = State(
        parent_state.rule_id,
        parent_state.sequence_id,
        parent_state.element_id + 1,
        parent_state.left_utf8_bytes,
        parent_state.element_in_string,
        parent_state.parent_pos
    );
    in_vec = std::find(history_states.back().begin(), history_states.back().end(), new_state) !=
             history_states.back().end();
    if (!in_vec) {
      history_states.back().push_back(new_state);
      queue.push(new_state);
    }
  }
  return;
}
inline void EarleyParser::Predict(const State& state) {
  auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
  auto crrent_element = grammar_->GetRuleExpr(cur_rule[state.element_id]);
  switch (cur_rule.type) {
    case RuleExprType::kByteString:
    case RuleExprType::kCharacterClass:
    case RuleExprType::kCharacterClassStar:
    case RuleExprType::kEmptyStr:
      return;
    // These types will never predict other new rules.
    case RuleExprType::kRuleRef: {
      const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
      bool in_vec = ptr != states.back().end();
      if (in_vec) {
        ptr->predictions->push_back(crrent_element[0]);
      } else {
        states.back().push_back(state);
        states.back().back().predictions = std::vector<int32_t>();
        states.back().back().predictions->push_back(crrent_element[0]);
      }
      history_states.back().emplace_back(
          crrent_element[0], kUnexpandedRuleStartSequenceId, 0, state.rule_id
      );
      break;
    }
    case RuleExprType::kSequence: {
      break;
    }
    case RuleExprType::kChoices: {
      break;
    }
    case RuleExprType::kTagDispatch: {
      break;
    }
  }
}
inline void MoveToNextPosition(State& state) {
  state.element_id += 1;
  state.element_in_string = 0;
  state.left_utf8_bytes = 0;
  return;
}
inline void EarleyParser::Scan(const State& state, const uint8_t& ch) {
  // TODO:
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
  for (const auto& state : history_states.back()) {
    queue.push(state);
  }
  while (!queue.empty()) {
    const auto& state = queue.front();
    if (visited.find(queue.back()) != visited.end()) {
      queue.pop();
      continue;
    }
    visited.insert(state);
    Complete(state);
    Predict(state);
    queue.pop();
  }
  return true;
}

EarleyParser::EarleyParser(const Grammar& grammar) : grammar_(grammar) {
  states.push_back(std::vector<State>());
  history_states.push_back(std::vector<State>());
  history_states.back().emplace_back(
      grammar_->GetRootRuleId(), kUnexpandedRuleStartSequenceId, 0, State::kNoParent
  );
  std::unordered_set<State, StateHash> visited;
  queue.push(history_states[0][0]);
  while (!queue.empty()) {
    const auto& state = queue.front();
    if (visited.find(queue.back()) != visited.end()) {
      queue.pop();
      continue;
    }
    visited.insert(state);
    Complete(state);
    Predict(state);
    queue.pop();
  }
}
}  // namespace xgrammar
