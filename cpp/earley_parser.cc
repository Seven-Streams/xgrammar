#include "earley_parser.h"

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "grammar_data_structure.h"
namespace xgrammar {

inline bool EarleyParser::IsEndOfGrammar(const State& stack_element) const {
  if (stack_element.parent_id != State::kNoParent) {
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
  const auto& current_states = states.back();
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

inline void EarleyParser::Complete(const State& state) const {}
bool EarleyParser::Advance(const uint8_t& ch) {
  std::vector<State> state_queue = states.back();
  class StateHash {
   public:
    size_t operator()(const State& state) const {
      return std::hash<int32_t>()(state.rule_id) ^ std::hash<int32_t>()(state.sequence_id) ^
             std::hash<int32_t>()(state.element_id) ^ std::hash<int32_t>()(state.parent_id);
    }
  };
  std::unordered_set<State, StateHash> visited;
  history_states.push_back(std::vector<State>());
  states.push_back(std::vector<State>());
  while (!state_queue.empty()) {
    auto state = state_queue.back();
    state_queue.pop_back();
    if (visited.find(state) != visited.end()) {
      continue;
    }
    visited.insert(state);
    Predict(state);
    Complete(state);
    Scan(state, ch);
  }
  if (history_states.back().empty()) {
    history_states.pop_back();
    states.pop_back();
    return false;
  }
  return true;
}
}  // namespace xgrammar
