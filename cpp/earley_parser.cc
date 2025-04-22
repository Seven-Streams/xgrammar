#include "earley_parser.h"

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "grammar_data_structure.h"
#include "support/logging.h"
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

inline void EarleyParser::Complete(const State& state) const {
  // TODO:
}
inline void EarleyParser::Predict(const State& state) const {
  // TODO:
}
inline bool EarleyParser::Scan(const State& state, const uint8_t& ch) const {
  // TODO:
  XGRAMMAR_LOG(FATAL) << "Scan is not implemented yet.";
}

bool EarleyParser::Advance(const uint8_t& ch) {
  current_states.clear();
  std::vector<State> state_queue = states.back();
  history_states.push_back(std::vector<State>());
  states.push_back(std::vector<State>());
  while (!state_queue.empty()) {
    auto state = state_queue.back();
    state_queue.pop_back();
    if (current_states.find(state) != current_states.end()) {
      continue;
    }
    current_states.insert(state);
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
