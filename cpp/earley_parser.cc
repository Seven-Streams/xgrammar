#include "earley_parser.h"

#include <cstdint>

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
}  // namespace xgrammar
