#include "earley_parser.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <unordered_set>
#include <utility>
#include <vector>

#include "grammar_data_structure.h"
#include "support/encoding.h"
#include "support/logging.h"
#include "xgrammar/grammar.h"
namespace xgrammar {
constexpr int32_t kUnexpandedRuleStartSequenceId = 128000;
constexpr int32_t kUnexpandedRuleFinishElementId = 128000;
using RuleExprType = Grammar::Impl::RuleExprType;
using RuleExpr = Grammar::Impl::RuleExpr;

inline bool EarleyParser::IsEndOfGrammar(const State& state) const {
  // The root rule.
  if (state.parent_pos == State::kNoParent) {
    return true;
  }
  if (state.sequence_id == kUnexpandedRuleStartSequenceId) {
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
  can_reach_end.erase(can_reach_end.end() - cnt, can_reach_end.end());
  return;
}

inline void EarleyParser::Complete(const State& state) {
  if (state.sequence_id == kUnexpandedRuleStartSequenceId) {
    if (state.element_id != kUnexpandedRuleFinishElementId) {
      return;
    }
  } else {
    auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
    // Xgrammar_LOG(INFO) << "Complete: " << state << ", type is " << int(cur_rule.type)
    //   << ", the element size is " << cur_rule.size() << std::endl;
    /*
      Case 1: The current state is the end of the rule.
      Case 2: The current rule can be empty.
      Case 3: The type is kCharacterClassStar, and the element_id == 0.
      If none of the above cases is true, we return.
    */
    if ((!((cur_rule.size() == state.element_id) || (state.parent_pos == State::kNoParent))) &&
        (cur_rule.type != RuleExprType::kEmptyStr) &&
        (!((cur_rule.type == RuleExprType::kCharacterClassStar) && (state.element_id == 0)))) {
      return;
    }
  }
  // Check all the possible parent states.
  const auto& parent_states_list = states[state.parent_pos];
  for (const auto& parent_state : parent_states_list) {
    // XGRAMMAR_DCHECK(parent_state.predictions.has_value());
    bool in_vec = std::find(
                      parent_state.predictions->begin(),
                      parent_state.predictions->end(),
                      std::make_pair(state.rule_id, state.sequence_id)
                  ) != parent_state.predictions->end();
    if (!in_vec) {
      continue;
    }
    // The parent state indeed has a prediction for the current state.

    if (parent_state.sequence_id == kUnexpandedRuleStartSequenceId) {
      queue.emplace(State{
          parent_state.rule_id,
          parent_state.sequence_id,
          kUnexpandedRuleFinishElementId,
          parent_state.parent_pos
      });
      return;
    }
    auto parent_expr = grammar_->GetRuleExpr(parent_state.sequence_id);
    switch (parent_expr.type) {
      // These types can never predict other new rules, Thus
      // They will never be a legal parent state.
      case RuleExprType::kByteString:
      case RuleExprType::kCharacterClass:
      case RuleExprType::kCharacterClassStar:
      case RuleExprType::kEmptyStr:
        XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in Complete: "
                            << static_cast<int>(parent_expr.type);
      // These two types can predict other new rules. We need to
      // to move to the next element.
      case RuleExprType::kRuleRef:
      case RuleExprType::kSequence: {
        queue.emplace(State{
            parent_state.rule_id,
            parent_state.sequence_id,
            parent_state.element_id + 1,
            parent_state.parent_pos
        });
        break;
      }
      // These two types can predict other new rules, and have a
      // completion is enough to complete the parent state.
      // We need to move to the end of the element. i.e. complete the parent state.
      case RuleExprType::kChoices:
      case RuleExprType::kTagDispatch: {
        queue.emplace(State{
            parent_state.rule_id,
            parent_state.sequence_id,
            parent_expr.size(),
            parent_state.parent_pos
        });
        break;
      }
    }
  }
  return;
}
inline void EarleyParser::Predict(const State& state) {
  // If it's an unexpanded rule, we need to expand it,
  // and add all the possible rules into the queue.
  if (state.sequence_id == kUnexpandedRuleStartSequenceId) {
    auto cur_rule_id = state.rule_id;
    auto cur_rule_body_id = grammar_->GetRule(cur_rule_id).body_expr_id;
    auto cur_rule_body = grammar_->GetRuleExpr(cur_rule_body_id);
    if (cur_rule_body.type == RuleExprType::kTagDispatch) {
      const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
      bool in_vec = ptr != states.back().end();
      if (in_vec) {
        ptr->predictions->emplace_back(std::make_pair(cur_rule_id, cur_rule_body_id));
      } else {
        states.back().push_back(state);
        states.back().back().predictions =
            std::vector<std::pair<int32_t, int32_t>>({std::make_pair(cur_rule_id, cur_rule_body_id)}
            );
      }
      queue.emplace(State(
          cur_rule_id,
          cur_rule_body_id,
          grammar_->root_tag_dispatch_fsm->StartNode(),
          states.size() - 1
      ));
      return;
    } else {
      // XGRAMMAR_DCHECK(cur_rule_body.type == RuleExprType::kChoices);
      bool in_vec =
          std::find(states.back().begin(), states.back().end(), state) != states.back().end();
      if (!in_vec) {
        states.back().push_back(state);
        states.back().back().predictions = std::vector<std::pair<int32_t, int32_t>>();
      }
      const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
      for (auto sequence_id : cur_rule_body) {
        ptr->predictions->emplace_back(std::make_pair(cur_rule_id, sequence_id));
        queue.push(State(cur_rule_id, sequence_id, 0, states.size() - 1));
      }
      return;
    }
  }
  auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
  // Xgrammar_LOG(INFO) << "Predict: " << state << ", type is " << int(cur_rule.type) << std::endl;
  //  If the current state is the end of the rule, we do not need to predict.
  if (state.element_id == cur_rule.size()) {
    return;
  }
  switch (cur_rule.type) {
    // These types will never predict other new rules.
    case RuleExprType::kByteString:
    case RuleExprType::kCharacterClass:
    case RuleExprType::kCharacterClassStar:
    case RuleExprType::kEmptyStr: {
      // Xgrammar_LOG(INFO) << "Predict: " << state << ", type is " << int(cur_rule.type)
      //   << ", and returned." << std::endl;
      return;
    }
    // If the type if kRuleRef, then:
    // 1. Add the new rule into the queue.
    // 2. Mark the new rule as a prediction of the current state.
    case RuleExprType::kRuleRef: {
      const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
      bool in_vec = ptr != states.back().end();
      if (in_vec) {
        ptr->predictions->emplace_back(std::make_pair(cur_rule[0], kUnexpandedRuleStartSequenceId));
      } else {
        states.back().push_back(state);
        states.back().back().predictions = std::vector<std::pair<int32_t, int32_t>>();
        states.back().back().predictions->emplace_back(
            std::make_pair(cur_rule[0], kUnexpandedRuleStartSequenceId)
        );
      }
      queue.emplace(State(cur_rule[0], kUnexpandedRuleStartSequenceId, 0, states.size() - 1));
      break;
    }
    // If the type if kSequence, then:
    // 1. Add the new rule_expr into the queue.
    // 2. Mark the new rule_expr as a prediction of the current state.
    case RuleExprType::kSequence: {
      const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
      bool in_vec = ptr != states.back().end();
      if (in_vec) {
        ptr->predictions->emplace_back(std::make_pair(state.rule_id, cur_rule[state.element_id]));
      } else {
        states.back().push_back(state);
        states.back().back().predictions = std::vector<std::pair<int32_t, int32_t>>();
        states.back().back().predictions->emplace_back(
            std::make_pair(state.rule_id, cur_rule[state.element_id])
        );
      }
      queue.emplace(state.rule_id, cur_rule[state.element_id], 0, states.size() - 1);
      return;
    }
      // If the type if kSequence, then:
    // 1. Add all the new rule_exprs into the queue.
    // 2. Mark the new rule_exprs as predictions of the current state.
    case RuleExprType::kChoices: {
      const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
      bool in_vec = ptr != states.back().end();
      if (!in_vec) {
        states.back().push_back(state);
        states.back().back().predictions = std::vector<std::pair<int32_t, int32_t>>();
        for (const auto& sequence_id : cur_rule) {
          states.back().back().predictions->emplace_back(std::make_pair(state.rule_id, sequence_id)
          );
          queue.emplace(state.rule_id, sequence_id, 0, states.size() - 1);
        }
      } else {
        for (const auto& sequence_id : cur_rule) {
          ptr->predictions->emplace_back(std::make_pair(state.rule_id, sequence_id));
          queue.emplace(state.rule_id, sequence_id, 0, states.size() - 1);
        }
      }
      return;
    }
    case RuleExprType::kTagDispatch: {
      const auto& root_tag_dispatch_fsm = grammar_->root_tag_dispatch_fsm;
      if (root_tag_dispatch_fsm->IsEndNode(state.element_id)) {
        // XGRAMMAR_DCHECK(grammar_->tag_dispatch_end_node_to_rule_id.count(state.element_id))
        // << "The end node of the tag dispatch fsm does not correspond to any rule id";
        auto refered_rule_id = grammar_->tag_dispatch_end_node_to_rule_id.at(state.element_id);
        const auto& ptr = std::find(states.back().begin(), states.back().end(), state);
        bool in_vec = ptr != states.back().end();
        if (!in_vec) {
          states.back().push_back(state);
          states.back().back().predictions = std::vector<std::pair<int32_t, int32_t>>(
              {std::make_pair(refered_rule_id, kUnexpandedRuleStartSequenceId)}
          );
        } else {
          ptr->predictions->emplace_back(
              std::make_pair(refered_rule_id, kUnexpandedRuleStartSequenceId)
          );
        }
        queue.emplace(State(refered_rule_id, kUnexpandedRuleStartSequenceId, 0, states.size() - 1));
      }
      return;
    }
  }
}
inline void EarleyParser::Scan(const State& state, const uint8_t& ch) {
  // Xgrammar_LOG(INFO) << "Scan: " << state << ", ch is " << ch << std::endl;
  //  An unexpanded rule cannot be scanned.
  if (state.sequence_id == kUnexpandedRuleStartSequenceId) {
    return;
  }
  auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
  // If the current state is the end of the rule, we do not need to scan.
  if (state.element_id == cur_rule.size()) {
    return;
  }
  switch (cur_rule.type) {
    // These types can never accept a character directly.
    case (RuleExprType::kRuleRef):
    case (RuleExprType::kChoices):
    case (RuleExprType::kSequence):
    case (RuleExprType::kEmptyStr):
      return;
    // In kCharacterClass, we use a negative integer to represent how many UTF-8 codes
    // are left.
    case (RuleExprType::kCharacterClass): {
      if (IsAccepted(state, ch)) {
        // It can still match at least one UTF-8 code.
        if (state.element_id < -1) {
          queue.emplace(
              State{state.rule_id, state.sequence_id, state.element_id + 1, state.parent_pos}
          );
          return;
        }
        // The UTF-8 codes have been matched. The sequence is complete.
        if (state.element_id == -1) {
          queue.emplace(State{state.rule_id, state.sequence_id, cur_rule.size(), state.parent_pos});
          return;
        }
        // It's a brand new sequence.
        auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
        if (num_bytes > 1) {
          queue.emplace(State{state.rule_id, state.sequence_id, -(num_bytes - 1), state.parent_pos}
          );
        } else {
          queue.emplace(State{state.rule_id, state.sequence_id, cur_rule.size(), state.parent_pos});
        }
      }
      return;
    }
    case (RuleExprType::kCharacterClassStar): {
      if (state.element_id <= -1) {
        queue.emplace(
            State{state.rule_id, state.sequence_id, state.element_id + 1, state.parent_pos}
        );
      } else {
        // XGRAMMAR_DCHECK(state.element_id == 0);
        auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
        // XGRAMMAR_DCHECK(accepted);
        queue.emplace(State{state.rule_id, state.sequence_id, num_bytes - 1, state.parent_pos});
      }
      return;
    }
    case (RuleExprType::kByteString): {
      if (cur_rule.size() >= state.element_id) {
        return;
      }
      if (ch == cur_rule[state.element_id]) {
        queue.emplace(state.rule_id, state.sequence_id, state.element_id + 1, state.parent_pos);
      }
      return;
    }
    case (RuleExprType::kTagDispatch): {
      auto root_tag_dispatch_fsm = grammar_->root_tag_dispatch_fsm;
      if (!root_tag_dispatch_fsm) {
        XGRAMMAR_LOG(FATAL
        ) << "The grammar does not have a root tag dispatch rule; it is not built.";
        XGRAMMAR_UNREACHABLE();
      }
      auto start_node = root_tag_dispatch_fsm->StartNode();
      auto next_node = root_tag_dispatch_fsm->Transition(state.element_id, ch);
      auto new_stack_element = state;
      if (next_node == CompactFSM::NO_TRANSITION) {
        // Case 1. The new char cannot continue to be accepted by the tag dispatch fsm.
        // We try to accept the new char from the start node. If accepted, we go to the target node.
        // If it still cannot be accepted, we stay at the start node.
        auto new_next_node = root_tag_dispatch_fsm->Transition(start_node, ch);
        new_stack_element.element_id =
            new_next_node == CompactFSM::NO_TRANSITION ? start_node : new_next_node;
        queue.emplace(new_stack_element);
        return;
      } else {
        // Case 2. The new char can continue to be accepted by the tag dispatch fsm.
        // We need to update the element id to the next node.
        new_stack_element.element_id = next_node;
        queue.emplace(new_stack_element);
        return;
      }
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
  // Xgrammar_LOG(INFO) << "Start Scan: " << ch << ", the " << history_states.size() << "th
  // character"
  //   << std::endl;
  for (const auto& state : latest_states) {
    Scan(state, ch);
  }
  if (queue.empty()) {
    return false;
  }
  // We need a copy of the states, since we need to rollback the states.
  history_states.push_back(std::vector<State>());
  states.push_back(std::vector<State>());
  std::unordered_set<State, StateHash> visited;
  while (!queue.empty()) {
    // Xgrammar_LOG(INFO) << "Now Size: " << queue.size() << std::endl;
    const auto& state = queue.front();
    if (visited.find(queue.back()) != visited.end()) {
      queue.pop();
      continue;
    }
    visited.insert(state);
    history_states.back().push_back(state);
    Complete(state);
    Predict(state);
    queue.pop();
  }
  can_reach_end.push_back(CanReachEnd());
  return true;
}

EarleyParser::EarleyParser(const Grammar& grammar) : grammar_(grammar) {
  PushInitialState(
      State(grammar_->GetRootRuleId(), kUnexpandedRuleStartSequenceId, 0, State::kNoParent)
  );
}
inline bool EarleyParser::IsAccepted(const State& state, uint8_t ch) const {
  auto current_sequence = grammar_->GetRuleExpr(state.sequence_id);
  if (current_sequence.type == Grammar::Impl::RuleExprType::kTagDispatch) {
    return true;
  }

  if (state.parent_pos == State::kNoParent && current_sequence.size() == state.element_id) {
    // This StackElement means previous elements has matched the complete rule.
    // But we are still need to accept a new character, so this stack will become invalid.
    return false;
  }

  if (current_sequence.type == RuleExprType::kCharacterClass ||
      current_sequence.type == RuleExprType::kCharacterClassStar) {
    if (state.element_id < 0) {
      return (ch & 0xC0) == 0x80;
    }
    auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
    if (!accepted) {
      return false;
    }
    bool is_negative = static_cast<bool>(current_sequence[0]);
    if (num_bytes > 1) {
      return is_negative;
    }
    for (int i = 1; i < current_sequence.size(); i += 2) {
      if (current_sequence[i] <= ch && ch <= current_sequence[i + 1]) {
        return !is_negative;
      }
    }
    return is_negative;
  } else if (current_sequence.type == RuleExprType::kByteString) {
    return current_sequence[state.element_id] == ch;
  } else {
    XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in CheckIfAccepted: "
                        << static_cast<int>(current_sequence.type);
  }
  return false;
}

void EarleyParser::PushInitialState(const State& stack_element) {
  history_states.push_back(std::vector<State>());
  states.push_back(std::vector<State>());
  queue.push(stack_element);
  std::unordered_set<State, StateHash> visited;
  while (!queue.empty()) {
    const auto& state = queue.front();
    // Xgrammar_LOG(INFO) << "Now Size: " << queue.size() << ", " << state << std::endl;
    if (visited.find(queue.back()) != visited.end()) {
      std::cout << state << " is skipped!" << std::endl;
      queue.pop();
      continue;
    }
    visited.insert(state);
    history_states.back().push_back(state);
    // Xgrammar_LOG(INFO) << "Completion" << std::endl;
    Complete(state);
    // Xgrammar_LOG(INFO) << "Prediction" << std::endl;
    Predict(state);
    // Xgrammar_LOG(INFO) << "Prediction Done" << std::endl;
    queue.pop();
  }
  // Xgrammar_LOG(INFO) << "Loop Done" << std::endl;
  can_reach_end.push_back(CanReachEnd());
  // Xgrammar_LOG(INFO) << "INIT END" << std::endl;
  return;
}
}  // namespace xgrammar
