#include "earley_parser.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "grammar_data_structure.h"
#include "support/encoding.h"
#include "support/logging.h"
#include "xgrammar/grammar.h"
namespace xgrammar {

using RuleExprType = Grammar::Impl::RuleExprType;

using RuleExpr = Grammar::Impl::RuleExpr;

bool EarleyParser::IsEndOfGrammar(const State& state) const {
  // The root rule.
  if (state.parent_pos != State::kNoParent) {
    if (state.sequence_id != State::kUnexpandedRuleStartSequenceId) {
      auto seq_expr = grammar_->GetRuleExpr(state.sequence_id);
      if (seq_expr.type == RuleExprType::kTagDispatch) {
        return state.element_id != State::kTagDispatchEndFlag;
      }
    }
    return false;
  }
  if (state.sequence_id == State::kUnexpandedRuleStartSequenceId) {
    return state.element_id == State::kUnexpanedRuleFinishFlag;
  }
  auto seq_expr = grammar_->GetRuleExpr(state.sequence_id);
  if (seq_expr.type == RuleExprType::kTagDispatch) {
    return state.element_id != State::kTagDispatchEndFlag;
  } else {
    return state.element_id == seq_expr.size();
  }
}

bool EarleyParser::CanReachEnd() const {
  const auto& current_states = history_states[history_states.Size() - 1];
  return std::any_of(current_states.begin(), current_states.end(), [&](const State& state) {
    return IsEndOfGrammar(state);
  });
}

void EarleyParser::PopBackStates(int32_t cnt) {
  if (cnt >= static_cast<int32_t>(states.size())) {
    XGRAMMAR_LOG(FATAL) << "The number of states to be popped is larger than the size of states.";
  }
  states.erase(states.end() - cnt, states.end());
  for (int i = 0; i < cnt; i++) {
    history_states.PopBack();
  }
  can_reach_end.erase(can_reach_end.end() - cnt, can_reach_end.end());
  return;
}

void EarleyParser::Complete(const State& state) {
  // Check if a rule is completed.
  if (state.parent_pos == State::kNoParent) {
    return;
  }
  if (state.sequence_id == State::kUnexpandedRuleStartSequenceId) {
    if (state.element_id != State::kUnexpanedRuleFinishFlag) {
      return;
    }
  } else {
    auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
    /*
      Case 1: The current state is the end of the rule.
      Case 2: The current rule can be empty.
      Case 3: The type is kCharacterClassStar, and the element_id == 0.
      If none of the above cases is true, we return.
    */
    if (cur_rule.type == RuleExprType::kTagDispatch) {
      if (state.element_id != State::kTagDispatchEndFlag) {
        return;
      }
      queue.PushBack(State{
          state.rule_id,
          state.sequence_id,
          grammar_->root_tag_dispatch_fsm->StartNode(),
          state.parent_pos,
          0
      });
    } else {
      if (((cur_rule.size() != state.element_id)) && (cur_rule.type != RuleExprType::kEmptyStr)) {
        return;
      }
    }
  }
  // Check all the possible parent states.
  const auto& parent_states_map = states[state.parent_pos];
  if (parent_states_map.find(std::make_pair(state.rule_id, state.sequence_id)) ==
      parent_states_map.end()) {
    return;
  }
  const auto& range =
      parent_states_map.equal_range(std::make_pair(state.rule_id, state.sequence_id));
  for (auto parent_state_iter = range.first; parent_state_iter != range.second;
       parent_state_iter++) {
    const auto& parent_state = parent_state_iter->second;
    if (parent_state.sequence_id == State::kUnexpandedRuleStartSequenceId) {
      queue.PushBack(State{
          parent_state.rule_id,
          parent_state.sequence_id,
          State::kUnexpanedRuleFinishFlag,
          parent_state.parent_pos,
          0
      });
      continue;
    }
    auto parent_expr = grammar_->GetRuleExpr(parent_state.sequence_id);
    switch (parent_expr.type) {
      // These types can never predict other new rules, Thus
      // They will never be a legal parent state.
      case RuleExprType::kByteString:
      case RuleExprType::kCharacterClass:
      case RuleExprType::kCharacterClassStar:
      case RuleExprType::kEmptyStr:
      case RuleExprType::kRuleRef:
        XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in Complete: "
                            << static_cast<int>(parent_expr.type);
      // These two types can predict other new rules. We need to
      // to move to the next element.
      case RuleExprType::kSequence: {
        queue.PushBack(State{
            parent_state.rule_id,
            parent_state.sequence_id,
            parent_state.element_id + 1,
            parent_state.parent_pos,
            0
        });
        break;
      }
      // These two types can predict other new rules, and have a
      // completion is enough to complete the parent state.
      // We need to move to the end of the element. i.e. complete the parent state.
      case RuleExprType::kChoices: {
        queue.PushBack(
            {parent_state.rule_id,
             parent_state.sequence_id,
             parent_expr.size(),
             parent_state.parent_pos,
             0}
        );
        break;
      }
      case RuleExprType::kTagDispatch: {
        queue.PushBack(State{
            parent_state.rule_id,
            parent_state.sequence_id,
            State::kTagDispatchEndFlag,
            parent_state.parent_pos,
            0
        });
        queue.PushBack(
            {parent_state.rule_id,
             parent_state.sequence_id,
             grammar_->root_tag_dispatch_fsm->StartNode(),
             parent_state.parent_pos,
             0}
        );
        break;
      }
    }
  }
  return;
}

void EarleyParser::Predict(const State& state) {
  // If it's an unexpanded rule, we need to expand it,
  // and add all the possible rules into the queue.
  if (state.sequence_id == State::kUnexpandedRuleStartSequenceId) {
    // The rule is already expanded, and finished.
    if (state.element_id == State::kUnexpanedRuleFinishFlag) {
      return;
    }
    auto cur_rule_id = state.rule_id;
    auto cur_rule_body_id = grammar_->GetRule(cur_rule_id).body_expr_id;
    auto cur_rule_body = grammar_->GetRuleExpr(cur_rule_body_id);
    // There are two types of an unexpanded rule:
    // 1. The rule is a tag dispatch rule.
    // 2. The rule is a choice, consisting of multiple sequences.
    if (cur_rule_body.type == RuleExprType::kTagDispatch) {
      auto& states_map = states.back();
      states_map.insert({std::make_pair(cur_rule_id, cur_rule_body_id), {state}});
      queue.PushBack(State{
          cur_rule_id,
          cur_rule_body_id,
          grammar_->root_tag_dispatch_fsm->StartNode(),
          int32_t(states.size()) - 1,
          0
      });
      return;
    }
    XGRAMMAR_DCHECK(cur_rule_body.type == RuleExprType::kChoices);
    auto& states_map = states.back();
    std::vector<std::pair<int32_t, int32_t>> prediction;
    for (auto sequence_id : cur_rule_body) {
      queue.PushBack(State{cur_rule_id, sequence_id, 0, int32_t(states.size()) - 1, 0});
      if (std::find_if(
              prediction.begin(),
              prediction.end(),
              [&](const std::pair<int32_t, int32_t> p) {
                return p == std::make_pair(cur_rule_id, sequence_id);
              }
          ) != prediction.end()) {
        continue;
      }
      prediction.emplace_back(cur_rule_id, sequence_id);
      states_map.insert({std::make_pair(cur_rule_id, sequence_id), state});
    }
    return;
  }
  auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
  //  If the current state is the end of the rule, we do not need to predict.
  if (cur_rule.type == RuleExprType::kTagDispatch) {
    if (state.element_id == State::kTagDispatchEndFlag) {
      return;
    }
  } else {
    // If the current state is the end of the rule, we do not need to predict,
    // since the rule is already completed.
    if (state.element_id == cur_rule.size()) {
      return;
    }
  }
  switch (cur_rule.type) {
    // These types will never be added into the queue.
    case RuleExprType::kByteString:
    case RuleExprType::kCharacterClass:
    case RuleExprType::kCharacterClassStar:
    case RuleExprType::kRuleRef: {
      XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in Predict: "
                          << static_cast<int>(cur_rule.type);
    }
    // If the type is kChoices, then it shoule have been expanded.
    case RuleExprType::kChoices: {
      return;
    }
    // If the type is kEmptyStr, then it should be a empty string,
    // it will predict nothing.
    case RuleExprType::kEmptyStr: {
      return;
    }
    /* If the type is kSequence, then it should be a sequence consisting of:
      - kByteString
      - kCharacterClass
      - kCharacterClassStar
      - RuleRef
      RuleRef need to be predicted, and the others only need to be completed or scanned.
    */
    case RuleExprType::kSequence: {
      const auto& element_expr = grammar_->GetRuleExpr(cur_rule[state.element_id]);
      if (element_expr.type == RuleExprType::kRuleRef) {
        auto& states_map = states.back();
        const auto& rule_ref_id = element_expr[0];
        const auto& ref_rule = grammar_->GetRule(rule_ref_id);
        const auto& ref_rule_body_id = ref_rule.body_expr_id;
        const auto& ref_rule_body = grammar_->GetRuleExpr(ref_rule_body_id);
        if (ref_rule_body.type == RuleExprType::kTagDispatch) {
          states_map.insert({std::make_pair(rule_ref_id, ref_rule_body_id), state});
          queue.PushBack(State{
              rule_ref_id,
              ref_rule_body_id,
              grammar_->root_tag_dispatch_fsm->StartNode(),
              int32_t(states.size()) - 1,
              0
          });
        } else {
          XGRAMMAR_DCHECK(ref_rule_body.type == RuleExprType::kChoices);
          for (const auto& sequence_id : ref_rule_body) {
            states_map.insert({std::make_pair(rule_ref_id, sequence_id), state});
            queue.PushBack(State{rule_ref_id, sequence_id, 0, int32_t(states.size()) - 1, 0});
          }
        }
      }
      if (element_expr.type == RuleExprType::kCharacterClassStar && state.sub_element_id == 0) {
        queue.PushBack(
            State{state.rule_id, state.sequence_id, state.element_id + 1, state.parent_pos, 0}
        );
        return;
      }
      return;
    }
    case RuleExprType::kTagDispatch: {
      if (state.element_id == State::kTagDispatchEndFlag) {
        // The tag has been dispatched.
        return;
      }
      const auto& root_tag_dispatch_fsm = grammar_->root_tag_dispatch_fsm;
      if (root_tag_dispatch_fsm->IsEndNode(state.element_id)) {
        XGRAMMAR_DCHECK(grammar_->tag_dispatch_end_node_to_rule_id.count(state.element_id))
            << "The end node of the tag dispatch fsm does not correspond to any rule id";
        auto refered_rule_id = grammar_->tag_dispatch_end_node_to_rule_id.at(state.element_id);
        auto& states_map = states.back();
        states_map.insert(
            {std::make_pair(refered_rule_id, State::kUnexpandedRuleStartSequenceId), state}
        );
        queue.PushBack(
            State(refered_rule_id, State::kUnexpandedRuleStartSequenceId, 0, states.size() - 1, 0)
        );
      }
      return;
    }
  }
}

void EarleyParser::Scan(const State& state, const uint8_t& ch) {
  //  An unexpanded rule cannot be scanned.
  if (state.sequence_id == State::kUnexpandedRuleStartSequenceId) {
    return;
  }
  auto cur_rule = grammar_->GetRuleExpr(state.sequence_id);
  // If the current state is the end of the rule, we do not need to scan.
  if (state.element_id == cur_rule.size() && cur_rule.type != RuleExprType::kTagDispatch) {
    return;
  }
  switch (cur_rule.type) {
    // These types can never accept a character directly.
    case (RuleExprType::kChoices):
    case (RuleExprType::kEmptyStr):
      return;
    // These types can never be added into the queue.
    case (RuleExprType::kRuleRef):
    case (RuleExprType::kByteString):
    case (RuleExprType::kCharacterClass):
    case (RuleExprType::kCharacterClassStar): {
      XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType in Scan: " << static_cast<int>(cur_rule.type);
    }
    case (RuleExprType::kSequence): {
      const auto& element_expr = grammar_->GetRuleExpr(cur_rule[state.element_id]);
      // The element is a rule reference, we do not need to scan it.
      switch (element_expr.type) {
        case (RuleExprType::kChoices):
        case (RuleExprType::kSequence):
        case (RuleExprType::kEmptyStr):
        case (RuleExprType::kTagDispatch): {
          XGRAMMAR_LOG(FATAL) << "Unexpected RuleExprType of sequence's element: "
                              << static_cast<int>(element_expr.type);
        }
        case (RuleExprType::kRuleRef):
          return;
        case (RuleExprType::kByteString): {
          // The rule has been completed.
          if (state.sub_element_id == element_expr.size()) {
            return;
          }
          if (element_expr[state.sub_element_id] == ch) {
            auto new_state = state;
            new_state.sub_element_id++;
            if (new_state.sub_element_id == element_expr.size()) {
              new_state.element_id++;
              new_state.sub_element_id = 0;
            }
            queue.PushBack(new_state);
          }
          return;
        }
        case (RuleExprType::kCharacterClass):
        case (RuleExprType::kCharacterClassStar): {
          if (!IsAccepted(state, ch)) {
            return;
          }
          if (state.sub_element_id > 0) {
            auto new_state = state;
            new_state.sub_element_id--;
            if (new_state.sub_element_id == 0) {
              if (element_expr.type == RuleExprType::kCharacterClassStar) {
                queue.PushBack(new_state);
              }
              new_state.element_id++;
              new_state.sub_element_id = 0;
            }
            queue.PushBack(new_state);
            return;
          }
          auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
          if (!accepted) {
            return;
          }
          // A single byte is accepted.
          if (num_bytes == 1) {
            auto new_state = state;
            new_state.element_id++;
            queue.PushBack(new_state);
            if (element_expr.type == RuleExprType::kCharacterClassStar) {
              // The element is a character class star, we need to add the old state again.
              queue.PushBack(state);
            }
            return;
          }
          // A UTF8 character is accepted.
          auto new_state = state;
          new_state.sub_element_id = num_bytes - 1;
          queue.PushBack(new_state);
          return;
        }
      }
      return;
    }
    case (RuleExprType::kTagDispatch): {
      if (state.element_id == State::kTagDispatchEndFlag) {
        // The tag has been dispatched.
        return;
      }
      auto root_tag_dispatch_fsm = grammar_->root_tag_dispatch_fsm;
      if (!root_tag_dispatch_fsm) {
        XGRAMMAR_LOG(FATAL
        ) << "The grammar does not have a root tag dispatch rule; it is not built.";
        XGRAMMAR_UNREACHABLE();
      }
      if (root_tag_dispatch_fsm->IsEndNode(state.element_id)) {
        // The tag has been dispatched.
        return;
      }
      auto start_node = root_tag_dispatch_fsm->StartNode();
      auto next_node = root_tag_dispatch_fsm->Transition(state.element_id, ch);
      auto new_state = state;
      if (next_node == CompactFSM::NO_TRANSITION) {
        // Case 1. The new char cannot continue to be accepted by the tag dispatch fsm.
        // We try to accept the new char from the start node. If accepted, we go to the target
        // node. If it still cannot be accepted, we stay at the start node.
        auto new_next_node = root_tag_dispatch_fsm->Transition(start_node, ch);
        new_state.element_id =
            new_next_node == CompactFSM::NO_TRANSITION ? start_node : new_next_node;
        queue.PushBack(new_state);
        return;
      } else {
        // Case 2. The new char can continue to be accepted by the tag dispatch fsm.
        // We need to update the element id to the next node.
        new_state.element_id = next_node;
        queue.PushBack(new_state);
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
  const auto& latest_states = history_states[history_states.Size() - 1];
  for (const auto& state : latest_states) {
    Scan(state, ch);
  }
  if (queue.begin() == queue.end()) {
    return false;
  }
  tmp_states.clear();
  states.emplace_back();
  std::vector<State> visited;

  while (queue.begin() != queue.end()) {
    const auto& state_iter = queue.begin();
    const auto state = *state_iter;
    if (std::find_if(visited.begin(), visited.end(), [&](const State& s) {
          return CheckingStateEqual()(state, s);
        }) != visited.end()) {
      queue.Erase(state_iter);
      continue;
    }
    visited.push_back(state);
    tmp_states.push_back(state);
    Complete(state);
    Predict(state);
    queue.Erase(state_iter);
  }
  history_states.Insert(tmp_states);
  can_reach_end.push_back(CanReachEnd());
  return true;
}

EarleyParser::EarleyParser(const Grammar& grammar, const State& init_state, const bool& need_expand)
    : grammar_(grammar) {
  if (init_state.IsInvalid()) {
    this->init_state = State(
        grammar_->GetRootRuleId(), State::kUnexpandedRuleStartSequenceId, 0, State::kNoParent, 0
    );
  } else {
    this->init_state = init_state;
  }
  if (need_expand) {
    PushInitialState(this->init_state);
    return;
  }
  states.emplace_back();
  history_states.Insert({this->init_state});
  can_reach_end.push_back(CanReachEnd());
}

bool EarleyParser::IsAccepted(const State& state, uint8_t ch) const {
  auto sequence_expr = grammar_->GetRuleExpr(state.sequence_id);
  auto element_expr = grammar_->GetRuleExpr(sequence_expr[state.element_id]);
  XGRAMMAR_DCHECK(
      element_expr.type == RuleExprType::kCharacterClass ||
      element_expr.type == RuleExprType::kCharacterClassStar
  ) << "The element type is not supported!";
  if (state.sub_element_id > 0) {
    return (ch & 0xC0) == 0x80;
  }
  auto [accepted, num_bytes, codepoint] = HandleUTF8FirstByte(ch);
  if (!accepted) {
    return false;
  }
  bool is_negative = static_cast<bool>(element_expr[0]);
  if (num_bytes > 1) {
    return is_negative;
  }
  for (int i = 1; i < element_expr.size(); i += 2) {
    if (element_expr[i] <= ch && ch <= element_expr[i + 1]) {
      return !is_negative;
    }
  }
  return is_negative;
}

void EarleyParser::PushInitialState(const State& state) {
  tmp_states.clear();
  states.emplace_back();
  if (state.IsInvalid()) {
    queue.PushBack(State(
        grammar_->GetRootRuleId(), State::kUnexpandedRuleStartSequenceId, 0, State::kNoParent, 0
    ));
  } else {
    queue.PushBack(state);
  }
  std::vector<State> visited;
  while (queue.begin() != queue.end()) {
    const auto& state_iter = queue.begin();
    const auto state = *state_iter;
    if (std::find_if(visited.begin(), visited.end(), [&](const State& s) {
          return CheckingStateEqual()(state, s);
        }) != visited.end()) {
      queue.Erase(state_iter);
      continue;
    }
    visited.push_back(state);
    tmp_states.push_back(state);
    Complete(state);
    Predict(state);
    queue.Erase(state_iter);
  }
  history_states.Insert(tmp_states);
  can_reach_end.push_back(CanReachEnd());
  return;
}

void EarleyParser::ParserReset() {
  states.clear();
  history_states = CSRArray<State>();
  queue.Clear();
  can_reach_end.clear();
  PushInitialState(State(
      grammar_->GetRootRuleId(), State::kUnexpandedRuleStartSequenceId, 0, State::kNoParent, 0
  ));
}

void EarleyParser::PopFrontStates(const int32_t& cnt) {
  if (cnt >= history_states.Size()) {
    return;
  }
  // TODO:
  return;
}

}  // namespace xgrammar
