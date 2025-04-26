/*!
 *  Copyright (c) 2025 by Contributors
 * \file xgrammar/earley_parser.h
 * \brief The header for the definition of the Earley parser.
 */
#ifndef XGRAMMAR_EARLEY_PARSER_H_
#define XGRAMMAR_EARLEY_PARSER_H_
#include <cstdint>
#include <functional>
#include <list>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgrammar/grammar.h"
namespace xgrammar {
// Be the same as StackElement.
struct State {
  /*! \brief A sequence_id value of kUnexpandedRuleStartSequenceId means a rule hasn't been
   * expanded.*/
  static constexpr int32_t kUnexpandedRuleStartSequenceId = 128000;
  /*! \brief A element_id value of kUnexpanedRuleFinishFlag means an unexpanded rule is ended.*/
  static constexpr int32_t kUnexpanedRuleFinishFlag = 128001;
  /*! \brief A element_id value of kTagDispatchEndFlag means the kTagDispatch pair is finished.*/
  static constexpr int32_t kTagDispatchEndFlag = -2;
  /*! \brief A parent_id value of kNoParent means this State is the root of the tree. */
  static constexpr int32_t kNoParent = -1;
  /*! \brief The rule's id. Used for debug purposes. */
  int32_t rule_id = -1;
  /*! \brief Which choice in this rule is selected. */
  int32_t sequence_id = -1;
  /*! \brief Which element of the choice sequence is to be visited. When the current sequence is
   * a tag dispatch rule, this element id the currently visited node. */
  int32_t element_id = -1;
  /*! \brief The id of the parent node in the Earley parser. i.e. from the k-th character, the
   * rule starts to match the string.*/
  int32_t parent_pos = -1;
  /*! \brief The id of the sub element in the current selement of the sequence. */
  int32_t sub_element_id = 0;

  constexpr State() = default;
  constexpr State(const State&) = default;
  State& operator=(const State&) = default;
  constexpr State(
      int32_t rule_id,
      int32_t sequence_id,
      int32_t element_id,
      int32_t parent_pos,
      int32_t sub_element_id
  )
      : rule_id(rule_id),
        sequence_id(sequence_id),
        element_id(element_id),
        parent_pos(parent_pos),
        sub_element_id(sub_element_id) {}

  // The element is invalid when sequence_id is -1.
  bool IsInvalid() const { return sequence_id == -1; }

  bool operator==(const State& other) const {
    return rule_id == other.rule_id && sequence_id == other.sequence_id &&
           element_id == other.element_id && sub_element_id == other.sub_element_id;
  }
  friend std::ostream& operator<<(std::ostream& os, const State& state) {
    os << "State(rule_id=" << state.rule_id << ", sequence_id=" << state.sequence_id
       << ", element_id=" << state.element_id << ", parent_pos=" << state.parent_pos
       << ", sub_element_id=" << state.sub_element_id << ")";
    return os;
  }
};
class StateHash {
 public:
  size_t operator()(const State& state) const {
    return std::hash<int32_t>()(state.rule_id) << 16 ^
           std::hash<int32_t>()(state.sequence_id) << 8 ^
           std::hash<int32_t>()(state.element_id) << 4 ^ std::hash<int32_t>()(state.sub_element_id);
  }
};

class CheckingStateEqual {
 public:
  bool operator()(const State& lhs, const State& rhs) const {
    return lhs.rule_id == rhs.rule_id && lhs.sequence_id == rhs.sequence_id &&
           lhs.element_id == rhs.element_id && lhs.parent_pos == rhs.parent_pos &&
           lhs.sub_element_id == rhs.sub_element_id;
  }
};

class CheckingStateHash {
 public:
  size_t operator()(const State& state) const {
    return std::hash<int32_t>()(state.rule_id) << 16 ^
           std::hash<int32_t>()(state.sequence_id) << 8 ^
           std::hash<int32_t>()(state.element_id) << 4 ^
           std::hash<int32_t>()(state.parent_pos) << 2 ^ std::hash<int32_t>()(state.sub_element_id);
  }
};
class EarleyParser {
 protected:
  /*! \brief The grammar to be parsed. */
  Grammar grammar_;
  /*! \brief The tree storing all states. It's used for completation. */
  std::vector<std::unordered_map<
      std::pair<int32_t, int32_t>,
      std::unordered_set<State, CheckingStateHash, CheckingStateEqual>>>
      states;
  /*!
      \brief The history of states. i.e. the i-th(0-base) vector
      will store the states after matching i characters. It's used
      for rollback.
   */
  std::vector<std::vector<State>> history_states;

  /*! \brief The vector stores whether at present, the grammar can reach the end. */
  std::vector<bool> can_reach_end;

  /*! \brief It's the processing queue of the earley parser.*/
  std::list<State> queue;

  State init_state;
  /*!
    \brief The scanning operation of the Earley parser.
  */
  void Scan(const State& state, const uint8_t& ch);

  /*!
      \brief The prediction operation of the Earley parser.
  */
  void Predict(const State& state);

  /*!
      \brief The completion operation of the Earley parser.
  */
  void Complete(const State& state);

  /*!
    \brief Check if the state is the end of the grammar.
    \param state The state to be checked.
    \return True if the state is the end of the grammar, false otherwise.
  */
  bool IsEndOfGrammar(const State& state) const;

  /*!
    \brief Check if a character can be accepted.
  */
  bool IsAccepted(const State& state, uint8_t ch) const;

 public:
  /*!
    \brief From the current states, advance to the next state.
    \param ch The character to be advanced.
    \return True if the character is accepted, false otherwise.
    \note If the character isn't accepted, then the states won't be changed.
  */
  bool Advance(const uint8_t& ch);

  /*!
    \brief Remove the newly added states.
    \param count The number of states to be removed.
  */
  void PopBackStates(int32_t count);

  /*!
    \brief Check if the root rule is completed.
    \return True if the root rule is completed, false otherwise.
  */
  bool CanReachEnd() const;

  /*!
    \brief Push the initial state into the Earley parser.
    \param state The initial state to be pushed.
  */
  void PushInitialState(const State& state);

  /*!
    \brief Constructor of the Earley parser.
    \param grammar The grammar to be parsed.
    \param initial_state The initial state to be pushed into the parser.
  */
  EarleyParser(const Grammar& grammar, const State& initial_state);

  /*!
    \brief Reset the parser.
    \note This function is used to reset the parser, and initialize the
    parser with the root rule.
  */
  void ParserReset();
};
}  // namespace xgrammar
#endif  // XGRAMMAR_EARLEY_PARSER_H_
