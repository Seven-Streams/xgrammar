/*!
 *  Copyright (c) 2025 by Contributors
 * \file xgrammar/earley_parser.h
 * \brief The header for the definition of the Earley parser.
 */
#ifndef XGRAMMAR_EARLEY_PARSER_H_
#define XGRAMMAR_EARLEY_PARSER_H_
#include <cstdint>
#include <optional>
#include <queue>
#include <vector>

#include "xgrammar/grammar.h"
namespace xgrammar {
// Be the same as StackElement.
struct State {
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
  /*! \brief Store all the possible predictions rule_ids. Used for completion.
   *   The first element is the rule_id, and  the second element is the sequence_id.*/
  std::optional<std::vector<std::pair<int32_t, int32_t>>> predictions = std::nullopt;
  /*! \brief A parent_id value of kNoParent means this StackElement is the root of the tree. */
  static constexpr int32_t kNoParent = -1;

  constexpr State() = default;
  constexpr State(const State&) = default;
  constexpr State(
      int32_t rule_id, int32_t sequence_id, int32_t element_id, int32_t parent_pos = kNoParent
  )
      : rule_id(rule_id),
        sequence_id(sequence_id),
        element_id(element_id),
        parent_pos(parent_pos) {}

  // The element is invalid when sequence_id is -1.
  bool IsInvalid() const { return sequence_id == -1; }

  bool operator==(const State& other) const {
    return rule_id == other.rule_id && sequence_id == other.sequence_id &&
           element_id == other.element_id && parent_pos == other.parent_pos;
  }
};
class StateHash {
 public:
  size_t operator()(const State& state) const {
    return std::hash<int32_t>()(state.rule_id) ^ std::hash<int32_t>()(state.sequence_id) ^
           std::hash<int32_t>()(state.element_id) ^ std::hash<int32_t>()(state.parent_pos);
  }
};
class EarleyParser {
 private:
  /*! \brief The grammar to be parsed. */
  Grammar grammar_;
  /*! \brief The tree storing all states. It's used for completation. */
  std::vector<std::vector<State>> states = {std::vector<State>()};
  /*!
      \brief The history of states. i.e. the i-th(0-base) vector
      will store the states after matching i characters. It's used
      for rollback.
   */
  std::vector<std::vector<State>> history_states = {
      std::vector<State>(),
  };

  std::queue<State> queue;
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

  EarleyParser(const Grammar& grammar);
};
}  // namespace xgrammar
#endif  // XGRAMMAR_EARLEY_PARSER_H_
